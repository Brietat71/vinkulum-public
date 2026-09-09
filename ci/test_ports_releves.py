"""Contre-épreuves des ports relevés et de leur assemblage énergétique.

Les oracles denses restent limités aux petites matrices. La grande chaîne
utilise sa solution analytique. Les tolérances de comparaison flottante
ne constituent pas des certificats machine des enveloppes conditionnelles.
"""
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import eigh
from scipy.sparse import csc_matrix, csr_matrix, diags

import condensation_energie
from condensation_energie import BudgetQR, CondensationEnergie
from condensation_energie_lu import CondensationEnergieLU
from energie_ports import energie_console
from modeles_ports import (borne_lambda_trace_console, console, metrique_ports,
                           reference_chaine)
from ports_krylov import InterieurKrylov
from ports_releves import AssemblagePorts, PortsReleves
from reference_ports_precision import schur_console, schur_facteur


def _petit(seed=81, ni=8, p=2):
    """Énergie rectangulaire et masse couplée, toutes deux coercives."""
    rng = np.random.default_rng(seed)
    n = ni+p
    d = np.vstack((2*np.eye(n), .35*rng.standard_normal((n+3, n))))
    z = .25*rng.standard_normal((n+2, n))
    m = z.T@z+.5*np.eye(n)
    k = d.T@d
    i, s = np.arange(ni), np.arange(ni, n)
    metric = k[ni:, ni:]-k[ni:, :ni]@np.linalg.solve(k[:ni, :ni], k[:ni, ni:])
    minimum = eigh(k[:ni, :ni], m[:ni, :ni], eigvals_only=True)[0]
    return d, m, i, s, metric, .7*minimum, .55*np.sqrt(minimum)


def _releve(data, **options):
    d, m, i, s, metric, lower, omega = data
    return PortsReleves(d, m, i, s, metric, lower, omega, **options)


def _schur(k, m, i, s, w, omega):
    """Élimination dense indépendante, utilisée seulement sur petits cas."""
    a = k-omega**2*m
    ai, b = a[np.ix_(i, i)], a[np.ix_(i, s)]
    return w.T@(a[np.ix_(s, s)]-b.T@np.linalg.solve(ai, b))@w


def _matrice_qr(qr):
    """Modèle exact des facteurs stockés, pour l'oracle de petite taille."""
    ni, p = len(qr.i), len(qr.s)
    a = qr.r.toarray()/qr.echelles[None, :]
    b = np.linalg.solve(qr.w.T, qr.p.T).T
    e = np.linalg.solve(qr.w.T, qr.e.T).T
    facteur = np.block([[a, b], [np.zeros((len(e), ni)), e]])
    ordre = np.r_[qr.i, qr.s]
    k = np.empty((ni+p, ni+p))
    k[np.ix_(ordre, ordre)] = facteur.T@facteur
    return k


def _chaine(n, k=1., m=1.):
    d = diags((np.full(n, np.sqrt(k)), np.full(n-1, -np.sqrt(k))),
              (0, -1), shape=(n, n), format="csr")
    masse = diags(np.full(n, m), format="csc")
    minimum = 4*k/m*np.sin(np.pi/(2*n))**2
    return (d, masse, np.arange(n-1), np.array([n-1]), np.array([[k/n]]),
            .9*minimum, .8*np.sqrt(minimum))


class CondensationEnergieTests(unittest.TestCase):
    def test_qr_rectangulaire_conserve_energie_et_dualite(self):
        d, _, i, s, metric, _, _ = _petit()
        qr = CondensationEnergie(d, i, s, metric)
        k = d.T@d
        kii = k[np.ix_(i, i)]
        attendu = -np.linalg.solve(kii, k[np.ix_(i, s)]@qr.w)
        assert_allclose(qr.psi, attendu, atol=2e-14, rtol=2e-13)
        assert_allclose(_matrice_qr(qr), k, atol=3e-14, rtol=3e-13)
        assert_allclose(qr.k0, _schur(k, np.zeros_like(k), i, s, qr.w, 0.),
                        atol=3e-14, rtol=3e-13)
        rhs = np.random.default_rng(118).standard_normal((len(i), 3))
        assert_allclose(qr.solve(rhs), np.linalg.solve(kii, rhs), atol=2e-14, rtol=3e-13)
        dual = qr.dual(rhs)
        assert_allclose(dual.T@dual, rhs.T@np.linalg.solve(kii, rhs),
                        atol=3e-14, rtol=3e-13)
        x = np.random.default_rng(119).standard_normal((len(i), 2))
        assert_allclose(qr.produit(x), kii@x, atol=5e-14, rtol=3e-13)
        assert_allclose(qr.gram(x, x), (d[:, i]@x).T@(d[:, i]@x),
                        atol=8e-14, rtol=3e-13)

    def test_defaut_du_solve_statique_reste_dans_energie_et_couplage(self):
        d, _, i, s, metric, _, _ = _petit(82)
        solve_initial = condensation_energie.spsolve_triangular

        def solve_perturbe(a, b, **kwargs):
            x = solve_initial(a, b, **kwargs)
            # Injection d'une erreur de résolution ; les oracles utilisent D.
            return x+1e-3*np.arange(1, x.size+1).reshape(x.shape)

        with patch.object(condensation_energie, "spsolve_triangular", solve_perturbe):
            qr = CondensationEnergie(d, i, s, metric)
        champ = np.zeros((d.shape[1], len(s)))
        champ[i], champ[s] = qr.psi, qr.w
        deformation = d@champ
        self.assertGreater(np.linalg.norm(qr.h), 1e-4)
        self.assertGreater(np.linalg.norm(qr.k0-qr.e.T@qr.e), 1e-8)
        assert_allclose(qr.k0, deformation.T@deformation, atol=4e-13, rtol=4e-13)
        assert_allclose(qr.r0, d[:, i].T@deformation, atol=4e-13, rtol=4e-13)

    def test_energie_factorisee_ne_certifie_pas_la_matrice_assemblee(self):
        # Le produit D.T D perd ici le terme 1 de son premier coefficient.
        # La vraie énergie de D au port normalisé demeure d'ordre un.
        d = np.array([[1., 0.], [1e8, 1.]])
        qr = CondensationEnergie(d, np.array([0]), np.array([1]), np.array([[1e-16]]))
        k_assemble = d.T@d
        soustraction = _schur(k_assemble, np.zeros((2, 2)),
                             np.array([0]), np.array([1]), qr.w, 0.)
        self.assertAlmostEqual(qr.k0[0, 0], 1., delta=3e-14)
        self.assertGreater(abs(qr.k0[0, 0]-soustraction[0, 0]), .5)
        audit = qr.audit()
        self.assertIs(audit["certification_machine"], False)
        self.assertTrue(np.isfinite(audit["ecart_energie_statique"]))
        self.assertLess(audit["ecart_energie_statique"], 1e-12)

    def test_budget_qr_refuse_sans_supprimer_des_coefficients(self):
        d, _, i, s, metric, _, _ = _chaine(20)
        with self.assertRaises(BudgetQR):
            CondensationEnergie(d, i, s, metric, budget=1)

    def test_semences_ecartees_restent_dans_le_residu_du_backend(self):
        d, m, i, s, metric, lower, maximum = _petit(84)
        qr = CondensationEnergie(d, i, s, metric)
        b = .01*np.random.default_rng(85).standard_normal((len(i), 2))
        interieur = InterieurKrylov(qr.matrice(), m[np.ix_(i, i)], b,
                                    lower, maximum, tolerance=1e-20,
                                    facteur_energie=qr, seuil_semence=1e3)
        self.assertEqual(interieur.taille_interieure, 0)
        self.assertNotEqual(interieur.statut, "tolerance_estimee")
        self.assertGreater(interieur.borne_uniforme, 1e-8)
        k = _matrice_qr(qr)[np.ix_(i, i)]
        for omega in (0., maximum):
            exact = b.T@np.linalg.solve(k-omega**2*m[np.ix_(i, i)], b)
            self.assertLessEqual(np.linalg.norm(exact, 2), interieur.borne_uniforme+1e-14)
            assert_allclose(interieur.transfert(omega), np.zeros((2, 2)), atol=0., rtol=0.)


class PortsRelevesTests(unittest.TestCase):
    def test_lu_energie_forces_gram_et_refus_du_gram_singulier(self):
        data = _petit(86)
        d, _, i, _, _, _, _ = data
        r = _releve(data, methode="lu_energie", tolerance=1e-12, max_directions=len(i))
        self.assertEqual(r.methode, "lu_energie")
        facteur = r.qr
        rhs = np.random.default_rng(186).standard_normal((len(i), 3))
        k = d[:, i].T@d[:, i]
        x = facteur.solve(rhs)
        assert_allclose(x, np.linalg.solve(k, rhs), atol=5e-14, rtol=5e-13)
        assert_allclose(d[:, i].T@(d[:, i]@x), rhs, atol=5e-13, rtol=5e-13)
        assert_allclose(facteur.gram(x, x), (d[:, i]@x).T@(d[:, i]@x),
                        atol=5e-14, rtol=5e-13)
        dual = facteur.dual(rhs)
        assert_allclose(dual.T@dual, rhs.T@np.linalg.solve(k, rhs),
                        atol=5e-14, rtol=5e-13)
        audit = facteur.audit()
        self.assertIs(audit["certification_machine"], False)
        self.assertEqual(audit["norme_duale"], "approchee_par_D_I_solve")
        self.assertEqual(audit["statut_lu"], "preconditionneur_admis_sous_audit")
        # D reste de rang plein, mais fl(D_I.T D_I) perd epsilon².
        # Le préconditionneur LU doit refuser ; aucune bascule QR implicite.
        depsilon = np.array([[1., 1., 0.], [0., 1e-9, 0.], [0., 0., 1.]])
        interieur, interface = np.array([0, 1]), np.array([2])
        qr = CondensationEnergie(depsilon, interieur, interface, np.eye(1))
        assert_allclose(qr.solve(np.ones((2, 1))), [[1.], [0.]], atol=1e-14, rtol=1e-14)
        assert_allclose(qr.k0, [[1.]], atol=1e-14, rtol=1e-14)
        with self.assertRaises((ValueError, ArithmeticError)):
            CondensationEnergieLU(depsilon, interieur, interface, np.eye(1))

    def test_lu_console_128_cible_decimal_70_et_audit_original(self):
        model = console(128)
        energie = energie_console(model.metadata)
        lower = borne_lambda_trace_console(model.metadata)
        maximum = .8*np.sqrt(lower)
        r = PortsReleves(energie.d, energie.m, energie.interieur, energie.interface,
                         metrique_ports(model), lower, maximum, tolerance=1e-10,
                         max_blocs=12, max_directions=128, methode="lu_energie")
        self.assertEqual(r.statut, "tolerance_estimee")
        self.assertLessEqual(r.borne_uniforme, 1e-10)
        self.assertEqual(r.borne_uniforme, r.audit_bande["majorant_total"])
        self.assertIs(r.audit_bande["certification_machine"], False)
        self.assertIs(r.audit_modele()["certification_machine"], False)
        for omega in (0., maximum):
            with self.subTest(omega=omega):
                # Les deux oracles ont des données d'entrée distinctes :
                # paramètres physiques, puis coefficients binaires exacts de D/M.
                physiques = schur_console(model.metadata, omega, dps=70)
                donnees = schur_facteur(energie.d, energie.m, energie.interieur,
                                       energie.interface, omega, dps=70)
                for oracle in (physiques, donnees):
                    normalise = r.qr.w.T@oracle@r.qr.w
                    self.assertLessEqual(np.linalg.norm(r.schur(omega)-normalise, 2), 1e-10)

    def test_masse_couplee_origine_et_reconstruction(self):
        data = _petit(87, ni=7)
        d, m, i, s, _, _, omega_max = data
        r = _releve(data, tolerance=1e-12, max_blocs=10, max_directions=len(i))
        self.assertEqual(r.borne_uniforme, r.audit_bande["majorant_total"])
        self.assertIs(r.audit_bande["certification_machine"], False)
        self.assertGreater(np.linalg.norm(m[np.ix_(i, s)]), .1)
        k = d.T@d
        sans_couplage = m.copy()
        sans_couplage[np.ix_(i, s)] = 0.
        sans_couplage[np.ix_(s, i)] = 0.
        ecart_ignore = 0.
        for omega in np.linspace(0., omega_max, 9):
            exact = _schur(k, m, i, s, r.qr.w, omega)
            assert_allclose(r.schur(omega), exact, atol=2e-11, rtol=2e-11)
            champ = r.reconstruire(omega, np.eye(len(s)))
            a = k-omega**2*m
            uexact = -np.linalg.solve(a[np.ix_(i, i)], a[np.ix_(i, s)]@r.qr.w)
            assert_allclose(champ[i], uexact, atol=2e-9, rtol=2e-9)
            assert_allclose(champ[s], r.qr.w, atol=0., rtol=0.)
            ecart_ignore = max(ecart_ignore, np.linalg.norm(
                exact-_schur(k, sans_couplage, i, s, r.qr.w, omega), 2))
        self.assertGreater(ecart_ignore, 1e-4)

    def test_correcteur_energetique_et_residu_non_nul(self):
        data = _petit(89, ni=10)
        _, m, i, s, _, _, omega_max = data
        r = _releve(data, tolerance=1e-20, max_blocs=1, max_directions=1)
        k = _matrice_qr(r.qr)
        erreurs = []
        for omega in np.linspace(0., omega_max, 7):
            a = k-omega**2*m
            champ = r.reconstruire(omega, np.eye(len(s)))
            residu = a[np.ix_(i, np.arange(len(a)))]@champ
            erreur_energie = residu.T@np.linalg.solve(a[np.ix_(i, i)], residu)
            exact = _schur(k, m, i, s, r.qr.w, omega)
            calcule = r.schur(omega)
            assert_allclose(calcule, champ.T@a@champ, atol=8e-13, rtol=8e-13)
            assert_allclose(calcule-exact, erreur_energie, atol=8e-13, rtol=2e-10)
            self.assertGreaterEqual(np.linalg.eigvalsh(calcule-exact)[0], -1e-12)
            self.assertLessEqual(np.linalg.norm(erreur_energie, 2), r.borne_uniforme+1e-12)
            p = np.vstack((np.eye(len(s)), -omega**2/r.lambda_min*np.eye(len(s))))
            borne = p.T@r.interieur.borne_ponctuelle(omega)@p
            self.assertGreaterEqual(np.linalg.eigvalsh(borne-erreur_energie)[0], -1e-12)
            erreurs.append(np.linalg.norm(erreur_energie, 2))
        self.assertGreater(max(erreurs), 1e-7)

    def test_audit_original_compte_le_defaut_du_fonctionnel_rapporte(self):
        data = _petit(90, ni=8)
        d, m, i, s, _, _, maximum = data
        r = _releve(data, tolerance=1e-12, max_directions=len(i))
        # Injection d'un défaut du fonctionnel, sans changer le champ candidat.
        # Un audit du seul résidu ne pourrait pas détecter ce déplacement.
        r.qr.k0 = r.qr.k0+2e-5*np.eye(len(s))
        r.interieur.g = r.interieur.g+1e-3*np.eye(r.taille_interieure)
        audit = r.audit_uniforme()
        self.assertIs(audit["certification_machine"], False)
        self.assertGreater(audit["ecart_modele"], 1e-5)
        erreur_max = 0.
        for omega in np.linspace(0., maximum, 9):
            exact = _schur(d.T@d, m, i, s, r.qr.w, omega)
            erreur = np.linalg.norm(r.schur(omega)-exact, 2)
            self.assertLessEqual(erreur, audit["majorant_total"]+2e-12)
            ponctuel = r.audit_energie(omega)
            self.assertLessEqual(erreur, ponctuel["majorant_total"]+2e-12)
            erreur_max = max(erreur_max, erreur)
        self.assertGreater(erreur_max, 1e-5)

    def test_permutations_physiques_et_lignes_energetiques(self):
        data = _petit(91)
        d, m, i, s, metric, lower, maximum = data
        initial = _releve(data, tolerance=1e-12, max_directions=len(i))
        rng = np.random.default_rng(92)
        permutation = rng.permutation(d.shape[1])
        inverse = np.argsort(permutation)
        lignes = rng.permutation(d.shape[0])
        signes = rng.choice((-1., 1.), d.shape[0])
        dp = signes[:, None]*d[lignes][:, permutation]
        mp = m[np.ix_(permutation, permutation)]
        permute = PortsReleves(dp, mp, inverse[i[::-1]], inverse[s], metric,
                              lower, maximum, tolerance=1e-12, max_directions=len(i))
        for omega in (0., .31*maximum, maximum):
            assert_allclose(permute.schur(omega), initial.schur(omega), atol=5e-11, rtol=5e-11)
            u = permute.reconstruire(omega, np.array([.7, -.2]))
            assert_allclose(u[inverse], initial.reconstruire(omega, np.array([.7, -.2])),
                            atol=5e-9, rtol=5e-9)

    def test_chaine_401_sans_conversion_dense_et_cible_analytique(self):
        n = 401
        data = _chaine(n)
        anciens = {cls: cls.toarray for cls in (csr_matrix, csc_matrix)}

        def garde(cls):
            def conversion(a, *args, **kwargs):
                if min(a.shape) > 64:
                    raise AssertionError("conversion dense d'une grande matrice intérieure")
                return anciens[cls](a, *args, **kwargs)
            return conversion

        with patch.object(csr_matrix, "toarray", garde(csr_matrix)), \
                patch.object(csc_matrix, "toarray", garde(csc_matrix)):
            r = _releve(data, tolerance=1e-10, max_blocs=12, max_directions=64)
            self.assertLessEqual(r.borne_uniforme, 1e-10)
            self.assertLess(r.taille_interieure, 32)
            for omega in np.linspace(0., r.omega_max, 17):
                oracle = n*reference_chaine(n, omega)["schur"]
                self.assertLessEqual(abs(r.schur(omega)[0, 0]-oracle), 1e-10)
            self.assertIs(r.audit_modele()["certification_machine"], False)

    def test_tolerance_sous_arrondi_reste_un_diagnostic(self):
        data = _petit(94, ni=5)
        r = _releve(data, tolerance=1e-28, max_blocs=8, max_directions=5)
        self.assertIs(r.audit_modele()["certification_machine"], False)
        self.assertTrue(np.isfinite(r.borne_uniforme))
        # Une comparaison flottante garde son propre plancher, même lorsque
        # le résidu théorique au carré est beaucoup plus petit.
        d, m, i, s, _, _, maximum = data
        exact = _schur(d.T@d, m, i, s, r.qr.w, maximum)
        assert_allclose(r.schur(maximum), exact, atol=2e-12, rtol=2e-12)


class AssemblagePortsTests(unittest.TestCase):
    def _assemblage(self):
        donnees = [_petit(103, ni=5), _petit(105, ni=4)]
        maximum = min(a[-1] for a in donnees)*.6
        reductions = []
        for d, m, i, s, metric, lower, _ in donnees:
            reductions.append(PortsReleves(d, m, i, s, metric, lower, maximum,
                                           tolerance=1e-20, max_blocs=1, max_directions=1))
        applications = [np.array([[1., 0., 0.], [0., 1., 0.]]),
                        np.array([[0., 1., 0.], [0., 0., 1.]])]
        return donnees, reductions, applications, maximum

    def test_deux_sous_structures_contre_elimination_globale(self):
        donnees, reductions, applications, maximum = self._assemblage()
        ni = sum(len(d[2]) for d in donnees)
        k, m = np.zeros((ni+3, ni+3)), np.zeros((ni+3, ni+3))
        decalage = 0
        for data, reduction, e in zip(donnees, reductions, applications):
            dl, ml, il, sl, _, _, _ = data
            injection = np.zeros((len(il)+len(sl), ni+3))
            injection[il, decalage+np.arange(len(il))] = 1.
            injection[np.ix_(sl, ni+np.arange(3))] = e
            k += injection.T@(dl.T@dl)@injection
            m += injection.T@ml@injection
            decalage += len(il)
        ke, me = np.diag([2., 3., 4.]), np.diag([.1, .2, .3])
        k[ni:, ni:] += ke
        m[ni:, ni:] += me
        assemblage = AssemblagePorts(reductions, applications, ke, me)
        force = np.array([.3, -.7, 1.1])
        for omega in (0., .4*maximum, maximum):
            exact = _schur(k, m, np.arange(ni), ni+np.arange(3), np.eye(3), omega)
            reduit = assemblage.schur(omega)
            ecart = reduit-exact
            # Le défaut QR→D_original n'a pas de signe imposé : l'enveloppe
            # d'assemblage est bilatérale pour le Schur effectivement renvoyé.
            self.assertGreaterEqual(np.linalg.eigvalsh(assemblage.erreur_matrice+ecart)[0], -3e-12)
            self.assertGreaterEqual(np.linalg.eigvalsh(assemblage.erreur_matrice-ecart)[0], -3e-12)
            reponse = assemblage.reponse(omega, force)
            self.assertIs(reponse["certification_machine"], False)
            self.assertGreater(reponse["marge"], 0.)
            erreur = np.linalg.norm(reponse["deplacement"]-np.linalg.solve(exact, force))
            self.assertLessEqual(erreur, reponse["borne_erreur_norme"]+3e-12)

    def test_coercivite_interieure_ne_suffit_pas_pres_resonance_globale(self):
        _, reductions, applications, maximum = self._assemblage()
        assemblage = AssemblagePorts(reductions, applications)
        omega = .7*maximum
        delta = .5*assemblage.borne_uniforme
        self.assertGreater(delta, 0.)
        masse_support = (assemblage.schur(omega)-delta*np.eye(3))/omega**2
        self.assertGreater(np.linalg.eigvalsh(masse_support)[0], 0.)
        proche = AssemblagePorts(reductions, applications, m_externe=masse_support)
        reponse = proche.reponse(omega, np.ones(3))
        self.assertLessEqual(reponse["marge"], 0.)
        self.assertIsNone(reponse["borne_erreur_norme"])
        self.assertIs(reponse["certification_machine"], False)


if __name__ == "__main__":
    unittest.main()
