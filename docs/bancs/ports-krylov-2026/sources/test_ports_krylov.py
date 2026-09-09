"""Contre-calculs des transferts par ports, indépendants des récurrences Krylov.

Les solveurs spectraux et directs denses ne servent qu'aux petits oracles.
Les marges d'arrondi des assertions sont séparées des bornes mathématiques
évaluées par le prototype ; ces tests ne constituent pas une certification
en arithmétique flottante.
"""
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import block_diag, eigh, solve
from scipy.sparse import csc_matrix, diags

from modeles_ports import (borne_lambda_trace_console, chaine, console,
                          metrique_ports, reference_chaine)
from ports_krylov import InterieurKrylov, reduire


def donnees_spd(n=17, ports=3, graine=807):
    rng = np.random.default_rng(graine)
    q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    p, _ = np.linalg.qr(rng.normal(size=(n, n)))
    k = (q*np.geomspace(1., 30., n))@q.T
    m = (p*np.linspace(.6, 1.4, n))@p.T
    b = rng.normal(size=(n, ports))/np.sqrt(n)
    # Oracle spectral petit, extérieur à l'algorithme essayé.
    lam = .9*float(eigh(k, m, eigvals_only=True)[0])
    return k, m, b, lam


def transfert_direct(k, m, b, omega):
    return b.T@solve(k-omega**2*m, b, assume_a='pos')


def marge_arrondi(*matrices):
    """Tolérance de comparaison numérique, distincte du majorant essayé."""
    return 512*np.finfo(float).eps*max(
        1., *(float(np.linalg.norm(a, 2)) for a in matrices))


class PortsKrylov(unittest.TestCase):
    def test_spd_multiport_reference_dense_et_directions_conjointes(self):
        k, m, b, lam = donnees_spd()
        modele = InterieurKrylov(k, m, b, lam, .7*np.sqrt(lam),
                                tolerance=2e-11, max_blocs=12,
                                max_directions=len(k))
        self.assertEqual(modele.statut, 'tolerance_estimee')
        self.assertLessEqual(modele.borne_uniforme, modele.tolerance)
        for omega in np.linspace(0., modele.omega_max, 9):
            with self.subTest(omega=omega):
                exact = transfert_direct(k, m, b, omega)
                obtenu = modele.transfert(omega)
                erreur = exact-obtenu
                arrondi = marge_arrondi(exact, obtenu)
                self.assertLessEqual(np.linalg.norm(erreur, 2),
                                     modele.borne_uniforme+arrondi)
                self.assertGreaterEqual(eigh(erreur, eigvals_only=True)[0], -arrondi)
                # Une combinaison des ports vérifie aussi l'énergie physique.
                poids = np.array([.3, -1.2, .7])
                charge = b@poids
                solution = solve(k-omega**2*m, charge, assume_a='pos')
                self.assertAlmostEqual(float(poids@exact@poids),
                                       float(charge@solution), places=13)

    def test_ports_dependants_et_colonne_nulle(self):
        k, m = np.eye(7), np.eye(7)
        b = np.column_stack((k[:, 0], k[:, 1], 2*k[:, 0],
                             np.zeros(7), k[:, 1]-k[:, 0]))
        modele = InterieurKrylov(k, m, b, 1., .6, tolerance=1e-20)
        self.assertEqual(modele.taille_interieure, 2)
        self.assertEqual(modele.statut, 'tolerance_estimee')
        for omega in (0., .3, .6):
            obtenu = modele.transfert(omega)
            np.testing.assert_allclose(obtenu, b.T@b/(1-omega**2),
                                       rtol=5e-14, atol=5e-14)
            np.testing.assert_array_equal(obtenu[3], np.zeros(5))

    def test_tous_ports_nuls_espace_vide(self):
        k, m, b, lam = donnees_spd(8, 4)
        modele = InterieurKrylov(k, m, np.zeros_like(b), lam, .7*np.sqrt(lam))
        self.assertEqual(modele.taille_interieure, 0)
        self.assertEqual(modele.borne_uniforme, 0.)
        self.assertEqual(modele.statut, 'tolerance_estimee')
        for omega in (0., modele.omega_max):
            np.testing.assert_array_equal(modele.transfert(omega), np.zeros((4, 4)))
            np.testing.assert_array_equal(modele.borne_ponctuelle(omega),
                                          np.zeros((4, 4)))

    def test_identite_energetique_apres_resolution_reduite_perturbee(self):
        k, m, b, lam = donnees_spd(11, 2, 713)
        modele = InterieurKrylov(k, m, b, lam, .8*np.sqrt(lam), max_blocs=2)
        omega = .73*np.sqrt(lam)
        y = modele.coefficients(omega)
        y += .04*np.random.default_rng(21).normal(size=y.shape)
        x = modele.w@y
        dynamique = k-omega**2*m
        residu = b-dynamique@x
        erreur_energetique = residu.T@solve(dynamique, residu, assume_a='pos')
        exact = transfert_direct(k, m, b, omega)
        # Injection d'une erreur de résolution, sans réutiliser l'identité
        # réduite : le fonctionnel corrigé doit rester valable pour tout X.
        with patch.object(modele, 'coefficients', return_value=y):
            obtenu = modele.transfert(omega)
        np.testing.assert_allclose(exact-obtenu, erreur_energetique,
                                   rtol=2e-12, atol=3e-14)
        self.assertGreater(np.linalg.norm(obtenu-modele.d.T@y, 2), 1e-4)

    def test_borne_uniforme_grille_dense_et_borne_ponctuelle(self):
        k, m, b, lam = donnees_spd(19, 3, 660)
        modele = InterieurKrylov(k, m, b, lam, .91*np.sqrt(lam),
                                tolerance=1e-18, max_blocs=3)
        maximum_observe = 0.
        for omega in np.linspace(0., modele.omega_max, 161):
            exact = transfert_direct(k, m, b, omega)
            obtenu = modele.transfert(omega)
            erreur = (exact-obtenu+exact.T-obtenu.T)*.5
            arrondi = marge_arrondi(exact, obtenu)
            maximum_observe = max(maximum_observe, np.linalg.norm(erreur, 2))
            self.assertLessEqual(np.linalg.norm(erreur, 2),
                                 modele.borne_uniforme+arrondi)
            self.assertGreaterEqual(eigh(erreur, eigvals_only=True)[0], -arrondi)
            ponctuelle = modele.borne_ponctuelle(omega)
            self.assertGreaterEqual(eigh(ponctuelle-erreur,
                                        eigvals_only=True)[0], -arrondi)
        self.assertGreater(maximum_observe, 1e-7)
        self.assertGreater(modele.borne_uniforme, modele.tolerance)

    def test_budget_de_directions_conserve_erreur_de_semence(self):
        modele = InterieurKrylov(np.eye(3), np.eye(3), np.eye(3), 1., .5,
                                max_directions=1, tolerance=1e-10)
        self.assertEqual(modele.taille_interieure, 1)
        self.assertEqual(modele.statut, 'budget_directions')
        erreur_statique = np.eye(3)-modele.transfert(0.)
        self.assertGreater(np.linalg.norm(erreur_statique, 2), .99)
        self.assertGreaterEqual(modele.borne_uniforme, .99)
        self.assertGreater(modele.borne_uniforme, modele.tolerance)

    def test_budget_de_blocs_ne_vaut_pas_tolerance(self):
        k, m, b, lam = donnees_spd(11, 2)
        modele = InterieurKrylov(k, m, b, lam, .85*np.sqrt(lam),
                                max_blocs=1, tolerance=1e-16)
        self.assertEqual(modele.statut, 'budget_blocs')
        self.assertGreater(modele.borne_uniforme, modele.tolerance)
        self.assertGreater(np.linalg.norm(transfert_direct(k, m, b, modele.omega_max)
                                          -modele.transfert(modele.omega_max)), 1e-4)

    def test_deflation_de_port_faible_ne_vaut_pas_tolerance(self):
        k = np.eye(4)
        b = np.column_stack((k[:, 0], k[:, 0]+1e-3*k[:, 1]))
        modele = InterieurKrylov(k, k, b, 1., .4, seuil_dependance=1e-3,
                                tolerance=1e-14)
        self.assertEqual(modele.taille_interieure, 1)
        self.assertEqual(modele.statut, 'stagnation')
        self.assertGreater(modele.borne_uniforme, 1e-8)
        self.assertGreater(np.linalg.norm(b.T@b-modele.transfert(0.), 2), 1e-8)

    def test_partition_permutation_et_reconstruction_physique(self):
        kii, mii, couplage, lam = donnees_spd(9, 3, 19)
        metrique = np.array([[3., .4, -.2], [.4, 2., .1], [-.2, .1, 1.]])
        kss = metrique+couplage.T@solve(kii, couplage, assume_a='pos')
        k = np.block([[kii, couplage], [couplage.T, kss]])
        m = block_diag(mii, np.diag([.7, 1.2, .9]))
        i, s = np.arange(9), np.arange(9, 12)
        options = dict(lambda_min=lam, omega_max=.7*np.sqrt(lam),
                       tolerance=1e-20, max_directions=9)
        original = reduire(k, m, i, s, metrique, **options)
        permutation = np.random.default_rng(17).permutation(12)
        inverse = np.argsort(permutation)
        ports = np.array([2, 0, 1])
        permute = reduire(k[np.ix_(permutation, permutation)],
                         m[np.ix_(permutation, permutation)],
                         inverse[i][::-1], inverse[s][ports],
                         metrique[np.ix_(ports, ports)], **options)
        for omega in (0., .3*np.sqrt(lam), .7*np.sqrt(lam)):
            for modele, ordre in ((original, np.arange(3)), (permute, ports)):
                # Retour aux unités des ports : S = W^-T S_normalise W^-1.
                w = modele.normalisation
                gauche = solve(w.T, modele.schur(omega))
                physique = solve(w.T, gauche.T).T
                exact = (kss-omega**2*m[9:, 9:]
                         -transfert_direct(kii, mii, couplage, omega))
                np.testing.assert_allclose(physique, exact[np.ix_(ordre, ordre)],
                                           rtol=2e-12, atol=2e-12)
            port_physique = np.array([.2, -.4, .7])
            champ = original.reconstruire(omega,
                                         solve(original.normalisation, port_physique))
            champ_permute = permute.reconstruire(
                omega, solve(permute.normalisation, port_physique[ports]))
            np.testing.assert_allclose(champ_permute, champ[permutation],
                                       rtol=2e-12, atol=2e-12)
            np.testing.assert_allclose((k-omega**2*m)@champ,
                                       np.r_[np.zeros(9), exact@port_physique],
                                       rtol=3e-12, atol=3e-12)

    def test_chaine_oracle_sinus_et_champ_impose(self):
        n, raideur, masse = 25, 7., .4
        fixture = chaine(n, raideur, masse)
        metrique = metrique_ports(fixture)
        lam = fixture.metadata['mu_min_interieur_s_moins_2']*(1-1e-12)
        modele = reduire(fixture.k, fixture.m, fixture.interieur, fixture.interface,
                         metrique, lambda_min=lam, omega_max=.88*np.sqrt(lam),
                         tolerance=2e-10, max_blocs=16)
        self.assertEqual(modele.interieur.statut, 'tolerance_estimee')
        for proportion in (0., .13, .37, .71, .88):
            omega = proportion*np.sqrt(lam)
            oracle = reference_chaine(n, omega, raideur, masse)
            schur_physique = float(modele.schur(omega)[0, 0]*metrique[0, 0])
            self.assertAlmostEqual(schur_physique, oracle['schur'], delta=2e-10)
            port = solve(modele.normalisation, np.ones(1))
            champ = modele.reconstruire(omega, port)
            np.testing.assert_allclose(champ, oracle['rapports_deplacement'],
                                       rtol=3e-8, atol=3e-8)

    def test_console_native_petite_metrique_et_schur_dense(self):
        fixture = console(4)
        i, s = fixture.interieur, fixture.interface
        metrique = metrique_ports(fixture)
        lam = borne_lambda_trace_console(fixture.metadata)*(1-1e-10)
        kii = fixture.k[i][:, i].toarray()
        mii = fixture.m[i][:, i].toarray()
        self.assertLessEqual(lam, float(eigh(kii, mii, eigvals_only=True)[0]))
        modele = reduire(fixture.k, fixture.m, i, s, metrique,
                         lambda_min=lam, omega_max=.5*np.sqrt(lam),
                         tolerance=2e-9, max_directions=len(i))
        for omega in (0., modele.interieur.omega_max):
            b = fixture.k[i][:, s].toarray()@modele.normalisation
            exact = modele.kss-omega**2*modele.mss-transfert_direct(kii, mii, b, omega)
            np.testing.assert_allclose(modele.schur(omega), exact,
                                       rtol=1e-9, atol=2e-9)
        np.testing.assert_allclose(modele.schur(0.), np.eye(6),
                                   rtol=2e-9, atol=2e-9)

    def test_bande_et_contradiction_spectrale_detectee(self):
        for omega in (1., 1.1):
            with self.subTest(omega=omega), self.assertRaisesRegex(ValueError, 'bande'):
                InterieurKrylov(np.eye(2), np.eye(2), np.ones((2, 1)),
                                1., omega)
        with self.assertRaisesRegex(ValueError, 'Rayleigh'):
            InterieurKrylov(np.diag([1., 3.]), np.eye(2), np.array([[1.], [0.]]),
                            1.25, .5)
        modele = InterieurKrylov(np.eye(2), np.eye(2), np.ones((2, 1)), 1., .5)
        for omega in (-.1, .50001, np.nan, np.inf):
            for methode in (modele.transfert, modele.coefficients, modele.borne_ponctuelle):
                with self.subTest(omega=omega, methode=methode.__name__):
                    with self.assertRaisesRegex(ValueError, 'bande'):
                        methode(omega)

    def test_domaines_non_pris_en_charge_refuses(self):
        k, m = np.eye(3), np.eye(3)
        with self.assertRaisesRegex(TypeError, 'quadrature'):
            InterieurKrylov(k, m, np.ones((3, 1)), 1., .5,
                            quadrature='radau')
        for i, s in (([0, 1], [1]), ([0], [2]), ([0., 1.], [2.])):
            with self.subTest(i=i, s=s), self.assertRaisesRegex(ValueError, 'partitionner'):
                reduire(k, m, i, s, np.eye(1), lambda_min=1., omega_max=.5)
        m[0, 2] = m[2, 0] = .1
        with self.assertRaisesRegex(ValueError, 'M_IS'):
            reduire(k, m, [0, 1], [2], np.eye(1), lambda_min=1., omega_max=.5)

    def test_grande_chaine_sans_densification_de_matrice_originale(self):
        n = 401
        k = diags((-np.ones(n-1), 2*np.ones(n), -np.ones(n-1)),
                  (-1, 0, 1), format='csc')
        m = diags(np.ones(n), format='csc')
        b = np.zeros((n, 1))
        b[-1, 0] = 1.
        lam = 4*np.sin(np.pi/(2*(n+1)))**2*(1-1e-12)
        with patch.object(csc_matrix, 'toarray', side_effect=AssertionError('densification')):
            modele = InterieurKrylov(k, m, b, lam, .4*np.sqrt(lam),
                                    tolerance=1e-9, max_blocs=8)
            valeur = modele.transfert(modele.omega_max)
        self.assertEqual(valeur.shape, (1, 1))
        self.assertLess(modele.taille_interieure, 20)
        self.assertLess(modele.nnz_facteurs, 8*n)


if __name__ == '__main__':
    unittest.main()
