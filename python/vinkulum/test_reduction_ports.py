"""Contrats physiques de l'API publique, exécutables depuis une roue installée.

Les références denses portent seulement sur de petits systèmes complets.
Les constantes spectrales fournies et calculées automatiquement sont testées
séparément ; les marges flottantes ne sont pas des certificats machine.
"""
import copy
from decimal import Decimal, localcontext
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import block_diag, eigh, solve
from scipy.sparse import csc_matrix, lil_matrix
from scipy.spatial.transform import Rotation

from vinkulum import Noyau
from vinkulum.reduction_ports import ReductionMaterielle, reduire_poutres


FLOTTANT = 3e-10
LONGUEUR = .6
EA, GAY, GAZ, GJ, EIY, EIZ = 3200., 180., 240., 24., 18., 9.
INERTIE = np.diag([.12, .16, .20])


def _csc(donnees):
    lignes, colonnes, ptr, ind, valeurs = donnees
    return csc_matrix((valeurs, ind, ptr), shape=(lignes, colonnes))


def _poids(formulation, longueur=LONGUEUR):
    gy, gz = GAY, GAZ
    if formulation == "integree":
        gy = 1 / (1 / gy + longueur**2 / (12 * EIZ))
        gz = 1 / (1 / gz + longueur**2 / (12 * EIY))
    return np.sqrt(longueur * np.array([EA, gy, gz, GJ, EIY, EIZ]))


def _console(nb=3, repere=None, formulation="integree"):
    q = np.eye(3) if repere is None else repere
    n = Noyau([0., 0., 0.])
    for j in range(nb + 1):
        n.corps(str(j), 1. + .1 * j, INERTIE.ravel().tolist(),
                (q @ [j * LONGUEUR, 0., 0.]).tolist(), rot=q.ravel().tolist())
    n.liaison("encastrement", None, 0, ra=q.ravel().tolist(),
              bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    for j in range(nb):
        n.poutre(str(j), j, j + 1, EA, GAY, GJ, EIY, ei3=EIZ, ga3=GAZ,
                 formulation=formulation)
    return n


def _facteur_droit(nb, formulation):
    """Dérivées des six déformations linéaires, écrites sans le noyau."""
    b = np.zeros((6, 12))
    b[:, :6], b[:, 6:] = -np.eye(6) / LONGUEUR, np.eye(6) / LONGUEUR
    b[1, [5, 11]] = -.5
    b[2, [4, 10]] = .5
    local = _poids(formulation)[:, None] * b
    d = np.zeros((6 * nb, 6 * (nb + 1)))
    for j in range(nb):
        d[6*j:6*j+6, 6*j:6*j+12] = local
    return d


def _donnees(n, libres=None, ports=None):
    natif = n.facteurs_materiels_poutres()
    d, m = (_csc(natif[k]).toarray() for k in ("d", "masse"))
    libres = np.arange(6, d.shape[1]) if libres is None else np.asarray(libres)
    ports = np.arange(d.shape[1] - 6, d.shape[1]) if ports is None else np.asarray(ports)
    positions = {int(v): j for j, v in enumerate(libres)}
    s = np.array([positions[int(v)] for v in ports])
    i = np.setdiff1d(np.arange(len(libres)), s)
    return d[:, libres], m[np.ix_(libres, libres)], libres, ports, i, s


def _bande(d, m, i):
    k = d.T @ d
    mini = eigh(k[np.ix_(i, i)], m[np.ix_(i, i)], eigvals_only=True)[0]
    globale = eigh(k, m, eigvals_only=True)[0]
    return .8 * mini, .25 * np.sqrt(globale)


def _reduire(n, libres=None, ports=None, metrique=None, **options):
    d, m, libres, ports, i, s = _donnees(n, libres, ports)
    lower, omega = _bande(d, m, i)
    if metrique is None:
        metrique = np.diag(np.linspace(.4, 2., len(s)))
    reglages = dict(lambda_min=lower, tolerance=1e-11, max_blocs=12,
                   max_directions=len(i))
    reglages.update(options)
    r = reduire_poutres(n, libres, ports, metrique, omega, **reglages)
    return r, d, m, s


def _exacte(d, m, s, omega, force):
    rhs = np.zeros(d.shape[1])
    rhs[s] = force
    return solve(d.T @ d - omega**2 * m, rhs, assume_a="pos")


class FacteursMaterielsPublicsTests(unittest.TestCase):
    def test_console_anisotrope_facteur_analytique_et_tangente_au_repos(self):
        for formulation in ("milieu", "integree"):
            with self.subTest(formulation=formulation):
                n = _console(formulation=formulation)
                f = n.facteurs_materiels_poutres()
                d = _csc(f["d"]).toarray()
                attendu = _facteur_droit(3, formulation)
                assert_allclose(d, attendu, rtol=2e-14, atol=2e-13)
                k, c, m, _, g = map(np.asarray, n.k_c_m_z())
                assert_allclose(d.T @ d, k, rtol=3e-13, atol=3e-11)
                assert_allclose(_csc(f["masse"]).toarray(), m, atol=2e-15)
                assert_allclose(_csc(f["contraintes"]).toarray(), g, atol=2e-15)
                assert_allclose(f["deformations"], 0., atol=2e-14)
                assert_allclose(c, 0., atol=2e-14)
                self.assertEqual(f["poutres"], [(str(j), j, j + 1) for j in range(3)])

    def test_repere_tourne_conserve_facteur_et_masse_physiques(self):
        q = Rotation.from_rotvec([.37, -.42, .28]).as_matrix()
        f0 = _console().facteurs_materiels_poutres()
        fq = _console(repere=q).facteurs_materiels_poutres()
        t = block_diag(*([q] * 8))
        d0, dq = (_csc(f["d"]).toarray() for f in (f0, fq))
        m0, mq = (_csc(f["masse"]).toarray() for f in (f0, fq))
        assert_allclose(dq @ t, d0, rtol=2e-13, atol=2e-12)
        assert_allclose(t.T @ mq @ t, m0, rtol=2e-14, atol=2e-15)
        self.assertGreater(np.linalg.norm(mq[3:6, 3:6] - np.diag(np.diag(mq[3:6, 3:6]))), .01)

    def test_branche_3d_assemblage_et_six_mouvements_rigides(self):
        positions = np.array([[.4, .7, -.2], [0., 0., 0.], [.8, -.3, .5], [-.4, .2, .8]])
        liens = [(2, 0), (1, 2), (3, 2)]
        n = Noyau([0., 0., 0.])
        for j, p in enumerate(positions):
            q = Rotation.from_rotvec([.07*j, -.04*j, .11*j]).as_matrix()
            n.corps(str(j), 1., INERTIE.ravel().tolist(), p.tolist(), rot=q.ravel().tolist())
        for j, (a, b) in enumerate(liens):
            n.poutre(str(j), a, b, EA, GAY, GJ, EIY, ei3=EIZ, ga3=GAZ,
                     formulation="integree")
        f = n.facteurs_materiels_poutres()
        d = _csc(f["d"]).toarray()
        assert_allclose(d.T @ d, n.k_c_m_z()[0], rtol=2e-12, atol=3e-11)
        for j, (a, b) in enumerate(liens):
            autres = np.setdiff1d(np.arange(24), np.r_[6*a:6*a+6, 6*b:6*b+6])
            assert_allclose(d[6*j:6*j+6, autres], 0., atol=0.)
        rigides = np.zeros((24, 6))
        for j, p in enumerate(positions):
            rigides[6*j:6*j+3, :3] = np.eye(3)
            rigides[6*j:6*j+3, 3:] = np.column_stack([np.cross(e, p) for e in np.eye(3)])
            rigides[6*j+3:6*j+6, 3:] = np.eye(3)
        assert_allclose(d @ rigides, 0., atol=3e-13)
        self.assertEqual(np.linalg.matrix_rank(d, tol=1e-10), 18)

    def test_precontrainte_energie_et_gradient_distincts_de_la_tangente(self):
        for formulation in ("milieu", "integree"):
            with self.subTest(formulation=formulation):
                n = _console(nb=1, formulation=formulation)
                state = n.etat()
                state[1][1][0] += .08
                state[1][1][1] += .04
                state[2][1] = Rotation.from_rotvec([.04, -.07, .11]).as_matrix().ravel().tolist()
                n.pose_etat(*state)
                f = n.facteurs_materiels_poutres()
                d, z = _csc(f["d"]).toarray(), np.asarray(f["deformations"])
                _, gamma, kappa = n.poutres()[0]
                attendu = _poids(formulation) * np.r_[gamma, kappa]
                assert_allclose(z, attendu, rtol=3e-14, atol=2e-14)
                self.assertAlmostEqual(.5 * z @ z, .5 * attendu @ attendu, places=12)
                assert_allclose(d.T @ z, n.residu_statique(), rtol=3e-13, atol=2e-11)
                k = np.asarray(n.k_c_m_z()[0])
                self.assertGreater(np.linalg.norm(k - d.T @ d) / np.linalg.norm(k), 1e-3)

    def test_extraction_copie_etat_sans_avancer_ni_modifier_le_noyau(self):
        n = _console(nb=1)
        state = n.etat()
        state[3][1][0], state[4][1][2] = .3, -.2
        n.pose_etat(*state)
        avant = copy.deepcopy(n.etat_precis())
        f = n.facteurs_materiels_poutres(t=.75)
        self.assertEqual(n.etat_precis(), avant)
        self.assertEqual(f["t"], .75)
        self.assertEqual(f["etat"]["temps_modele"], avant[0])
        self.assertTrue(any("vitesses" in x for x in f["infos_domaine"]))
        reference = copy.deepcopy(n.facteurs_materiels_poutres(t=.75))
        f["d"][4][0] += 7.
        f["etat"]["corps"][1]["position"][0] += 9.
        f["deformations"][0] = 12.
        self.assertEqual(n.facteurs_materiels_poutres(t=.75), reference)
        self.assertEqual(n.etat_precis(), avant)

    def test_super_element_signale_et_exclu_du_facteur_materiel(self):
        n = _console(nb=1)
        avant = n.facteurs_materiels_poutres()
        supplement = np.zeros((12, 12))
        supplement[np.ix_([0, 6], [0, 6])] = 37. * np.array([[1., -1.], [-1., 1.]])
        n.superelement("ressort_externe", [0, 1], supplement.ravel().tolist())
        apres = n.facteurs_materiels_poutres()
        self.assertEqual(apres["d"], avant["d"])
        self.assertTrue(any("superéléments" in x for x in apres["infos_domaine"]))
        self.assertIn("ressort_externe", apres["etat"]["superelements"])
        d = _csc(apres["d"]).toarray()
        assert_allclose(np.asarray(n.k_c_m_z()[0]) - d.T @ d, supplement, atol=2e-11)


class ReductionMateriellePubliqueTests(unittest.TestCase):
    def _verifie_reponse(self, r, d, m, s, omega, force):
        result = r.reponse(omega, force)
        exact = _exacte(d, m, s, omega, force)
        candidat = result["champ_physique"]
        assert_allclose(candidat, exact, rtol=FLOTTANT, atol=FLOTTANT)
        assert_allclose(result["deplacement_ports_physique"], exact[s],
                        rtol=FLOTTANT, atol=FLOTTANT)
        err = candidat - exact
        for nom, norme in (("masse", lambda x: np.sqrt(x @ m @ x)),
                           ("deformation", lambda x: np.linalg.norm(d @ x))):
            mesure = result[nom]
            self.assertIsNotNone(mesure["borne_absolue"])
            self.assertLessEqual(norme(err), mesure["borne_absolue"] + FLOTTANT)
            self.assertAlmostEqual(mesure["norme_candidate"], norme(candidat), delta=FLOTTANT)
            if mesure["borne_relative"] is not None:
                self.assertLessEqual(norme(err) / norme(exact), mesure["borne_relative"] + FLOTTANT)
        self.assertIs(result["certification_machine"], False)
        return result

    def test_reduction_native_qr_et_lu_contre_systeme_complet(self):
        force = np.array([.8, -.3, .6, .04, -.07, .09])
        for methode in ("qr", "lu_energie"):
            with self.subTest(methode=methode):
                r, d, m, s = _reduire(_console(), methode=methode)
                self.assertIsNone(r.certificat_spectral)
                self.assertFalse(r.origine["tangente_globale"])
                for omega in (0., .4 * r.omega_max, r.omega_max):
                    self._verifie_reponse(r, d, m, s, omega, force)

    def test_unites_metrique_et_ordre_des_ports_preservent_le_travail(self):
        n = _console()
        ordre = np.array([4, 0, 5, 2, 1, 3])
        ports = np.arange(18, 24)[ordre]
        force = np.array([.06, .5, -.03, .7, -.2, .04])
        z = np.array([[2., .3, 0., 0., 0., 0.], [0., .7, .1, 0., 0., 0.],
                      [0., 0., 3., .2, 0., 0.], [0., 0., 0., .4, .1, 0.],
                      [0., 0., 0., 0., 4., .3], [0., 0., 0., 0., 0., .8]])
        resultats = []
        for metric in (np.eye(6), z.T @ z):
            r, d, m, s = _reduire(n, ports=ports, metrique=metric)
            result = self._verifie_reponse(r, d, m, s, r.omega_max, force)
            y = result["coordonnees_ports_normalisees"]
            u = result["deplacement_ports_physique"]
            self.assertAlmostEqual(force @ u, (r.qr.w.T @ force) @ y, delta=2e-13)
            self.assertAlmostEqual(u @ metric @ u, y @ y, delta=2e-13)
            resultats.append(result["champ_physique"])
        assert_allclose(*resultats, rtol=FLOTTANT, atol=FLOTTANT)

    def test_selection_permutee_conserve_sous_matrice_masse_couplee(self):
        q = Rotation.from_rotvec([.7, -.4, .3]).as_matrix()
        n = _console(nb=2, repere=q)
        libres = np.array([17, 8, 10, 6, 15, 12, 9, 16, 13, 11, 7, 14])
        ports = np.array([16, 12, 14])
        r, d, m, s = _reduire(n, libres=libres, ports=ports)
        self.assertGreater(np.linalg.norm(m[np.ix_(r.qr.i, s)]), .01)
        assert_allclose(r.m.toarray(), m, atol=0.)
        np.testing.assert_array_equal(r.coordonnees_physiques, libres)
        result = self._verifie_reponse(r, d, m, s, r.omega_max, np.array([.05, .4, -.2]))
        np.testing.assert_array_equal(result["coordonnees_physiques"], libres)
        libres[0] = 6
        self.assertEqual(r.coordonnees_physiques[0], 17)

    def test_constructeur_matriciel_garde_tous_les_couplages_masse(self):
        rng = np.random.default_rng(706)
        d = np.vstack((2 * np.eye(9), .4 * rng.normal(size=(12, 9))))
        z = .4 * rng.normal(size=(11, 9))
        m = .7 * np.eye(9) + z.T @ z
        s, i = np.array([7, 2]), np.array([0, 1, 3, 4, 5, 6, 8])
        lower, omega = _bande(d, m, i)
        r = ReductionMaterielle(d, m, i, s, np.array([[2., .3], [.3, .7]]),
                                omega, lambda_min=lower, tolerance=1e-12,
                                max_directions=len(i))
        force = np.array([.7, -.4])
        self._verifie_reponse(r, d, m, s, omega, force)
        self.assertGreater(np.linalg.norm(m[np.ix_(i, s)]), .1)
        mauvaise = m.copy()
        mauvaise[np.ix_(i, s)] = 0.
        mauvaise[np.ix_(s, i)] = 0.
        self.assertGreater(np.linalg.norm(_exacte(d, m, s, omega, force)
                                          - _exacte(d, mauvaise, s, omega, force)), 1e-4)

    def test_mutation_noyau_ne_change_ni_modele_extrait_ni_resultat_passe(self):
        n = _console(nb=2)
        r, d, m, s = _reduire(n)
        force = np.array([.3, -.6, .2, .08, .03, -.04])
        result = r.reponse(r.omega_max, force)
        sauvegarde = copy.deepcopy(result)
        origine = copy.deepcopy(r.origine)
        state = n.etat()
        state[1][2][1] += .2
        state[2][2] = Rotation.from_rotvec([.1, .2, -.15]).as_matrix().ravel().tolist()
        n.pose_etat(*state)
        self.assertNotEqual(n.facteurs_materiels_poutres()["deformations"],
                            origine["donnees_natives"]["deformations"])
        self.assertEqual(r.origine, origine)
        for nom in ("champ_physique", "deplacement_ports_physique", "coordonnees_physiques"):
            assert_allclose(result[nom], sauvegarde[nom], rtol=0., atol=0.)
        nouveau = self._verifie_reponse(r, d, m, s, r.omega_max, force)
        assert_allclose(nouveau["champ_physique"], sauvegarde["champ_physique"], atol=0., rtol=0.)
        nouveau["coordonnees_physiques"][0] = -1
        self.assertEqual(r.coordonnees_physiques[0], 6)

    def test_refus_selections_non_admissibles_et_interfaces_invalides(self):
        n = _console(nb=2)
        bons = np.arange(6, 18)
        cas = [(np.arange(18), [12]), (np.r_[bons, 7], [12]),
               (bons, [12, 12]), (bons, [0]), (bons, list(bons)),
               (bons.astype(float), [12])]
        for libres, ports in cas:
            with self.subTest(libres=list(libres), ports=ports):
                with self.assertRaises(ValueError):
                    reduire_poutres(n, libres, ports, np.eye(len(ports)), .01, lambda_min=1.)
        # La distance lie deux translations : leur simple conservation
        # exige une application admissible, pas la suppression de la racine.
        n.distance("diagonale", 1, 2, [0., 0., 0.], [0., .2, 0.])
        with self.assertRaisesRegex(ValueError, "non admissibles"):
            reduire_poutres(n, bons, [12], np.eye(1), .01, lambda_min=1.)
        for nature in ("non_holonome", "pilotee", "phi_non_nul"):
            with self.subTest(contrainte=nature):
                n = _console(nb=2)
                if nature == "non_holonome":
                    n.liaison("vitesse", None, 0, bloque_t=[0], bloque_r=[], nh=True)
                elif nature == "pilotee":
                    n.liaison("commande", None, 0,
                              cible_t=([1., 0., 0.], ("lineaire", [0., 1.])))
                else:
                    state = n.etat()
                    state[1][0][0] += .001
                    n.pose_etat(*state)
                with self.assertRaises(ValueError):
                    reduire_poutres(n, bons, [12], np.eye(1), .01, lambda_min=1.)

    def test_force_et_frequence_hors_contrat_sont_refusees(self):
        r, _, _, _ = _reduire(_console(nb=2))
        for force in (np.ones(5), np.ones((6, 1)), np.r_[np.ones(5), np.nan]):
            with self.subTest(force=np.asarray(force).shape):
                with self.assertRaises(ValueError):
                    r.reponse(0., force)
        for omega in (-.1, 1.01 * r.omega_max, np.nan):
            with self.subTest(omega=omega):
                with self.assertRaises(ValueError):
                    r.reponse(omega, np.ones(6))
        nul = r.reponse(0., np.zeros(6))
        assert_allclose(nul["champ_physique"], 0., atol=0.)
        self.assertIsNone(nul["masse"]["borne_relative"])
        self.assertIsNone(nul["deformation"]["borne_relative"])

    def test_lambda_automatique_natif_tourne_vise_les_donnees_originales(self):
        q = Rotation.from_rotvec([.7, -.4, .3]).as_matrix()
        n = _console(nb=2, repere=q)
        libres = np.array([17, 8, 10, 6, 15, 12, 9, 16, 13, 11, 7, 14])
        r, d, m, s = _reduire(n, libres=libres, ports=[16, 12, 14], lambda_min=None)
        certificat = r.certificat_spectral
        self.assertIs(certificat["certification_machine"], True)
        self.assertIs(certificat["minoration_D_original_etablie"], True)
        self.assertGreater(certificat["operations_decimal"], 0)
        self.assertLess(certificat["eta_superieur"], 1.)
        self.assertEqual(r.lambda_min, certificat["lambda_min"])
        self.assertGreater(r.lambda_min, r.omega_max**2)
        kii = (d.T @ d)[np.ix_(r.qr.i, r.qr.i)]
        mii = m[np.ix_(r.qr.i, r.qr.i)]
        minimum = eigh(kii, mii, eigvals_only=True)[0]
        self.assertLessEqual(r.lambda_min, minimum * (1 + 2e-12))
        resultat = self._verifie_reponse(r, d, m, s, r.omega_max, np.array([.05, .4, -.2]))
        self.assertIs(resultat["certification_machine"], False)
        self.assertIs(r.origine["certification_champ_machine"], False)

    def test_lambda_automatique_independante_du_contexte_decimal_ambiant(self):
        d = np.array([[1., 0., 0.], [0., 2., 0.], [0., 0., 3.], [.25, .5, .125]])
        m = np.array([[2., .25, 0.], [.25, 1., 0.], [0., 0., 1.]])
        resultats = []
        for precision in (6, 100):
            with localcontext() as contexte:
                contexte.prec = precision
                r = ReductionMaterielle(d, m, [0, 1], [2], [[1.]], .1,
                                        precision_spectrale=70, tolerance=1e-12)
                self.assertEqual(contexte.prec, precision)
                resultats.append(r.certificat_spectral)
                self._verifie_reponse(r, d, m, np.array([2]), .1, np.array([1.]))
        for cle in ("lambda_min", "trace_inferieure_decimal", "trace_superieure_decimal",
                    "eta_superieur_decimal", "lambda_inferieur_decimal"):
            self.assertEqual(resultats[0][cle], resultats[1][cle])
        minimum = eigh((d.T @ d)[:2, :2], m[:2, :2], eigvals_only=True)[0]
        self.assertGreater(resultats[0]["lambda_min"], 0.)
        self.assertLessEqual(resultats[0]["lambda_min"], minimum)
        with self.assertRaisesRegex(ValueError, "requiert.*qr"):
            ReductionMaterielle(d, m, [0, 1], [2], [[1.]], .1, methode="lu_energie")
        # Le certificat automatique a un domaine de masse plus restreint
        # que le constructeur avec constante spectrale fournie.
        masse_connexe = np.eye(8) + .125 * np.ones((8, 8))
        with self.assertRaisesRegex(RuntimeError, "bloc connexe de masse"):
            ReductionMaterielle(np.eye(8), masse_connexe, np.arange(7), [7], [[1.]], .1)

    def test_tolerance_contraintes_explicite_accepte_repos_tourne_ou_refuse_strictement(self):
        q = Rotation.from_rotvec([.7, -.4, .3]).as_matrix()
        n = _console(nb=2, repere=q)
        phi = np.asarray(n.facteurs_materiels_poutres()["phi"])
        residu = float(np.max(np.abs(phi)))
        self.assertGreater(residu, 0.)
        self.assertLess(residu, 1e-12)
        r, d, m, s = _reduire(n)
        self.assertEqual(r.origine["residu_contraintes_initial"], residu)
        self.assertEqual(r.origine["tolerance_contraintes"], 1e-12)
        self.assertIs(r.origine["certification_champ_machine"], False)
        self._verifie_reponse(r, d, m, s, r.omega_max,
                              np.array([.3, -.6, .2, .08, .03, -.04]))
        with self.assertRaisesRegex(ValueError, "contraintes doivent être satisfaites"):
            _reduire(n, tolerance_contraintes=0.)
        for tolerance in (-1., np.nan, np.inf):
            with self.subTest(tolerance=tolerance):
                with self.assertRaises(ValueError):
                    _reduire(n, tolerance_contraintes=tolerance)

    def test_lambda_automatique_refuse_doublons_lil_avant_normalisation_qr(self):
        for nom in ("D", "M"):
            with self.subTest(matrice=nom):
                d, m = lil_matrix(np.eye(3)), lil_matrix(np.eye(3))
                matrice = d if nom == "D" else m
                matrice.rows[0], matrice.data[0] = [0, 0], [1., .125]
                avant = (copy.deepcopy(matrice.rows.tolist()),
                         copy.deepcopy(matrice.data.tolist()))
                # Un budget QR presque nul distingue le refus d'entrée
                # d'une erreur tardive après sommation des doublons.
                with self.assertRaisesRegex(ValueError, "canoniques?"):
                    ReductionMaterielle(d, m, [0, 1], [2], [[1.]], .1, budget_qr=1)
                self.assertEqual((matrice.rows.tolist(), matrice.data.tolist()), avant)

    def test_adapter_natif_enrichit_en_force_physique_sans_refaire_facteur_ni_certificat(self):
        metric = np.diag([16., .25, 9., .04, 4., .09])
        r, d, m, s = _reduire(_console(nb=4), metrique=metric,
                              lambda_min=None, max_blocs=1)
        force = np.array([.2, 1., -.4, .03, .05, -.08])
        frequences = np.array([0., .4 * r.omega_max, r.omega_max])
        facteur, certificat = r.qr, r.certificat_spectral
        r_initial = facteur.r.copy()
        certificat_initial = copy.deepcopy(certificat)
        dimension = r.taille_interieure
        avant = r.reponse(r.omega_max, force)
        champ_initial = avant["champ_physique"].copy()
        tolerance = 1e-8
        self.assertGreater(avant["deformation"]["borne_relative"], tolerance)
        resultat = r.adapter(frequences, force, tolerance_relative=tolerance, max_etapes=5)
        self.assertEqual(resultat["statut"], "tolerance_aux_frequences_demandees")
        self.assertFalse(resultat["historique"][0]["accepte"])
        self.assertTrue(resultat["historique"][-1]["accepte"])
        self.assertGreater(r.taille_interieure, dimension)
        self.assertIs(r.qr, facteur)
        self.assertIs(r.certificat_spectral, certificat)
        self.assertEqual(certificat, certificat_initial)
        assert_allclose(facteur.r.toarray(), r_initial.toarray(), rtol=0., atol=0.)
        assert_allclose(avant["champ_physique"], champ_initial, rtol=0., atol=0.)
        self.assertIs(certificat["certification_machine"], True)
        self.assertIs(resultat["certification_machine"], False)
        self.assertNotIn("controle", resultat)
        self.assertEqual(len(resultat["reponses"]), len(frequences))
        for omega, rep in zip(frequences, resultat["reponses"]):
            with self.subTest(omega=omega):
                frais = self._verifie_reponse(r, d, m, s, omega, force)
                for cle in ("champ_physique", "deplacement_ports_physique",
                            "coordonnees_ports_normalisees", "normalisation_ports"):
                    assert_allclose(rep[cle], frais[cle], rtol=0., atol=0.)
                for norme in ("masse", "deformation"):
                    self.assertLessEqual(rep[norme]["borne_relative"], tolerance)
                    self.assertEqual(rep[norme]["borne_relative"], frais[norme]["borne_relative"])
                self.assertIs(rep["certification_machine"], False)
                self.assertNotIn("deplacement", rep)

    def test_borne_ports_normalises_et_budget_adaptation_sans_fausse_acceptation(self):
        metric = np.diag([16., .25, 9., .04, 4., .09])
        r, d, m, s = _reduire(_console(nb=4), metrique=metric,
                              max_blocs=1, max_directions=1)
        force = np.array([.2, 1., -.4, .03, .05, -.08])
        rep = r.reponse(r.omega_max, force)
        exact = _exacte(d, m, s, r.omega_max, force)
        erreur_physique = rep["deplacement_ports_physique"] - exact[s]
        w = rep["normalisation_ports"]
        erreur_normalisee = np.linalg.norm(np.linalg.solve(w, erreur_physique))
        self.assertGreater(erreur_normalisee, 1e-6)
        self.assertAlmostEqual(erreur_normalisee,
                               np.sqrt(erreur_physique @ metric @ erreur_physique), delta=2e-14)
        self.assertGreater(abs(erreur_normalisee - np.linalg.norm(erreur_physique)), 1e-6)
        self.assertLessEqual(erreur_normalisee, rep["borne_ports_normalises"] + FLOTTANT)
        self.assertNotIn("borne_ports", rep)
        self.assertIs(rep["certification_machine"], False)
        resultat = r.adapter([r.omega_max], force, tolerance_relative=1e-12, max_etapes=2)
        self.assertEqual(resultat["statut"], "stagnation_ou_budget_directions")
        self.assertFalse(resultat["historique"][-1]["accepte"])
        self.assertEqual(r.taille_interieure, 1)
        self.assertGreater(resultat["reponses"][0]["deformation"]["borne_relative"], 1e-12)
        # La copie de W renvoyée est descriptive, pas un accès mutable au facteur.
        w[0, 0] *= 2
        assert_allclose(r.reponse(r.omega_max, force)["champ_physique"],
                        rep["champ_physique"], rtol=0., atol=0.)
        complexe = force.astype(complex) + 1j * np.ones(len(force))
        for appel in (lambda: r.reponse(r.omega_max, complexe),
                      lambda: r.adapter([r.omega_max], complexe),
                      lambda: r.adapter(np.array([r.omega_max + 1j]), force)):
            with self.assertRaises((TypeError, ValueError)):
                appel()

    def test_flexion_statique_console_contre_oracle_decimal_timoshenko(self):
        n = _console(nb=3)
        r, _, _, _ = _reduire(n)
        force = np.array([0., 1., 0., 0., 0., 0.])
        champ = r.reponse(0., force)["champ_physique"].reshape(3, 6)
        # Solution physique continue, indépendante de D et du Schur :
        # effort tranchant constant et moment linéaire dans une console.
        with localcontext() as contexte:
            contexte.prec = 70
            h, ei, ga = map(Decimal, (str(LONGUEUR), str(EIZ), str(GAY)))
            longueur = 3 * h
            for j in range(1, 4):
                x = j * h
                fleche = x*x*(3*longueur-x)/(6*ei) + x/ga
                rotation = x*(2*longueur-x)/(2*ei)
                self.assertAlmostEqual(champ[j-1, 1], float(fleche), delta=2e-12)
                self.assertAlmostEqual(champ[j-1, 5], float(rotation), delta=2e-12)
        assert_allclose(champ[:, [0, 2, 3, 4]], 0., atol=2e-13)


if __name__ == "__main__":
    unittest.main()
