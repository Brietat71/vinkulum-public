"""Contre-épreuves des champs et observables dans le modèle énergétique initial.

Les références résolvent seulement de petits systèmes denses complets. La
marge FLOTTANT couvre leurs comparaisons arrondies ; elle ne transforme pas
les majorants conditionnels du prototype en certificats machine.
"""
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import solve

from champ_interieur import ControleChamp, adapter_champ, borne_relative, observable_locale
from ports_releves import AssemblagePorts
from test_ports_releves import _petit, _chaine, _releve


FLOTTANT = 5e-12


def _assemblage_petit(**options):
    """Deux intérieurs disjoints, trois ports globaux et applications physiques."""
    donnees = [_petit(411, ni=7), _petit(412, ni=9)]
    reglages = dict(tolerance=1e-12, max_blocs=1, max_directions=16)
    reglages.update(options)
    structures = [_releve(d, **reglages) for d in donnees]
    applications = [np.array([[1., .2, 0.], [0., 1., -.25]]),
                    np.array([[.1, 0., 1.], [0., 1., .4]])]
    k_ext = np.diag([.12, .15, .08])
    m_ext = np.array([[.03, .004, 0.], [.004, .04, .002], [0., .002, .02]])
    assemblage = AssemblagePorts(structures, applications, k_ext, m_ext)
    taille = 3 + sum(len(d[2]) for d in donnees)
    k, m = np.zeros((taille, taille)), np.zeros((taille, taille))
    k[:3, :3], m[:3, :3] = k_ext, m_ext
    restrictions, debut = [], 3
    for data, application in zip(donnees, applications):
        d, masse, i, s = data[:4]
        restriction = np.zeros((d.shape[1], taille))
        restriction[np.ix_(i, np.arange(debut, debut + len(i)))] = np.eye(len(i))
        restriction[s, :3] = application
        # Assemblage indépendant des matrices réduites et du relèvement.
        dr = d @ restriction
        k += dr.T @ dr
        m += restriction.T @ masse @ restriction
        restrictions.append(restriction)
        debut += len(i)
    return assemblage, k, m, restrictions


def _reponse_dense_et_candidate(a, k, m, omega, force, y):
    rhs = np.zeros(len(k))
    rhs[:len(force)] = force
    exacte = solve(k - omega**2 * m, rhs, assume_a="sym")
    candidate = np.empty(len(k))
    candidate[:len(force)] = y
    debut = len(force)
    for r, e in zip(a.sous_structures, a.applications):
        local = r.reconstruire(omega, e @ y)
        candidate[debut:debut + len(r.qr.i)] = local[r.qr.i]
        debut += len(r.qr.i)
    return exacte, candidate


class ControleChampTests(unittest.TestCase):
    def test_borne_relative_vise_la_solution_exacte_et_est_optimale(self):
        candidate = np.array([3., 4.])
        exacte = .6 * candidate
        erreur = np.linalg.norm(candidate - exacte)
        relative = erreur / np.linalg.norm(exacte)
        self.assertAlmostEqual(borne_relative(erreur, np.linalg.norm(candidate)), relative)
        self.assertGreater(relative, erreur / np.linalg.norm(candidate))
        self.assertIsNone(borne_relative(None, 5.))
        self.assertIsNone(borne_relative(5., 5.))
        self.assertIsNone(borne_relative(6., 5.))
        self.assertIsNone(borne_relative(0., 0.))
        self.assertEqual(borne_relative(0., 5.), 0.)

    def test_assemblage_masse_couplee_borne_le_champ_complet(self):
        a, k, m, restrictions = _assemblage_petit()
        controle = ControleChamp(a)
        force = np.array([.7, -.4, 1.1])
        erreurs = []
        for omega in np.linspace(.2 * a.omega_max, a.omega_max, 7):
            with self.subTest(omega=omega):
                result = controle.reponse(omega, force)
                exact, candidat = _reponse_dense_et_candidate(a, k, m, omega, force,
                                                             result["deplacement"])
                defaut = exact - candidat
                erreur = np.sqrt(defaut @ m @ defaut)
                erreurs.append(erreur)
                mesure = result["masse"]
                self.assertIsNotNone(mesure["borne_absolue"])
                self.assertLessEqual(erreur, mesure["borne_absolue"] + FLOTTANT)
                self.assertLessEqual(np.linalg.norm(defaut[:3]), result["borne_ports"] + FLOTTANT)
                self.assertAlmostEqual(mesure["norme_candidate"],
                                       np.sqrt(candidat @ m @ candidat), delta=FLOTTANT)
                if mesure["borne_relative"] is not None:
                    self.assertLessEqual(erreur / np.sqrt(exact @ m @ exact),
                                         mesure["borne_relative"] + FLOTTANT)
                for r, e, restriction in zip(a.sous_structures, a.applications, restrictions):
                    assert_allclose(r.reconstruire(omega, e @ result["deplacement"]),
                                    restriction @ candidat, atol=FLOTTANT, rtol=FLOTTANT)
                self.assertIs(result["certification_machine"], False)
        self.assertGreater(max(erreurs), 1e-6)
        # Ce cas échoue si les blocs masse intérieur/interface sont oubliés.
        self.assertGreater(np.linalg.norm(m[:3, 3:]), .1)
        mauvaise_masse = m.copy()
        mauvaise_masse[:3, 3:] = 0.
        mauvaise_masse[3:, :3] = 0.
        rhs = np.r_[force, np.zeros(len(k) - 3)]
        incorrect = solve(k - a.omega_max**2 * mauvaise_masse, rhs, assume_a="sym")
        self.assertGreater(np.linalg.norm(incorrect - exact), 1e-3)

    def test_deformation_energetique_bornee_dans_D_original(self):
        a, k, m, restrictions = _assemblage_petit()
        controle = ControleChamp(a)
        force = np.array([-.2, 1.2, .6])
        for omega in (0., .55 * a.omega_max, a.omega_max):
            with self.subTest(omega=omega):
                result = controle.reponse(omega, force)
                exact, candidat = _reponse_dense_et_candidate(a, k, m, omega, force,
                                                             result["deplacement"])
                defaut = exact - candidat
                # Énergie élémentaire : aucun Schur réduit dans cet oracle.
                carre = defaut[:3] @ a.k_externe @ defaut[:3]
                for r, restriction in zip(a.sous_structures, restrictions):
                    deformation = r.qr.d @ (restriction @ defaut)
                    carre += deformation @ deformation
                erreur = np.sqrt(carre)
                mesure = result["deformation"]
                self.assertAlmostEqual(erreur, np.sqrt(defaut @ k @ defaut), delta=FLOTTANT)
                self.assertLessEqual(erreur, mesure["borne_absolue"] + FLOTTANT)
                self.assertAlmostEqual(mesure["norme_candidate"],
                                       np.sqrt(candidat @ k @ candidat), delta=FLOTTANT)
                if mesure["borne_relative"] is not None:
                    self.assertLessEqual(erreur / np.sqrt(exact @ k @ exact),
                                         mesure["borne_relative"] + FLOTTANT)

    def test_observable_primal_dual_et_reste_exact_a_port_fixe(self):
        data = _petit(421, ni=9, p=2)
        d, m, i, s, _, lower, omega = data
        r = _releve(data, tolerance=1e-12, max_blocs=1, max_directions=len(i))
        port = np.array([.8, -.3])
        ell = np.random.default_rng(422).standard_normal(len(i))
        a = d.T @ d - omega**2 * m
        aii = a[np.ix_(i, i)]
        exact = -solve(aii, a[np.ix_(i, s)] @ r.qr.w @ port, assume_a="pos")
        candidat = r.reconstruire(omega, port)[i]
        dual_exact = solve(aii, ell, assume_a="pos")
        dual = .4 * dual_exact + .02 * np.arange(len(i))
        result = observable_locale(r, omega, port, ell, dual)
        residu = aii @ (candidat - exact)
        residu_dual = ell - aii @ dual
        reste = -residu_dual @ solve(aii, residu, assume_a="pos")
        self.assertGreater(abs(ell @ (exact - candidat)), 1e-7)
        self.assertAlmostEqual(result["valeur_candidate"], ell @ candidat, delta=FLOTTANT)
        self.assertAlmostEqual(ell @ exact - result["valeur_corrigee"], reste, delta=FLOTTANT)
        self.assertLessEqual(abs(ell @ (exact - candidat)),
                             result["borne_sans_correction"] + FLOTTANT)
        self.assertLessEqual(abs(reste), result["borne_apres_correction"] + FLOTTANT)
        self.assertGreater(lower - omega**2, 0.)
        self.assertIs(result["certification_machine"], False)

    def test_dual_exact_corrige_observable_sans_enrichir_le_primal(self):
        data = _petit(423, ni=8, p=1)
        d, m, i, s, _, _, omega = data
        r = _releve(data, tolerance=1e-12, max_blocs=1)
        port = np.array([.9])
        ell = np.linspace(-1., 1., len(i))
        a = d.T @ d - omega**2 * m
        dual = solve(a[np.ix_(i, i)], ell, assume_a="pos")
        exact = -solve(a[np.ix_(i, i)], a[np.ix_(i, s)] @ r.qr.w @ port, assume_a="pos")
        dimension = r.taille_interieure
        result = observable_locale(r, omega, port, ell, dual)
        self.assertGreater(abs(result["valeur_candidate"] - ell @ exact), 1e-7)
        self.assertAlmostEqual(result["valeur_corrigee"], ell @ exact, delta=FLOTTANT)
        self.assertLess(result["borne_apres_correction"], 1e-12 * result["borne_sans_correction"])
        self.assertEqual(r.taille_interieure, dimension)

    def test_marge_globale_insuffisante_refuse_meme_interieur_coercif(self):
        r = _releve(_petit(424, ni=8, p=1), tolerance=1e-12, max_blocs=1)
        base = AssemblagePorts([r], [np.ones((1, 1))])
        omega = .8 * r.omega_max
        cible = min(base.borne_uniforme / 4., base.schur(omega)[0, 0] / 4.)
        self.assertGreater(cible, 0.)
        masse = (base.schur(omega)[0, 0] - cible) / omega**2
        a = AssemblagePorts([r], [np.ones((1, 1))], m_externe=np.array([[masse]]))
        result = ControleChamp(a).reponse(omega, np.ones(1))
        self.assertGreater(r.lambda_min - omega**2, 0.)
        self.assertLess(result["marge"], 0.)
        self.assertIsNone(result["borne_ports"])
        for nom in ("masse", "deformation"):
            self.assertIsNone(result[nom]["borne_absolue"])
            self.assertIsNone(result[nom]["borne_relative"])
        adapte = adapter_champ(a, [omega], np.ones(1), max_etapes=0)
        self.assertEqual(adapte["statut"], "budget_etapes")
        self.assertFalse(adapte["historique"][-1]["accepte"])

    def test_controle_perime_refuse_et_nouveau_controle_recalcule_les_bornes(self):
        a, _, _, _ = _assemblage_petit()
        ancien = ControleChamp(a)
        borne = ancien.borne_schur
        self.assertTrue(a.sous_structures[0].enrichir())
        with self.assertRaisesRegex(ValueError, "périmé"):
            ancien.reponse(.8 * a.omega_max, np.ones(3))
        nouveau = ControleChamp(a)
        self.assertLess(nouveau.borne_schur, borne)
        result = nouveau.reponse(.8 * a.omega_max, np.ones(3))
        self.assertEqual(result["borne_schur"], nouveau.borne_schur)
        self.assertIsNotNone(result["masse"]["borne_absolue"])

    def test_enrichissement_reutilise_facteur_et_conserve_identite_energetique(self):
        for methode in ("qr", "lu_energie"):
            with self.subTest(methode=methode):
                data = _petit(425, ni=10, p=2)
                d, m, i, s, _, _, omega = data
                r = _releve(data, methode=methode, tolerance=1e-12, max_blocs=1,
                            max_directions=len(i))
                facteur, original, dimension = r.qr, r.qr.d, r.taille_interieure
                with patch("ports_releves.CondensationEnergie", side_effect=AssertionError("refactorisation")), \
                     patch("condensation_energie_lu.CondensationEnergieLU", side_effect=AssertionError("refactorisation")):
                    self.assertTrue(r.enrichir())
                self.assertIs(r.qr, facteur)
                self.assertIs(r.qr.d, original)
                self.assertIs(r.interieur.facteur_energie, facteur)
                self.assertGreater(r.taille_interieure, dimension)
                v = r.interieur.w
                dv = d[:, i] @ v
                assert_allclose(facteur.gram(v, v), dv.T @ dv, atol=FLOTTANT, rtol=FLOTTANT)
                assert_allclose(dv.T @ dv, np.eye(v.shape[1]), atol=FLOTTANT, rtol=FLOTTANT)
                x = r.reconstruire(omega, np.eye(len(s)))
                energie = (d @ x).T @ (d @ x) - omega**2 * (x.T @ m @ x)
                a = d.T @ d - omega**2 * m
                aii, b = a[np.ix_(i, i)], a[np.ix_(i, s)] @ r.qr.w
                exact = r.qr.w.T @ a[np.ix_(s, s)] @ r.qr.w - b.T @ solve(aii, b, assume_a="pos")
                residu = aii @ x[i] + b
                assert_allclose(energie - exact, residu.T @ solve(aii, residu, assume_a="pos"),
                                atol=FLOTTANT, rtol=FLOTTANT)
                self.assertLessEqual(np.linalg.norm(r.schur(omega) - exact, 2),
                                     r.borne_uniforme + FLOTTANT)
                self.assertEqual(r.borne_uniforme, r.audit_bande["majorant_total"])

    def test_adaptation_utilise_criteres_physiques_frais_et_actualise_assemblage(self):
        a, k, m, _ = _assemblage_petit()
        force = np.array([.8, -.5, .3])
        frequences = [0., .6 * a.omega_max, a.omega_max]
        dimensions = sum(r.taille_interieure for r in a.sous_structures)
        result = adapter_champ(a, frequences, force, tolerance_relative=1e-7, max_etapes=8)
        self.assertEqual(result["statut"], "tolerance_aux_frequences_demandees")
        self.assertFalse(result["historique"][0]["accepte"])
        self.assertTrue(result["historique"][-1]["accepte"])
        self.assertGreater(result["historique"][-1]["directions"], dimensions)
        self.assertEqual(a.borne_uniforme, result["controle"].borne_schur)
        assert_allclose(a.erreur_matrice, result["controle"].erreur_schur, atol=0., rtol=0.)
        for omega, rep in zip(frequences, result["reponses"]):
            exact, candidat = _reponse_dense_et_candidate(a, k, m, omega, force, rep["deplacement"])
            erreur = exact - candidat
            for nom, gram in (("masse", m), ("deformation", k)):
                self.assertLessEqual(rep[nom]["borne_relative"], 1e-7)
                relative = np.sqrt(erreur @ gram @ erreur / (exact @ gram @ exact))
                self.assertLessEqual(relative, rep[nom]["borne_relative"] + FLOTTANT)
        self.assertIs(result["certification_machine"], False)

    def test_budget_directions_ne_se_fait_pas_passer_pour_une_tolerance(self):
        r = _releve(_chaine(24), tolerance=1e-12, max_blocs=1, max_directions=1)
        a = AssemblagePorts([r], [np.ones((1, 1))])
        dimension = r.taille_interieure
        result = adapter_champ(a, [.3 * a.omega_max], np.ones(1),
                               tolerance_relative=1e-14, max_etapes=3)
        self.assertEqual(result["statut"], "stagnation_ou_budget_directions")
        self.assertFalse(result["historique"][-1]["accepte"])
        self.assertEqual(r.taille_interieure, dimension)
        self.assertFalse(r.enrichir())
        self.assertEqual(r.interieur.statut, "budget_directions")

    def test_force_nulle_et_energie_externe_indefinie_ne_creent_pas_de_garantie_relative(self):
        a, _, _, _ = _assemblage_petit()
        result = ControleChamp(a).reponse(.5 * a.omega_max, np.zeros(3))
        for nom in ("masse", "deformation"):
            self.assertEqual(result[nom]["norme_candidate"], 0.)
            self.assertEqual(result[nom]["borne_absolue"], 0.)
            self.assertIsNone(result[nom]["borne_relative"])
        for nom in ("m_externe", "k_externe"):
            with self.subTest(nom=nom):
                r = a.sous_structures[0]
                options = {nom: np.diag([1., -1e-6])}
                indefini = AssemblagePorts([r], [np.eye(2)], **options)
                with self.assertRaisesRegex(ValueError, "non positive"):
                    ControleChamp(indefini)


if __name__ == "__main__":
    unittest.main(verbosity=2)
