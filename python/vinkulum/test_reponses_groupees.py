"""Réponses partagées : physique, annulations, invalidation et indépendance."""
import copy
from concurrent.futures import ThreadPoolExecutor
import gc
from threading import Event
import unittest
from unittest.mock import patch
import weakref

import numpy as np
from numpy.testing import assert_allclose

from .reduction_ports import ReductionMaterielle, AssemblagePorts, ControleChamp
from .test_reduction_ports import _console, _reduire, _exacte, _bande, FLOTTANT


class ReponsesGroupeesTests(unittest.TestCase):
    def test_deux_frequences_concurrentes_gardent_leur_propre_etat(self):
        r, _, _, _ = _reduire(_console(nb=4))
        force = np.array([.2, 1., -.4, .03, .05, -.08])
        bas = r.reponse(0., force)
        haut = r.reponse(r.omega_max, force)
        self.assertGreater(np.linalg.norm(bas["champ_physique"]-haut["champ_physique"]), 1e-7)
        r._frequence_reponses = None
        pret, liberer = Event(), Event()

        def affecter(objet, nom, valeur):
            object.__setattr__(objet, nom, valeur)
            if (objet is r and nom == "_frequence_reponses"
                    and valeur is not None and valeur.omega == 0.):
                # Suspend A après publication de son état ; B remplace le
                # cache et termine avant le retour de A. Aucun aléa de
                # temporisation ne décide si la course est exercée.
                pret.set()
                if not liberer.wait(10):
                    raise AssertionError("entrelacement non libéré")

        with patch.object(ReductionMaterielle, "__setattr__", affecter):
            with ThreadPoolExecutor(max_workers=2) as pool:
                a = pool.submit(r.reponse, 0., force)
                try:
                    self.assertTrue(pret.wait(10))
                    b = pool.submit(r.reponse, r.omega_max, force).result(timeout=10)
                finally:
                    liberer.set()
                a = a.result(timeout=10)
        assert_allclose(a["champ_physique"], bas["champ_physique"], atol=0., rtol=0.)
        assert_allclose(b["champ_physique"], haut["champ_physique"], atol=0., rtol=0.)
        self.assertEqual(a["deformation"], bas["deformation"])
        self.assertEqual(b["masse"], haut["masse"])

    def test_masse_couplee_charges_combinees_et_oracle_complet(self):
        rng = np.random.default_rng(815)
        d = np.vstack((2*np.eye(9), .3*rng.normal(size=(12, 9))))
        z = rng.normal(size=(11, 9))
        m = np.eye(9)+z.T@z
        i, s = np.array([0, 1, 3, 4, 5, 6, 8]), np.array([7, 2])
        lower, omega = _bande(d, m, i)
        r = ReductionMaterielle(d, m, i, s, [[2., .3], [.3, .7]], omega,
                                lambda_min=lower, tolerance=1e-12, max_directions=len(i))
        forces = rng.normal(size=(2, 7))
        forces[:, 3] = 0.
        avant = forces.copy()
        resultats = r.reponses(omega, forces)
        np.testing.assert_array_equal(forces, avant)
        self.assertEqual(len(resultats), forces.shape[1])
        for j, rep in enumerate(resultats):
            reference = _exacte(d, m, s, omega, forces[:, j])
            champ = rep["champ_physique"]
            assert_allclose(champ, reference, atol=FLOTTANT, rtol=FLOTTANT)
            unique = r.reponse(omega, forces[:, j])
            assert_allclose(champ, unique["champ_physique"], atol=FLOTTANT, rtol=FLOTTANT)
            for nom, norme in (("masse", lambda x: np.sqrt(x@m@x)),
                               ("deformation", lambda x: np.linalg.norm(d@x))):
                self.assertAlmostEqual(rep[nom]["norme_candidate"], norme(champ), delta=FLOTTANT)
                self.assertLessEqual(norme(champ-reference), rep[nom]["borne_absolue"]+FLOTTANT)
            self.assertIs(rep["certification_machine"], False)
        for nom in ("masse", "deformation"):
            self.assertEqual(resultats[3][nom]["norme_candidate"], 0.)
            self.assertEqual(resultats[3][nom]["borne_absolue"], 0.)
            self.assertIsNone(resultats[3][nom]["borne_relative"])

    def test_reutilisation_pour_charges_successives_et_memoire_une_frequence(self):
        r, d, m, s = _reduire(_console(nb=4))
        forces = np.eye(6)
        with patch.object(r, "reconstruire", wraps=r.reconstruire) as construire:
            a = r.reponse(0., forces[:, 0])
            b = r.reponses(0., forces)
            c = r.reponse(0., forces[:, 1])
            self.assertEqual(construire.call_count, 1)
            assert_allclose(a["champ_physique"], b[0]["champ_physique"], atol=FLOTTANT, rtol=FLOTTANT)
            assert_allclose(c["champ_physique"], b[1]["champ_physique"], atol=FLOTTANT, rtol=FLOTTANT)
            ancien = weakref.ref(r._frequence_reponses)
            autre = r.reponse(r.omega_max, forces[:, 0])
            self.assertEqual(construire.call_count, 2)
        gc.collect()
        self.assertIsNone(ancien(), "les anciens champs de fréquence ne doivent pas s'accumuler")
        assert_allclose(autre["champ_physique"], _exacte(d, m, s, r.omega_max, forces[:, 0]),
                        atol=FLOTTANT, rtol=FLOTTANT)

    def test_resultats_et_appels_independants_des_mutations_de_sortie(self):
        r, _, _, _ = _reduire(_console(nb=3))
        forces = np.eye(6)
        resultats = r.reponses(r.omega_max, forces)
        avant = copy.deepcopy(resultats)
        for cle in ("champ_physique", "deplacement_ports_physique",
                    "coordonnees_ports_normalisees", "normalisation_ports", "coordonnees_physiques"):
            resultats[0][cle][...] = -123.
        resultats[0]["deformation"]["borne_absolue"] = -4.
        forces[...] = 0.
        nouveau = r.reponses(r.omega_max, np.eye(6))
        for j in range(6):
            assert_allclose(nouveau[j]["champ_physique"], avant[j]["champ_physique"], atol=0., rtol=0.)
            self.assertEqual(nouveau[j]["deformation"], avant[j]["deformation"])
            if j:
                assert_allclose(resultats[j]["champ_physique"], avant[j]["champ_physique"], atol=0., rtol=0.)

    def test_enrichissement_invalide_champs_et_bornes_sans_refactoriser(self):
        r, _, _, _ = _reduire(_console(nb=4), lambda_min=None, max_blocs=1)
        forces = np.eye(6)
        avant = r.reponses(r.omega_max, forces)
        ancien = r._frequence_reponses
        facteur, certificat = r.qr, r.certificat_spectral
        self.assertTrue(r.enrichir())
        with self.assertRaisesRegex(ValueError, "périmé"):
            ancien._evaluer(forces)
        apres = r.reponses(r.omega_max, forces)
        self.assertIsNot(ancien, r._frequence_reponses)
        self.assertIs(r.qr, facteur)
        self.assertIs(r.certificat_spectral, certificat)
        self.assertLess(apres[0]["borne_schur"], avant[0]["borne_schur"])
        for j, rep in enumerate(apres):
            unique = r.reponse(r.omega_max, forces[:, j])
            assert_allclose(rep["champ_physique"], unique["champ_physique"], atol=FLOTTANT, rtol=FLOTTANT)

    def test_arguments_invalides_ne_polluent_pas_le_cache(self):
        r, _, _, _ = _reduire(_console(nb=2))
        avant = r.reponses(0., np.eye(6))
        etat = r._frequence_reponses
        for f in (np.ones(6), np.empty((6, 0)), np.ones((5, 6)),
                  np.full((6, 2), np.nan), np.full((6, 2), np.inf), np.eye(6).astype(complex)):
            with self.assertRaises(ValueError):
                r.reponses(0., f)
        for w in (-1., np.nan, np.inf, 1.01*r.omega_max, 0j, np.array([0.])):
            with self.assertRaises(ValueError):
                r.reponses(w, np.eye(6))
        self.assertIs(r._frequence_reponses, etat)
        apres = r.reponses(0., np.eye(6))
        for a, b in zip(avant, apres):
            assert_allclose(a["champ_physique"], b["champ_physique"], atol=0., rtol=0.)

    def test_petites_deformations_annulantes_ne_passent_pas_par_un_gram(self):
        eps = 2.**-30
        metrique = np.array([[1., -1/eps], [-1/eps, 2/eps**2]])
        r = ReductionMaterielle(np.eye(3), np.eye(3), [0], [1, 2], metrique,
                                .6, lambda_min=.9)
        a = AssemblagePorts([r], [r.qr.w], m_externe=4*np.eye(2))
        w = .5
        y = np.array([1., -1.])
        forces = np.column_stack((a.schur(w)@y, -a.schur(w)@y, np.zeros(2)))
        resultats = ControleChamp(a).reponses(w, forces)
        x = r.reconstruire(w, np.eye(2))
        gram = x.T@x
        self.assertEqual(float(y@gram@y), 0., "le témoin doit perdre sa petite direction dans le Gram")
        for rep in resultats[:2]:
            physique = r.reconstruire(w, rep["deplacement"])
            norme = np.linalg.norm(r.qr.d@physique)
            self.assertGreater(rep["deformation"]["norme_candidate"], .9*eps)
            self.assertAlmostEqual(rep["deformation"]["norme_candidate"], norme, delta=1e-15)
        self.assertEqual(resultats[2]["deformation"]["norme_candidate"], 0.)

    def test_assemblage_modifie_requiert_un_controle_frais(self):
        r, _, _, _ = _reduire(_console(nb=2))
        for cible in ("k_externe", "m_externe", "applications"):
            a = AssemblagePorts([r], [r.qr.w])
            c = ControleChamp(a)
            c.reponses(0., np.eye(6))
            if cible == "applications":
                a.applications[0][0, 0] *= 2.
            else:
                getattr(a, cible)[0, 0] = 1.
            with self.assertRaisesRegex(ValueError, "périmé"):
                c.reponses(0., np.eye(6))
            self.assertEqual(len(ControleChamp(a).reponses(0., np.eye(6))), 6)


if __name__ == "__main__":
    unittest.main()
