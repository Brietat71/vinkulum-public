"""Témoins du juge commun ; aucun appel ni code source concurrent requis."""
import unittest

import numpy as np
from scipy.sparse import diags, eye

from confronte_ports_exudyn import CIBLE, juger


class TestJugePorts(unittest.TestCase):
    def test_champs_exacts_et_couverture_des_resultats(self):
        champ = np.vstack((.2*np.eye(6), np.eye(6)))
        references = np.stack((champ, 2*champ))
        resultat = juger(eye(12, format="csr"), diags(np.arange(1., 13.)),
                         np.eye(6), references.copy(), references)
        self.assertTrue(resultat["accepte"])
        for nom in ("masse", "deformation", "port"):
            self.assertEqual(resultat["erreurs"][nom], [[0.]*6, [0.]*6])
            self.assertEqual(resultat["operateurs"][nom], [0., 0.])
            self.assertEqual(resultat["maxima"][nom], 0.)
            self.assertEqual(resultat["maxima_operateurs"][nom], 0.)

    def test_combinaison_annulante_refusee_malgre_six_colonnes_precises(self):
        reference = np.eye(6)
        epsilon = 1e-4
        reference[:, 1] = reference[:, 0]
        reference[1, 1] = epsilon
        candidat = reference.copy()
        candidat[1, 1] += .5*CIBLE
        resultat = juger(eye(6, format="csr"), eye(6, format="csr"), np.eye(6),
                         candidat[None], reference[None])
        self.assertTrue(all(v < CIBLE for v in resultat["maxima"].values()))
        self.assertFalse(resultat["accepte"])
        attendu = (candidat[1, 1]-reference[1, 1])/epsilon
        for nom in ("masse", "deformation", "port"):
            self.assertAlmostEqual(resultat["maxima_operateurs"][nom]/attendu, 1., places=10)

    def test_normalisation_port_et_transpose_cholesky(self):
        l = np.eye(6)
        l[:2, :2] = [[2., 0.], [1., 3.]]
        metrique = l@l.T
        reference = np.eye(6)
        candidat = reference.copy()
        delta = .2*CIBLE
        candidat[0, 1] += delta
        resultat = juger(eye(6, format="csr"), eye(6, format="csr"), metrique,
                         candidat[None], reference[None])
        self.assertAlmostEqual(resultat["maxima_operateurs"]["port"]/(2*delta/3), 1., places=12)
        self.assertAlmostEqual(resultat["erreurs"]["port"][0][1]/(2*delta/np.sqrt(10)), 1., places=12)

    def test_formes_nan_et_masse_hors_contrat(self):
        d, m, metrique = eye(6, format="csr"), eye(6, format="csr"), np.eye(6)
        valide = np.eye(6)[None]
        invalides = (np.eye(6), np.ones((1, 6, 5)), np.ones((0, 6, 6)),
                     np.ones((2, 6, 6)), np.full((1, 6, 6), np.nan),
                     np.full((1, 6, 6), np.inf), valide.astype(complex))
        for invalide in invalides:
            with self.assertRaises(ValueError):
                juger(d, m, metrique, invalide, valide)
        with self.assertRaises(ValueError):
            juger(d, m, metrique, valide, np.full((1, 6, 6), np.nan))
        non_diagonale = np.eye(6)
        non_diagonale[0, 1] = non_diagonale[1, 0] = .1
        with self.assertRaises(ValueError):
            juger(d, non_diagonale, metrique, valide, valide)
        with self.assertRaises(ValueError):
            juger(d, m, np.eye(5), valide, valide)

    def test_reference_singuliere_ou_numériquement_dependante_refusee(self):
        for epsilon in (0., 1e-18):
            reference = np.eye(6)
            reference[:, 1] = reference[:, 0]
            reference[1, 1] = epsilon
            with self.assertRaisesRegex(ArithmeticError, "rang inférieur à six"):
                juger(eye(6, format="csr"), eye(6, format="csr"), np.eye(6),
                      reference[None], reference[None])
        with self.assertRaisesRegex(ArithmeticError, "norme nulle"):
            juger(eye(6, format="csr"), eye(6, format="csr"), np.eye(6),
                  np.zeros((1, 6, 6)), np.zeros((1, 6, 6)))

    def test_amplitudes_disparates_et_changement_de_base_des_charges(self):
        reference = np.diag(np.geomspace(1e-12, 1e12, 6))
        candidat = (1+.1*CIBLE)*reference
        resultat = juger(eye(6, format="csr"), eye(6, format="csr"), np.eye(6),
                         candidat[None], reference[None])
        self.assertTrue(resultat["accepte"])
        rng = np.random.default_rng(45)
        reference = np.vstack((rng.normal(size=(6, 6)), np.eye(6)))
        candidat = reference+(.1*CIBLE)*rng.normal(size=(12, 6))
        changement = 2*np.eye(6)+np.triu(np.ones((6, 6)), 1)/8
        d, m = eye(12, format="csr"), diags(np.linspace(.5, 2., 12))
        a = juger(d, m, np.eye(6), candidat[None], reference[None])
        b = juger(d, m, np.eye(6), (candidat@changement)[None], (reference@changement)[None])
        for nom in ("masse", "deformation", "port"):
            self.assertAlmostEqual(a["maxima_operateurs"][nom]/b["maxima_operateurs"][nom], 1., places=8)


if __name__ == "__main__":
    unittest.main()
