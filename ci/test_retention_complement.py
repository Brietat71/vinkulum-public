"""Composition indépendante : oracle Fraction, KKT, contraintes et singularités.

Aucun test n'affirme un gain de temps ou une certification machine du témoin.
La base Z n'apparaît que dans les petites contre-épreuves exactes.
"""
from fractions import Fraction as F
from itertools import combinations
import unittest

import numpy as np

from retention_complement import TemoinRetention
from test_trace_complement import (
    matrice, transpose, produit, difference, echelle, identite, inverse,
    determinant, gram, restreindre, inverse_contrainte,
)


def tableau(a):
    return np.asarray(a, dtype=float)


class RetentionComplement(unittest.TestCase):
    def donnees(self):
        d = matrice([[2, 1, 0, 1], [0, 3, 1, 0],
                     [0, 0, 2, 1], [0, 0, 0, 4]])
        m = gram(matrice([[2, 1, 0, 0], [0, 1, 1, 0],
                          [0, 0, 2, 1], [0, 0, 0, 1]]))
        i, s = [2, 0, 1], [3]
        b, phi = [[1], [2], [-1]], [[1], [1], [0]]
        psi = [[.5], [-.25], [.75]]
        args = (tableau(d), tableau(m), i, s, b, phi, .01, .05)
        return d, m, args, dict(psi=psi, normalisation=[[2]])

    def assert_psd(self, a):
        self.assertEqual(a, transpose(a))
        for taille in range(1, len(a)+1):
            for ids in combinations(range(len(a)), taille):
                self.assertGreaterEqual(
                    determinant([[a[i][j] for j in ids] for i in ids]), 0)

    def test_masse_couplee_phi_non_modal_oracle_fraction(self):
        d, m, args, options = self.donnees()
        r = TemoinRetention(*args, **options)
        omega = .05
        forces = [[1., -2., 0.]]
        calcul = r.reponses(omega, forces)
        # L'oracle construit l'énergie exacte des valeurs D, sans K flottant.
        z = F.from_float(omega)**2
        a = difference(gram(d), echelle(m, z))
        f = matrice([[0, 0, 0], [0, 0, 0], [0, 0, 0], [1, -2, 0]])
        reference = tableau(produit(inverse(a), f))
        np.testing.assert_allclose(calcul["champ"], reference, rtol=3e-13, atol=3e-14)
        np.testing.assert_allclose(calcul["residu_physique"], 0, atol=2e-13)
        np.testing.assert_allclose(calcul["champ"][r.s], r.w@calcul["coordonnees_conservees"][:r.p],
                                   rtol=0, atol=0)
        self.assertFalse(np.allclose(r.b, r.mi@r.phi))
        self.assertGreater(np.linalg.norm(r.di.T@(r.di@r.phi)-r.mi@r.phi), 1.)
        self.assertLess(calcul["defaut_contrainte"], 5e-14)
        self.assertFalse(calcul["certification_machine"])
        self.assertEqual(calcul["taille_conservee"], 2)
        self.assertEqual(calcul["taille_kkt"], 4)
        self.assertEqual(calcul["convention_b"], "valeurs binary64 stockées")

    def test_changement_de_representants_affines_invariant(self):
        _, _, args, options = self.donnees()
        reference = TemoinRetention(*args, **options).reponses(.04, [[1.]])["champ"]
        autre = list(args)
        autre[4] = 7*np.asarray(args[4])
        autre[5] = 3*np.asarray(args[5])
        candidat = TemoinRetention(*autre, psi=[[-1.], [2.], [1.]],
                                  normalisation=[[-.5]])
        resultat = candidat.reponses(.04, [[1.]])
        np.testing.assert_allclose(resultat["champ"], reference, rtol=2e-12, atol=2e-13)

    def test_pole_interieur_traverse_sans_pole_global(self):
        d = [[1., 0., 1.], [0., 2., 0.], [0., 0., 1.]]
        r = TemoinRetention(d, np.eye(3), [0, 1], [2],
                           [[1], [0]], [[1], [0]], 4., 1.5)
        self.assertEqual(np.linalg.det(r.ki.toarray()-r.mi.toarray()), 0.)
        resultat = r.reponses(1., [[1.]])
        np.testing.assert_array_equal(resultat["schur"], [[1., 1.], [1., 0.]])
        np.testing.assert_array_equal(resultat["champ"], [[1.], [0.], [0.]])
        self.assertGreater(resultat["sigma_min_bloc_conserve"], .6)
        np.testing.assert_array_equal(resultat["residu_physique"], np.zeros((3, 1)))

    def test_pole_global_refuse_malgre_bloc_ports_inversible(self):
        r = TemoinRetention(np.diag([1., 2., 2.]), np.eye(3), [0, 1], [2],
                           [[1], [0]], [[1], [0]], 4., 1.5)
        etat = r.condensation(1.)
        self.assertEqual(etat["schur"][0, 0], 3.)  # Le port seul est inversible.
        self.assertEqual(etat["sigma_min_bloc_conserve"], 0.)
        with self.assertRaisesRegex(np.linalg.LinAlgError, "bloc conservé singulier"):
            r.reponses(1., [[1.]])

    def test_reactions_ne_sont_pas_erreurs_et_bornes_exactes(self):
        k = gram(matrice([[2, 1, 0], [0, 3, 1], [0, 0, 2]]))
        m = gram(matrice([[1, 1, 0], [0, 2, 1], [0, 0, 1]]))
        b = matrice([[1], [2], [-1]])
        zbase = matrice([[-2, 1], [1, 0], [0, 1]])
        _, tau, _, _ = inverse_contrainte(k, m, b)
        lam, freq2 = 1/tau, 1/(2*tau)
        alpha = lam-freq2
        a = difference(k, echelle(m, freq2))
        sm, _, _, _ = inverse_contrainte(m, m, b)
        erreur = zbase
        residu = produit(a, erreur)
        reaction = produit(b, matrice([[3, -5]]))
        avec_reaction = difference(residu, echelle(reaction, -1))
        gr = restreindre(sm, residu)
        self.assertEqual(restreindre(sm, avec_reaction), gr)
        self.assertEqual(produit(sm, b), matrice([[0], [0], [0]]))
        self.assert_psd(difference(echelle(gr, 1/alpha**2), restreindre(m, erreur)))
        self.assert_psd(difference(echelle(gr, 1/alpha), restreindre(a, erreur)))
        self.assert_psd(difference(echelle(gr, lam/alpha**2), restreindre(k, erreur)))

    def test_resolvante_contrainte_exacte_avec_interieur_singulier(self):
        k, m = matrice([[1, 0], [0, 4]]), identite(2)
        b, zbase = matrice([[1], [0]]), matrice([[0], [1]])
        s, _, _, _ = inverse_contrainte(k, m, b)
        a = difference(k, m)
        self.assertEqual(determinant(a), 0)
        operateur = produit(inverse(difference(identite(2), produit(s, m))), s)
        reference = produit(produit(zbase, inverse(restreindre(a, zbase))), transpose(zbase))
        self.assertEqual(operateur, reference)
        self.assertEqual(operateur, matrice([[0, 0], [0, F(1, 3)]]))

    def test_couplages_de_raideur_et_de_masse_non_facultatifs(self):
        k = matrice([[2, 1], [1, 3]])
        a = difference(k, identite(2))
        exact = a[0][0]-a[0][1]**2/a[1][1]
        self.assertEqual(exact, F(1, 2))
        self.assertEqual(a[0][0], 1)  # Faux découplage du mode approché.
        m = matrice([[2, 1], [1, 2]])
        a = difference(k, echelle(m, F(1, 2)))
        exact = a[0][0]-a[0][1]**2/a[1][1]
        sans_masse_croisee = a[0][0]-k[0][1]**2/a[1][1]
        self.assertEqual(exact, F(7, 8))
        self.assertEqual(sans_masse_croisee, F(1, 2))
        # B=e1 et Phi=e1, mais M Phi=(2,1) : annuler la masse croisée serait faux.

    def test_defaut_contrainte_interdit_borne_du_complement(self):
        epsilon = F(1, 100)
        k, m = matrice([[epsilon, 0], [0, 1]]), identite(2)
        b, erreur = matrice([[1], [0]]), matrice([[1], [0]])
        _, tau, _, _ = inverse_contrainte(k, m, b)
        lam = 1/tau
        residu = produit(k, erreur)
        annonce_fausse = restreindre(inverse(m), residu)[0][0]/lam**2
        self.assertEqual(lam, 1)
        self.assertEqual(restreindre(m, erreur)[0][0], 1)
        self.assertLess(annonce_fausse, 1)
        self.assertNotEqual(produit(transpose(b), erreur), matrice([[0]]))

    def test_rejet_paire_de_rangs_individuels_insuffisante(self):
        with self.assertRaisesRegex(ValueError, "B.T Phi inversible"):
            TemoinRetention(np.eye(3), np.eye(3), [0, 1], [2],
                            [[1], [0]], [[0], [1]], 1., .5)
        with self.assertRaisesRegex(ValueError, "B de rang plein"):
            TemoinRetention(np.eye(4), np.eye(4), [0, 1, 2], [3],
                            [[1, 2], [0, 0], [0, 0]],
                            [[1, 0], [0, 1], [0, 0]], 1., .5)

    def test_bande_entrees_et_forces_invalides(self):
        _, _, args, options = self.donnees()
        r = TemoinRetention(*args, **options)
        for omega in (-1, .06, np.nan, np.inf, 1j, [0]):
            with self.subTest(omega=omega), self.assertRaises(ValueError):
                r.condensation(omega)
        for forces in ([], [1], [[np.nan]], [[1j]], np.empty((1, 0)), [[1], [2]]):
            with self.subTest(forces=str(forces)), self.assertRaises(ValueError):
                r.reponses(0., forces)
        invalide = list(args)
        invalide[7] = np.sqrt(invalide[6])
        with self.assertRaisesRegex(ValueError, "omega_max²"):
            TemoinRetention(*invalide, **options)


if __name__ == "__main__":
    unittest.main()
