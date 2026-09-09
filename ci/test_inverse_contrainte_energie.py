"""Oracles rationnels et contre-épreuves du projecteur en énergie."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, diags

import inverse_contrainte_energie as module
from inverse_contrainte_energie import InverseContrainteEnergie
from test_trace_complement import (
    gram, identite, inverse, inverse_contrainte, produit, transpose)


def rationnel(a):
    return [[F(float(v)) for v in row] for row in a]


def flottant(a):
    return np.array([[float(v) for v in row] for row in a])


class InverseContrainteEnEnergie(unittest.TestCase):
    def donnees(self):
        r = np.array([[2., 1., 0., 1.], [0., 3., 1., 0.],
                      [0., 0., 2., 1.], [0., 0., 0., 4.]])
        e = np.array([.5, 2., .25, 4.])
        b = np.array([[1., 0.], [0., 1.], [2., -1.], [1., 3.]])
        k = gram([[F(float(v))/F(float(e[j])) for j, v in enumerate(row)] for row in r])
        return r, e, b, k

    def test_oblique_compare_inverse_et_relevement_rationnels(self):
        r, e, b, k = self.donnees()
        f = np.array([[1., -2., 0.], [3., 1., 0.], [-1., 4., 0.], [2., -3., 0.]])
        s, _, _, _ = inverse_contrainte(k, identite(4), rationnel(b))
        u = produit(inverse(k), rationnel(b))
        c = produit(u, inverse(produit(transpose(rationnel(b)), u)))
        a = InverseContrainteEnergie(csr_matrix(r), e, b)
        x = a.appliquer(f)
        np.testing.assert_allclose(x, flottant(produit(s, rationnel(f))), rtol=3e-12, atol=3e-14)
        np.testing.assert_allclose(a.canonique(), flottant(c), rtol=3e-12, atol=3e-14)
        np.testing.assert_allclose(b.T@a.canonique(), np.eye(2), rtol=0, atol=3e-14)
        controle = a.diagnostic(f, x, d_i=csr_matrix(r/e))
        self.assertLess(controle["defaut_contrainte"], 3e-13)
        self.assertLess(controle["residu_Q_dual_contraint_relatif"], 3e-13)
        self.assertLess(controle["residu_D_dual_contraint_relatif"], 3e-13)
        self.assertFalse(controle["certification_machine"])

    def test_cancellation_ne_soustrait_pas_deux_reponses_geantes(self):
        r = np.diag([2.**-50, 1., 2.])
        b = np.array([[1.], [1.], [0.]])
        f = np.array([[1.], [0.], [0.]])
        k = gram(rationnel(r))
        s, _, _, _ = inverse_contrainte(k, identite(3), rationnel(b))
        exact = flottant(produit(s, rationnel(f)))
        a = InverseContrainteEnergie(r, np.ones(3), b)
        x = a.appliquer(f)
        np.testing.assert_allclose(x, exact, rtol=3e-15, atol=1e-30)
        np.testing.assert_allclose(b.T@x, np.zeros((1, 1)), rtol=0, atol=1e-15)
        # La formule en différence est mathématiquement correcte, mais ses
        # deux termes perdent ici le premier coefficient physique, proche de 1.
        ki = np.diag([2.**100, 1., .25])
        u = ki@b
        naive = ki@f-u@np.linalg.solve(b.T@u, u.T@f)
        self.assertGreater(abs(naive[0, 0]-exact[0, 0]), .5)
        self.assertLess(abs(x[0, 0]-exact[0, 0]), 1e-14)

    def test_dual_contraint_compare_gram_rationnel_sans_rotation_retour(self):
        r, e, b, k = self.donnees()
        f = np.array([[1., -2., 0.], [3., 1., 0.], [-1., 4., 0.], [2., -3., 0.]])
        s, _, _, _ = inverse_contrainte(k, identite(4), rationnel(b))
        exact = flottant(produit(transpose(rationnel(f)), produit(s, rationnel(f))))
        a = InverseContrainteEnergie(r, e, b)
        rotations = []
        original = a._rotation
        def rotation(x, transpose=False):
            rotations.append(transpose)
            return original(x, transpose=transpose)
        with patch.object(a, "_rotation", rotation), \
             patch.object(a, "_physique", side_effect=AssertionError("résolution R inutile")):
            y = a.dual_contraint(f)
        self.assertEqual(y.shape, (2, 3))
        self.assertEqual(rotations, [True])
        np.testing.assert_allclose(y.T@y, exact, rtol=3e-12, atol=3e-14)
        # La composante duale pertinente reste d'ordre 1 malgré une énergie
        # non contrainte de l'ordre de 2^100 pour cette même charge.
        petit = InverseContrainteEnergie(np.diag([2.**-50, 1., 2.]), np.ones(3),
                                        np.array([[1.], [1.], [0.]]))
        ff = np.array([[1.], [0.], [0.]])
        yy = petit.dual_contraint(ff)
        attendu = float(F(1)/(1+F(1, 2**100)))
        self.assertAlmostEqual(float(np.sum(yy*yy)), attendu, places=14)

    def test_reechelonnement_extreme_et_permutation_conservent_le_complement(self):
        r, e, b, _ = self.donnees()
        f = np.array([[1., 3.], [-2., 1.], [4., -1.], [0., 2.]])
        a = InverseContrainteEnergie(r, e, b)
        reference, c = a.appliquer(f), a.canonique()
        for exposants in ((-400, 400), (400, -400)):
            with self.subTest(exposants=exposants):
                echelles = np.array([2.**p for p in exposants])
                bs = (b*echelles[None, :])[:, ::-1]
                autre = InverseContrainteEnergie(r, e, bs)
                x = autre.appliquer(f)
                np.testing.assert_allclose(x, reference, rtol=3e-12, atol=3e-14)
                # Comparer dans les unités originales : les petits défauts
                # hors diagonale de B.T C peuvent être énormément amplifiés
                # par le rapport 2^800 des unités des contraintes.
                np.testing.assert_allclose(autre.canonique()[:, ::-1]*echelles[None, :],
                                           c, rtol=3e-12, atol=3e-14)
                controle = autre.diagnostic(f, x, d_i=r/e)
                self.assertLess(controle["defaut_contrainte_colonnes_normalisees"], 3e-13)

    def test_reactions_et_defaut_du_modele_D_sont_distingues(self):
        r, e, b, _ = self.donnees()
        a = InverseContrainteEnergie(r, e, b)
        reaction = b@np.array([[2., -1.], [-3., 4.]])
        x = a.appliquer(reaction)
        np.testing.assert_allclose(x, np.zeros_like(x), rtol=0, atol=3e-13)
        controle = a.diagnostic(reaction, np.zeros_like(x), d_i=r/e)
        self.assertLess(controle["residu_D_dual_contraint_relatif"], 3e-14)
        self.assertGreater(np.linalg.norm(reaction), 1.)
        f = np.array([[1.], [-2.], [3.], [4.]])
        x = a.appliquer(f)
        d = np.vstack((r/e, [1., -2., 3., 4.]))
        controle = a.diagnostic(f, x, d_i=d)
        self.assertLess(controle["residu_Q_dual_contraint_relatif"], 3e-13)
        self.assertGreater(controle["residu_D_dual_contraint"], .01)

    def test_chainons_256_et_3066_sans_inverse_ni_base_dense(self):
        def pas_de_grand_dense(original):
            def convertir(a, *args, **kwargs):
                if min(a.shape) > 4:
                    raise AssertionError("conversion d'une grande matrice dense interdite")
                return original(a, *args, **kwargs)
            return convertir
        qr_original = module.qr
        def verifier_qr(a, *args, **kwargs):
            self.assertEqual(kwargs.get("mode"), "raw")
            self.assertEqual(a.shape[1], 2)
            return qr_original(a, *args, **kwargs)
        for n in (256, 3066):
            with self.subTest(n=n):
                r = diags([np.ones(n), .25*np.ones(n-1)], [0, 1], format="csr")
                e = 2.**(np.arange(n)%5-2)
                b = np.zeros((n, 2)); b[0, 0] = b[1, 1] = 1.
                f = np.column_stack((np.sin(np.arange(n)/7), np.cos(np.arange(n)/11)))
                with patch.object(csr_matrix, "toarray", pas_de_grand_dense(csr_matrix.toarray)), \
                     patch.object(csc_matrix, "toarray", pas_de_grand_dense(csc_matrix.toarray)), \
                     patch.object(np.linalg, "inv", side_effect=AssertionError("inverse dense interdite")), \
                     patch.object(module, "qr", verifier_qr):
                    a = InverseContrainteEnergie(r, e, b)
                    x, c = a.appliquer(f), a.canonique()
                    controle = a.diagnostic(f, x, d_i=r@diags(1/e))
                self.assertEqual(a._reflecteurs.shape, (n, 2))
                self.assertFalse(any(isinstance(v, np.ndarray) and v.shape == (n, n)
                                     for v in vars(a).values()))
                np.testing.assert_allclose(b.T@x, np.zeros((2, 2)), rtol=0, atol=3e-13)
                np.testing.assert_allclose(b.T@c, np.eye(2), rtol=0, atol=3e-13)
                self.assertLess(controle["residu_D_dual_contraint_relatif"], 3e-13)

    def test_entrees_et_rang_numerique_invalides_refusent(self):
        r, e, b, _ = self.donnees()
        for mauvais in (np.zeros((4, 1)), np.column_stack((b[:, 0], b[:, 0])),
                        b.astype(complex), np.full_like(b, np.nan), np.ones((3, 1)), np.eye(4)):
            with self.subTest(forme=mauvais.shape):
                with self.assertRaises(ValueError):
                    InverseContrainteEnergie(r, e, mauvais)
        quasi = np.array([[1., 1.], [0., 2.**-52], [0., 0.], [0., 0.]])
        with self.assertRaisesRegex(ValueError, "rang numérique"):
            InverseContrainteEnergie(r, e, quasi)
        rr = r.copy(); rr[1, 0] = .25
        with self.assertRaisesRegex(ValueError, "triangulaire"):
            InverseContrainteEnergie(rr, e, b)
        with self.assertRaises(ValueError):
            InverseContrainteEnergie(r, np.zeros(4), b)
        a = InverseContrainteEnergie(r, e, b)
        for mauvais in (np.ones(4), np.empty((4, 0)), np.ones((3, 2)),
                        np.ones((4, 2), dtype=complex), np.full((4, 2), np.inf)):
            with self.assertRaises(ValueError):
                a.appliquer(mauvais)

    def test_entrees_et_resultats_ne_sont_pas_des_caches_mutables(self):
        r, e, b, _ = self.donnees()
        f = np.arange(8.).reshape(4, 2)
        a = InverseContrainteEnergie(r, e, b)
        attendu = a.appliquer(f)
        r[:] = 0.; e[:] = 0.; b[:] = 0.
        obtenu = a.appliquer(f)
        np.testing.assert_array_equal(obtenu, attendu)
        obtenu[:] = 1e30
        np.testing.assert_array_equal(a.appliquer(f), attendu)
        np.testing.assert_array_equal(f, np.arange(8.).reshape(4, 2))


if __name__ == "__main__":
    unittest.main()
