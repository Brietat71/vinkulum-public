"""Contrôles numériques de racines massiques creuses par blocs dispersés."""
import unittest
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.linalg import block_diag
from scipy.sparse import csc_matrix, csr_matrix, diags, tril

from racine_masse_blocs import RacineMasseBlocs


class RacineMasseBlocsTests(unittest.TestCase):
    @staticmethod
    def masse():
        r1 = np.array([[2., .5, -.25], [0., 1.5, .3], [0., 0., .75]])
        r2 = np.array([[1., -.4], [0., 3.]])
        return block_diag(r1.T@r1, r2.T@r2, [[4.]])

    def verifier(self, m):
        racine = RacineMasseBlocs(m)
        dense = m.toarray()
        self.assertEqual(racine.r.format, "csr")
        self.assertTrue(racine.r.has_canonical_format)
        self.assertEqual(tril(racine.r, k=-1).nnz, 0)
        self.assertFalse(racine.certification_machine)
        assert_allclose((racine.r.T@racine.r).toarray(), dense, rtol=3e-15, atol=1e-15)
        rng = np.random.default_rng(2026)
        rhs = rng.normal(size=(m.shape[0], 4))
        original = rhs.copy()
        y = racine.dual(rhs)
        assert_allclose(racine.r.T@y, rhs, rtol=2e-14, atol=2e-15)
        reference = np.linalg.solve(racine.r.toarray().T, rhs)
        assert_allclose(y, reference, rtol=2e-14, atol=2e-15)
        # Dualité physique : ||R^-T f||² = f.T M^-1 f.
        dual2 = np.sum(rhs*np.linalg.solve(dense, rhs), axis=0)
        assert_allclose(np.sum(y*y, axis=0), dual2, rtol=2e-14, atol=2e-15)
        norme = np.sqrt(np.sum(rhs*(dense@rhs), axis=0))
        assert_allclose(racine.norme(rhs), norme, rtol=2e-14, atol=2e-15)
        assert_allclose(racine.norme(rhs[:, 2]), norme[2], rtol=2e-14)
        assert_allclose(racine.dual(rhs[:, 2]), y[:, 2], rtol=2e-14, atol=2e-15)
        self.assertIsInstance(racine.norme(rhs[:, 0]), float)
        assert_array_equal(rhs, original)

    def test_masse_couplee_csr_et_csc(self):
        for classe in (csr_matrix, csc_matrix):
            with self.subTest(classe=classe.__name__):
                self.verifier(classe(self.masse()))

    def test_permutation_physique_et_normes_invariantes(self):
        m = self.masse()
        p = np.array([4, 1, 5, 2, 0, 3])
        self.verifier(csr_matrix(m[np.ix_(p, p)]))
        original = RacineMasseBlocs(csr_matrix(m))
        permute = RacineMasseBlocs(csr_matrix(m[np.ix_(p, p)]))
        f = np.arange(1., 7.)
        assert_allclose(permute.norme(f[p]), original.norme(f), rtol=2e-15)
        assert_allclose(np.linalg.norm(permute.dual(f[p])),
                        np.linalg.norm(original.dual(f)), rtol=2e-15)

    def test_diagonale_sans_composantes_ni_cholesky(self):
        m = diags(np.array([1., 4., 9., 16.]), format="csc")
        with patch("racine_masse_blocs.connected_components", side_effect=AssertionError), \
                patch("racine_masse_blocs.cholesky", side_effect=AssertionError), \
                patch("racine_masse_blocs.spsolve_triangular", side_effect=AssertionError):
            r = RacineMasseBlocs(m, taille_bloc_max=1)
            assert_array_equal(r.r.diagonal(), [1., 2., 3., 4.])
            assert_array_equal(r.dual(np.ones(4)), [1., .5, 1/3, .25])
            assert_array_equal(r.dual(np.ones((4, 2))),
                               np.array([1., .5, 1/3, .25])[:, None]*np.ones((1, 2)))
            self.assertEqual(r.norme(np.zeros(4)), 0.)

    def test_pas_de_conversion_dense_globale(self):
        bloc = np.array([[3., .2, -.1], [.2, 2., .3], [-.1, .3, 1.]])
        dense = block_diag(*([bloc]*32))
        p = np.random.default_rng(27).permutation(len(dense))
        m = csr_matrix(dense[np.ix_(p, p)])
        ancienne = csr_matrix.toarray
        tailles = []

        def petit_bloc(a, *args, **kwargs):
            self.assertLessEqual(max(a.shape), 6, "conversion dense globale")
            tailles.append(a.shape)
            return ancienne(a, *args, **kwargs)

        with patch.object(csr_matrix, "toarray", petit_bloc), \
                patch.object(csc_matrix, "toarray", side_effect=AssertionError):
            r = RacineMasseBlocs(m)
            y = r.dual(np.ones((len(dense), 2)))
            normes = r.norme(y)
        self.assertEqual(tailles, [(3, 3)]*32)
        self.assertEqual(y.shape, (96, 2))
        self.assertEqual(normes.shape, (2,))
        assert_allclose(r.r.T@y, 1., rtol=2e-14, atol=2e-15)

    def test_zeros_explicites_et_copie_entree(self):
        m = csr_matrix(([4., 0., 0., 9.], ([0, 0, 1, 1], [0, 1, 0, 1])), shape=(2, 2))
        r = RacineMasseBlocs(m, taille_bloc_max=1)
        self.assertEqual(m.nnz, 4)
        self.assertEqual(r.r.nnz, 2)
        m.data[:] = 0.
        assert_array_equal(r.r.diagonal(), [2., 3.])
        assert_array_equal(r.dual(np.array([2., 3.])), [1., 1.])

    def test_budget_et_bloc_trop_grand_avant_conversion(self):
        m = csr_matrix(np.eye(7)+np.ones((7, 7)))
        with patch.object(csr_matrix, "toarray", side_effect=AssertionError):
            with self.assertRaisesRegex(ValueError, "taille de bloc"):
                RacineMasseBlocs(m)
            with self.assertRaisesRegex(ValueError, "taille de bloc"):
                RacineMasseBlocs(csr_matrix(self.masse()), taille_bloc_max=2)
        for budget in (0, -1, 7, 2.5, True):
            with self.subTest(budget=budget), self.assertRaises(ValueError):
                RacineMasseBlocs(csr_matrix(np.eye(2)), taille_bloc_max=budget)

    def test_masses_invalides(self):
        m = self.masse()
        non_sym = m.copy()
        non_sym[1, 0] = np.nextafter(non_sym[1, 0], np.inf)
        cas = [m, csr_matrix(m).tocoo(), csr_matrix(m.astype(np.float32)),
               csr_matrix(m.astype(complex)), csr_matrix(np.eye(2, dtype=np.int64)),
               csr_matrix((0, 0), dtype=float), csr_matrix(np.ones((2, 3))),
               csr_matrix([[1., 0.], [0., 0.]]), csr_matrix([[-1.]]),
               csr_matrix([[1., 2.], [2., 1.]]), csr_matrix(non_sym),
               csr_matrix([[np.nan]]), csr_matrix([[np.inf]])]
        for i, masse in enumerate(cas):
            with self.subTest(cas=i), self.assertRaises(ValueError):
                RacineMasseBlocs(masse)

    def test_stockage_non_canonique_et_drapeaux_falsifies(self):
        doublons = csr_matrix(([1., 2., 3.], [0, 0, 1], [0, 2, 3]), shape=(2, 2))
        non_trie = csr_matrix(([1., 2., 1.], [1, 0, 1], [0, 2, 3]), shape=(2, 2))
        hors_limite = csr_matrix(np.eye(2))
        hors_limite.indices[1] = 2
        pointeur = csr_matrix(np.eye(2))
        pointeur.indptr[1] = 3
        cas = [doublons, non_trie, hors_limite, pointeur]
        cas.append(csc_matrix((doublons.data.copy(), doublons.indices.copy(),
                               doublons.indptr.copy()), shape=(2, 2)))
        for masse in cas:
            masse.has_sorted_indices = True
            masse.has_canonical_format = True
            with self.subTest(format=masse.format, indices=masse.indices.tolist()), \
                    self.assertRaisesRegex(ValueError, "canonique|stockage"):
                RacineMasseBlocs(masse)

    def test_seconds_membres_invalides_et_norme_echelle_extreme(self):
        r = RacineMasseBlocs(csr_matrix(np.eye(2)))
        for rhs in (1., [1.], np.ones((2, 0)), np.ones((2, 1, 1)),
                    [1., np.nan], [np.inf, 0.], [1.+0j, 2.+0j]):
            for methode in (r.dual, r.norme):
                with self.subTest(rhs=rhs, methode=methode.__name__), self.assertRaises(ValueError):
                    methode(rhs)
        self.assertTrue(np.isfinite(r.norme(np.array([1.e200, -1.e200]))))
        assert_allclose(r.norme(np.array([1.e200, -1.e200])), np.sqrt(2.)*1.e200)
        un = RacineMasseBlocs(csr_matrix([[4.]]))
        self.assertEqual(un.norme(np.array([-3.])), 6.)


if __name__ == "__main__":
    unittest.main()
