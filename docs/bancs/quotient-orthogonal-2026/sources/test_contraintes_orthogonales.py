"""Références indépendantes du quotient local : géométrie, SVD et forces."""
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import block_diag, null_space
from scipy.sparse import csr_matrix, diags, vstack
from scipy.spatial.transform import Rotation

from contraintes_orthogonales import BudgetDepasse, QuotientOrthogonal, RangAmbigu


def chaine(n, redondante=True):
    g = diags((np.ones(n-1), -np.ones(n-1)), (0, 1), shape=(n-1, n), format='csr')
    return vstack((g, g[::2]+g[1::2]), format='csr') if redondante and n % 2 else g


class ContraintesOrthogonales(unittest.TestCase):
    def verifie(self, a, *, atol=3e-12):
        a = np.asarray(a)
        q = QuotientOrthogonal(a)
        # SVD complète seulement dans l'oracle des petits cas.
        z = null_space(q.g.toarray(), rcond=q.seuil)
        self.assertEqual(q.rang, a.shape[1]-z.shape[1])
        rng = np.random.default_rng(469)
        v = rng.normal(size=(a.shape[1], 4))
        projected = q.projette(v)
        np.testing.assert_allclose(projected, z@(z.T@v), atol=atol)
        np.testing.assert_allclose(q.projette(projected), projected, atol=atol)
        np.testing.assert_allclose(q.orthogonal(q.orthogonal(v), transpose=True), v, atol=atol)
        np.testing.assert_allclose(q.injecte(q.restreint(v)), projected, atol=atol)
        np.testing.assert_allclose(q.g@projected, 0., atol=atol)
        np.testing.assert_allclose(a.T@q.reactions(v), v-projected, atol=atol)
        # Reconstruction de toutes les colonnes et toutes les lignes.
        transformed = q.g@q.orthogonal(np.eye(q.n))
        np.testing.assert_allclose(transformed[:, q.pivots], q.r.toarray(), atol=atol)
        np.testing.assert_allclose(transformed[:, q.libres], 0., atol=atol)
        np.testing.assert_allclose(q.projette(v[:, 0]), projected[:, 0], atol=atol)
        np.testing.assert_allclose(q.reactions(v[:, 0]), q.reactions(v)[:, 0], atol=atol)
        return q

    def test_chaine_redondances_generales_et_permutations(self):
        a = chaine(81).toarray()
        rng = np.random.default_rng(813)
        for _ in range(3):
            self.verifie(a[rng.permutation(len(a))][:, rng.permutation(a.shape[1])])

    def test_noyau_connu_rotations_locales(self):
        rng = np.random.default_rng(621)
        rotations = Rotation.from_rotvec(rng.normal(size=(9, 3))).as_matrix()
        r = block_diag(*rotations)
        a = np.kron(chaine(9).toarray(), np.eye(3))@r
        q = self.verifie(a)
        self.assertEqual(q.rang, 24)
        v = rng.normal(size=(27, 2))
        w = (r@v).reshape(9, 3, 2)
        ref = r.T@np.broadcast_to(w.mean(axis=0), w.shape).reshape(27, 2)
        np.testing.assert_allclose(q.projette(v), ref, atol=3e-13)

    def test_sous_espace_aleatoire_et_dependances_denses(self):
        rng = np.random.default_rng(945)
        z, _ = np.linalg.qr(rng.normal(size=(21, 21)))
        b = rng.normal(size=(30, 13))@z[:, :13].T
        self.verifie(b)

    def test_changements_echelles_des_equations(self):
        b = chaine(11).toarray()
        scales = (-1.)**np.arange(len(b))*np.geomspace(1e-120, 1e120, len(b))
        self.verifie(b*scales[:, None])

    def test_metrique_de_masse_projection_ponderee(self):
        mass = np.geomspace(1e-3, 1e3, 15)
        w = 1/np.sqrt(mass)
        q = QuotientOrthogonal(chaine(15)@diags(w))
        v = np.random.default_rng(721).normal(size=15)
        p = w*q.projette(v/w)
        np.testing.assert_allclose(p, np.full(15, np.dot(mass, v)/sum(mass)), atol=2e-12)

    def test_cascade_native_et_oracle_cinematique(self):
        from mesure_quotient_orthogonal import modele, reference
        for family in ('cascade', 'cascade_tournee'):
            g, rotation = modele(family, 3)
            q = self.verifie(g.toarray())
            v = np.random.default_rng(936).normal(size=(g.shape[1], 3))
            ref, rank = reference(family, 3, rotation, v)
            self.assertEqual(q.rang, rank)
            np.testing.assert_allclose(q.projette(v), ref, atol=1e-12)

    def test_vide_zero_rang_plein_et_composantes(self):
        for b in (np.zeros((0, 7)), np.zeros((4, 7)), np.eye(7),
                  np.zeros((4, 0)), block_diag(chaine(7).toarray(), np.eye(3), np.zeros((2, 4)))):
            self.verifie(b)

    def test_grande_chaine_sans_conversion_dense(self):
        n = 2049
        g = chaine(n)
        v = np.random.default_rng(825).normal(size=(n, 3))
        with patch.object(csr_matrix, 'toarray', side_effect=AssertionError('densification')):
            q = QuotientOrthogonal(g)
            p = q.projette(v)
            np.testing.assert_allclose(p, np.broadcast_to(v.mean(axis=0), v.shape), atol=1e-12)
            np.testing.assert_allclose(g.T@q.reactions(v), v-p, atol=1e-11)
        self.assertEqual(q.rang, n-1)
        self.assertEqual(q.diagnostic['support_max'], 2)
        self.assertLess(q.diagnostic['coefficients_reflexions'], 2*n)
        self.assertLess(q.diagnostic['coefficients_transformes'], 12*n)

    def test_reactions_pas_de_norme_minimale(self):
        a = np.array([[1., 0.], [2., 0.]])
        q = self.verifie(a)
        force = np.array([1., 0.])
        eta = q.reactions(force)
        ref = np.linalg.lstsq(a.T, force, rcond=None)[0]
        self.assertGreater(np.linalg.norm(eta), 1.5*np.linalg.norm(ref))
        np.testing.assert_allclose(a.T@eta, a.T@ref, atol=1e-14)

    def test_refus_au_seuil_et_budget_de_remplissage(self):
        for delta in (1e-9/3, 1e-9, 3e-9):
            with self.assertRaises(RangAmbigu):
                QuotientOrthogonal([[1., 0.], [1., delta]], seuil=1e-9)
        self.assertEqual(QuotientOrthogonal([[1., 0.], [1., 1e-11]], seuil=1e-9).rang, 1)
        self.assertEqual(QuotientOrthogonal([[1., 0.], [1., 1e-7]], seuil=1e-9).rang, 2)
        with self.assertRaises(BudgetDepasse):
            QuotientOrthogonal(np.ones((20, 20)), budget=50)
        for b in ([[np.nan]], [[np.inf]], [[1.+1j]]):
            with self.assertRaises(ValueError):
                QuotientOrthogonal(b)

    def test_petit_residu_ne_certifie_pas_la_mobilite(self):
        # det(G)=delta ≠ 0 : le véritable noyau est {0}, même pour delta
        # minuscule. La troncature numérique crée une direction entière.
        q = QuotientOrthogonal([[1., 0.], [1., 1e-14]], seuil=1e-10)
        p = q.projette([0., 1.])
        self.assertEqual(q.rang, 1)
        self.assertLess(np.linalg.norm(q.g@p), 2e-14)
        self.assertEqual(np.linalg.norm(p), 1.)


if __name__ == '__main__':
    unittest.main()
