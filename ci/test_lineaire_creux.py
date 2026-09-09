"""Contre-épreuves exactes des certificats par facteurs creux."""
from fractions import Fraction as F
from itertools import product
import unittest
from unittest.mock import patch

from prototype_lineaire_creux import certifier, haut


def dense_exact(a, b):
    """Oracle d'élimination rationnelle, uniquement pour petits tests."""
    n = len(b)
    m = [[F(a[i].get(j, 0)) for j in range(n)]+[F(b[i])] for i in range(n)]
    for j in range(n):
        k = next(i for i in range(j, n) if m[i][j])
        m[j], m[k] = m[k], m[j]
        d = m[j][j]
        m[j] = [v/d for v in m[j]]
        for i in range(n):
            if i != j:
                d = m[i][j]
                m[i] = [v-d*w for v, w in zip(m[i], m[j], strict=True)]
    return [l[-1] for l in m]


class Creux(unittest.TestCase):
    def test_arrondi_relatif_sans_plancher_absolu(self):
        for e in (-1000, -200, -20, 0, 20, 200, 1000):
            for q in (F(1, 3), F(7, 5), F(1)):
                x = q*(F(2)**e)
                h = haut(x)
                self.assertLessEqual(x, h)
                self.assertLessEqual(h-x, x/F(2**127))
        with self.assertRaises(ValueError): haut(-1)

    def test_changement_unites_extreme(self):
        for facteur in (F(10)**-200, F(1), F(10)**200):
            a = [{0: 2*facteur, 1: facteur}, {0: facteur, 1: 3*facteur}]
            l = [{0: 1}, {0: F(1, 2), 1: 1}]
            u = [{0: 2*facteur, 1: facteur}, {1: F(5, 2)*facteur}]
            r = certifier(a, [facteur]*2, [F(1, 5)]*2, l, u)
            self.assertEqual(r['contraction'], 0)
            for borne in r['bornes_composantes']: self.assertLess(borne, F(1, 2))

    def test_solution_exacte_et_residu_non_nul(self):
        a = [{0: 4, 1: 1}, {0: 2, 1: 3}]
        l, u = [{0: 1}, {0: F(1, 2), 1: 1}], [{0: 4, 1: 1}, {1: F(5, 2)}]
        b = [1, 2]
        exact = dense_exact(a, b)
        r = certifier(a, b, exact, l, u)
        self.assertEqual(r['bornes_composantes'], [0, 0])
        x = [F(1, 8), F(3, 5)]
        r = certifier(a, b, x, l, u)
        for i in range(2): self.assertLessEqual(abs(x[i]-exact[i]), r['bornes_composantes'][i])

    def test_famille_complete_aux_sommets_contre_oracle(self):
        a = [{0: 4, 1: 1}, {0: 2, 1: 3}]
        l, u = [{0: 1}, {0: F(1, 2), 1: 1}], [{0: 4, 1: 1}, {1: F(5, 2)}]
        b, x, d = [1, 2], [0, F(1, 2)], F(1, 100)
        da = [{0: d, 1: d}, {0: d, 1: d}]
        r = certifier(a, b, x, l, u, delta_a=da, delta_b=[d, d])
        for signes in product((-1, 1), repeat=6):
            aa = [{j: a[i][j]+d*signes[2*i+j] for j in range(2)} for i in range(2)]
            bb = [b[i]+d*signes[4+i] for i in range(2)]
            exact = dense_exact(aa, bb)
            for i in range(2): self.assertLessEqual(abs(exact[i]-x[i]), r['bornes_composantes'][i])

    def test_facteurs_approches_et_poids(self):
        a = [{0: 2, 1: 1}, {0: 1, 1: 3}]
        l, u = [{0: 1}, {0: F(1, 2), 1: 1}], [{0: 2, 1: 1}, {1: F(251, 100)}]
        b, x = [1, 1], [F(1, 5), F(1, 5)]
        r = certifier(a, b, x, l, u, poids=[1, 2])
        self.assertGreater(r['contraction'], 0)
        exact = dense_exact(a, b)
        for i in range(2): self.assertLessEqual(abs(exact[i]-x[i]), r['bornes_composantes'][i])

    def test_kkt_indefini(self):
        # M=I, C=(1,0) ; le pivot négatif du KKT n'est pas une anomalie.
        a = [{0: 1, 2: 1}, {1: 1}, {0: 1}]
        l = [{0: 1}, {1: 1}, {0: 1, 2: 1}]
        u = [{0: 1, 2: 1}, {1: 1}, {2: -1}]
        r = certifier(a, [0, 1, 0], [0, 1, 0], l, u,
                      delta_a=[{0: F(1, 1000)}, {}, {}])
        self.assertLess(r['contraction'], 1)
        self.assertEqual(r['bornes_composantes'], [0, 0, 0])

    def test_dimension_4096_sans_inverse(self):
        n = 4096
        a = [{i: 2} for i in range(n)]
        l = [{i: 1} for i in range(n)]
        r = certifier(a, [1]*n, [F(1, 2)]*n, l, a)
        self.assertEqual(r['produits_facteurs'], n)
        self.assertEqual(r['bornes_composantes'], [0]*n)

    def test_refus_singularite_enveloppe_et_budgets(self):
        with self.assertRaises(ValueError): certifier([{}], [0], [0], [{0: 1}], [{0: 1}])
        with self.assertRaises(ValueError): certifier([{0: 1}], [0], [0], [{0: 1}], [{0: 1}], delta_a=[{0: 1}])
        with self.assertRaises(ValueError): certifier([{0: 1}], [0], [0], [{0: 0}], [{0: 1}])
        with self.assertRaises(ValueError): certifier([{0: 1}], [0], [0], [{0: 1}], [{0: 1}], poids=[0])
        with patch('prototype_lineaire_creux.PRODUITS_MAX', 0):
            with self.assertRaisesRegex(ValueError, 'budget'): certifier([{0: 1}], [0], [0], [{0: 1}], [{0: 1}])


if __name__ == '__main__': unittest.main()
