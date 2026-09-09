"""Matrices exportées par le noyau installé, pas seulement archives."""
from fractions import Fraction as F
import unittest

from qualifie_lineaire_creux import chaine, proposition, lignes
from prototype_lineaire_creux import certifier
from test_lineaire_creux import dense_exact


class Natif(unittest.TestCase):
    def test_oracle_rationnel_et_permutations(self):
        a, b = chaine(4)
        e, pr, pc, _ = proposition(a, b)
        preuve = certifier(**e)
        oracle = dense_exact(lignes(a), list(map(F, b)))
        for i, j in enumerate(pc):
            self.assertLessEqual(abs(oracle[j]-e['x'][i]), preuve['bornes_composantes'][i])
        self.assertEqual(oracle[:24], [0]*24)

    def test_au_dela_du_budget_dense(self):
        a, b = chaine(64)
        e, _, _, _ = proposition(a, b)
        preuve = certifier(**e)
        self.assertEqual(len(a.indptr)-1, 704)
        self.assertLess(preuve['contraction'], F(1, 10**8))
        self.assertLess(max(preuve['bornes_composantes']), F(1, 10**8))


if __name__ == '__main__': unittest.main()
