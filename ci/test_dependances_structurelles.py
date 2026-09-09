"""Rang préservé par structure, et impossibilité sans ce contrat."""
from fractions import Fraction as F
from itertools import product
import unittest

from prototype_lineaire_creux import certifier_dependances_structurelles, kkt
from test_lineaire_creux import dense_exact


class Dependances(unittest.TestCase):
    def test_famille_perturbee_et_reactions(self):
        m, c, t = [{0: 1}, {1: 1}], [{0: 1}], [{0: 1}, {0: 2}, {0: -3}]
        l, u = [{0: 1}, {1: 1}, {0: 1, 2: 1}], [{0: 1, 2: 1}, {1: 1}, {2: -1}]
        eps = F(1, 100)
        preuve = certifier_dependances_structurelles(m, c, t, [0], [0, 1], [0], [0, 1, 0], l, u,
                    delta_c=[{0: eps, 1: eps}], delta_t=[{}, {0: 10}, {0: 10}],
                    delta_force=[eps, eps], delta_d=[eps])
        self.assertEqual(preuve['rang_uniforme'], 1)
        self.assertEqual(preuve['contraintes_originales'], 3)
        for signes in product((-1, 1), repeat=5):
            cc = [{0: 1+eps*signes[0], 1: eps*signes[1]}]
            f, d = [eps*signes[2], 1+eps*signes[3]], [eps*signes[4]]
            sol = dense_exact(kkt(m, cc), f+d)
            for i, cible in enumerate([0, 1, 0]):
                self.assertLessEqual(abs(sol[i]-cible), preuve['systeme_reduit']['bornes_composantes'][i])
            for i in range(2):
                reaction = cc[0][i]*sol[2]
                self.assertLessEqual(abs(reaction-preuve['reaction_proposee'][i]), preuve['bornes_reaction'][i])
            # Même lorsque les lignes dépendantes deviennent nulles ou
            # changent de signe, la ligne de base demeure exactement C'.
            for alpha in (-8, 0, 12):
                self.assertEqual(alpha*sum(cc[0][i]*sol[i] for i in range(2)), alpha*d[0])

    def test_refus_identite_perturbee_et_rang_perdu(self):
        m, c = [{0: 1}, {1: 1}], [{0: 1}]
        l, u = [{0: 1}, {1: 1}, {0: 1, 2: 1}], [{0: 1, 2: 1}, {1: 1}, {2: -1}]
        for t, dt in [([{0: 2}], None), ([{0: 1}], [{0: F(1, 100)}])]:
            with self.assertRaisesRegex(ValueError, 'identité'):
                certifier_dependances_structurelles(m, c, t, [0], [0, 1], [0], [0, 1, 0], l, u, delta_t=dt)
        with self.assertRaisesRegex(ValueError, 'contraction'):
            certifier_dependances_structurelles(m, c, [{0: 1}], [0], [0, 1], [0], [0, 1, 0], l, u,
                                                delta_c=[{0: 1}])

    def test_contre_exemple_perturbation_arbitraire(self):
        # Au rang un, acceleration=(0,1). Pour toute perturbation non nulle,
        # les deux contraintes bloquent tout : écart un, indépendamment d'eps.
        m = [{0: 1}, {1: 1}]
        for eps in (F(1, 10), F(1, 10**12), F(1, 10**100)):
            g = [{0: 1}, {0: 1, 1: eps}]
            sol = dense_exact(kkt(m, g), [0, 1, 0, 0])
            self.assertEqual(sol[:2], [0, 0])
            self.assertEqual(sol[3], 1/eps)
            self.assertEqual(sol[2], -1/eps)


if __name__ == '__main__': unittest.main()
