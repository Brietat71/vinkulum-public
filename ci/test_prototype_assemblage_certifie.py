"""Qualification du prototype avant export natif et avant API publique."""
import copy
from fractions import Fraction as F
import itertools
import unittest

from prototype_assemblage_certifie import (
    Intervalle as I, evaluation, inclusion, proposer_inverse, rotation,
)


def pendule():
    return {'corps': 1, 'contraintes': [
        {'type': 'liaison', 'a': None, 'b': 0, 'pa': [0, 0, 0], 'pb': [0, 0, 1],
         'bt': [0, 1, 2], 'br': [0, 2]}], 'jauges': [(5, 0)]}


CENTRE = [0, 0, -1, 1, 0, 0, 0]


class Geometrie(unittest.TestCase):
    def test_produits_intervalles_oracle_sommets_et_interieurs(self):
        boites = [I(a, b) for a in range(-3, 4) for b in range(a, 4)]
        for a, b in itertools.product(boites, repeat=2):
            c = a*b
            points = [x*y for x in (a.bas, (a.bas+a.haut)/2, a.haut)
                      for y in (b.bas, (b.bas+b.haut)/2, b.haut)]
            self.assertEqual(c.bas, min(points))
            self.assertEqual(c.haut, max(points))

    def test_rotation_unitaire_rationnelle(self):
        # Quaternion rationnel unitaire, rotation non triviale.
        r = rotation([F(3, 5), 0, F(4, 5), 0])
        for i, j in itertools.product(range(3), repeat=2):
            self.assertEqual(sum(row[i]*row[j] for row in r), F(i == j))

    def test_pendule_solution_analytique_et_jacobien(self):
        m = pendule()
        f, j, bassins = evaluation(m, CENTRE)
        self.assertTrue(all(v.bas == v.haut == 0 for v in f))
        # Dérivées indépendantes des résidus (r+R(q)*ez).
        self.assertEqual(j[0][5], I(2, 2))
        self.assertEqual(j[1][4], I(-2, -2))
        self.assertEqual(j[2][3], I(2, 2))
        self.assertEqual(j[5][3], I(2, 2))
        self.assertGreater(bassins[0].bas, 0)
        resultat = inclusion(m, CENTRE, [F(1, 10**6)]*7, proposer_inverse(m, CENTRE))
        self.assertLess(resultat['contraction'], 1)

    def test_pose_et_norme_perturbees(self):
        m = pendule()
        centre = [F(v) + F(i+1, 10**9) for i, v in enumerate(CENTRE)]
        r = inclusion(m, centre, [F(1, 10**6)]*7, proposer_inverse(m, centre))
        for x, exact, rayon in zip(centre, CENTRE, r['rayons'], strict=True):
            self.assertLess(abs(x-exact), rayon)
        self.assertTrue(any(v.bas != 0 for v in r['residu']))

    def test_pendule_angle_non_nul(self):
        m = pendule()
        q = [F(3, 5), 0, F(4, 5), 0]
        rot = rotation(q)
        x = [-row[2] for row in rot] + q
        m['jauges'] = [(5, q[2])]
        r = inclusion(m, x, [F(1, 10**6)]*7, proposer_inverse(m, x))
        self.assertTrue(all(v.bas == 0 for v in r['residu']))

    def test_toutes_les_lignes_meme_incompatibles(self):
        m = pendule()
        autre = copy.deepcopy(m['contraintes'][0])
        autre['pa'] = [1, 0, 0]
        m['contraintes'].append(autre)
        with self.assertRaisesRegex(ValueError, 'système carré'):
            evaluation(m, CENTRE)

    def test_dependance_dans_un_systeme_carre(self):
        m = pendule()
        # La jauge y reproduit ici une direction déjà fixée au premier ordre.
        m['jauges'] = [(1, 0)]
        import numpy as np
        with self.assertRaises(np.linalg.LinAlgError):
            proposer_inverse(m, CENTRE)
        with self.assertRaisesRegex(ValueError, 'inclusion'):
            inclusion(m, CENTRE, [F(1, 100)]*7, np.eye(7).tolist())

    def test_refus_branche_retournee(self):
        m = pendule()
        x = [0, 0, 1, 0, 1, 0, 0]
        import numpy as np
        with self.assertRaisesRegex(ValueError, 'bassin'):
            inclusion(m, x, [F(1, 100)]*7, np.eye(7).tolist())

    def test_refus_boite_trop_petite_et_inverse_altere(self):
        m = pendule()
        x = list(CENTRE)
        x[0] = F(1, 1000)
        with self.assertRaisesRegex(ValueError, 'inclusion'):
            inclusion(m, x, [F(1, 10**6)]*7, proposer_inverse(m, x))
        with self.assertRaisesRegex(ValueError, 'inclusion'):
            inclusion(m, CENTRE, [F(1, 100)]*7, [[0]*7 for _ in range(7)])

    def test_jauge_absente_et_parametre_inconnu_refuses(self):
        m = pendule()
        m['jauges'] = []
        with self.assertRaisesRegex(ValueError, 'système carré'):
            evaluation(m, CENTRE)
        m = pendule()
        m['contraintes'][0]['cible_r'] = 0
        with self.assertRaisesRegex(ValueError, 'non qualifié'):
            evaluation(m, CENTRE)


if __name__ == '__main__':
    unittest.main()
