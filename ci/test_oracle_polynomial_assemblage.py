"""Confrontation des deux calculs rationnels du Jacobien et de l'inclusion."""
from fractions import Fraction as F
import random
import unittest

from prototype_assemblage_certifie import evaluation, inclusion, proposer_inverse, rotation
from oracle_polynomial_assemblage import polynomes, verifier
from test_prototype_assemblage_certifie import pendule, CENTRE


class OraclePolynomial(unittest.TestCase):
    def test_derivation_symbolique_contre_jets(self):
        modele = pendule()
        modele['corps'] = 2
        modele['contraintes'].append({'type': 'liaison', 'a': 0, 'b': 1,
            'pa': [0, 0, -1], 'pb': [0, 0, 1], 'bt': [0, 1, 2], 'br': [0, 2]})
        modele['jauges'].append((12, 0))
        for c in modele['contraintes']:
            c['ra'] = rotation([F(3,5), 0, F(4,5), 0])
            c['rb'] = rotation([F(3,5), F(4,5), 0, 0])
        f, bassins = polynomes(modele)
        derivees = [[p.derivee(k) for k in range(14)] for p in f]
        rng = random.Random(0xCE4704)
        for _ in range(20):
            x = [F(rng.randrange(-20, 21), 7) for _ in range(14)]
            fa, ja, ba = evaluation(modele, x)
            self.assertEqual(fa, [p.evalue(x) for p in f])
            self.assertEqual(ja, [[p.evalue(x) for p in row] for row in derivees])
            self.assertEqual(ba, [p.evalue(x) for p in bassins])

    def test_inclusion_refaite_et_inverse_altere_refuse(self):
        m = pendule()
        x = [F(v)+F(i+1, 10**9) for i, v in enumerate(CENTRE)]
        rayons = [F(1, 10**6)]*7
        inverse = proposer_inverse(m, x)
        inclusion(m, x, rayons, inverse)
        bornes = verifier(m, x, rayons, inverse)
        self.assertTrue(all(b < r for b, r in zip(bornes, rayons, strict=True)))
        inverse[0] = [0]*7
        with self.assertRaisesRegex(ValueError, 'inclusion'):
            verifier(m, x, rayons, inverse)

    def test_singularite_et_branche_retournee_refusees(self):
        import numpy as np
        m = pendule()
        rayons = [F(1, 10**6)]*7
        m['jauges'] = [(1, 0)]
        with self.assertRaisesRegex(ValueError, 'inclusion'):
            verifier(m, CENTRE, rayons, np.eye(7).tolist())
        m = pendule()
        with self.assertRaisesRegex(ValueError, 'bassin'):
            verifier(m, [0, 0, 1, 0, 1, 0, 0], rayons, np.eye(7).tolist())


if __name__ == '__main__':
    unittest.main()
