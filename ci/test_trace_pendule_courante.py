"""Contre-épreuve temporelle du noyau installé, avec limites gelées du banc."""
from fractions import Fraction as F
import unittest

from qualifie_trajectoire_certifiee import generer


class TraceCourante(unittest.TestCase):
    def test_pendule_64(self):
        doc = generer(64)
        # Limites publiées pour le banc 0.15.0, pas une estimation à partir
        # d'une comparaison entre deux exécutions du même intégrateur.
        limites = {'position': F(1084, 10**7), 'rotation': F(1084, 10**7),
                   'vitesse': F(1580, 10**7), 'omega': F(1438, 10**7)}
        for nom, limite in limites.items():
            self.assertLessEqual(F(doc['majorants'][nom]), limite, nom)


if __name__ == '__main__':
    unittest.main()
