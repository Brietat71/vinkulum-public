from fractions import Fraction as F
from math import factorial
import unittest
from intervalle_temporel import I, bas, haut, sincos, GRILLE


class Intervalles(unittest.TestCase):
    def test_arrondis_diriges_signes(self):
        for a in range(-20, 21):
            for b in range(1, 20):
                q = F(a,b)
                self.assertLessEqual(bas(q), q)
                self.assertGreaterEqual(haut(q), q)
                self.assertLessEqual(haut(q)-bas(q), F(1,GRILLE))

    def test_trigonometrie_contre_somme_rationnelle_longue(self):
        for x in [F(k,10) for k in range(-20,21)]:
            s,c=sincos(x)
            so=sum(((-1)**j*x**(2*j+1)/factorial(2*j+1) for j in range(35)), F(0))
            co=sum(((-1)**j*x**(2*j)/factorial(2*j) for j in range(36)), F(0))
            es=abs(x)**71/factorial(71)
            ec=abs(x)**72/factorial(72)
            self.assertLessEqual(s.lo,so-es)
            self.assertGreaterEqual(s.hi,so+es)
            self.assertLessEqual(c.lo,co-ec)
            self.assertGreaterEqual(c.hi,co+ec)

    def test_domaines_refuses(self):
        with self.assertRaises(ValueError): sincos(I(-3,1))
        with self.assertRaises(ValueError): I(1,0)
        with self.assertRaises(ValueError): I.de(1)/0


if __name__ == '__main__': unittest.main()
