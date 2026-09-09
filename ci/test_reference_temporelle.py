"""Contre-épreuves indépendantes : mouvement libre et séparatrice exacte."""
from fractions import Fraction as F
from math import factorial
import unittest
from intervalle_temporel import I, sincos
from prototype_trajectoire_certifiee import pas


def div_positive(a,b):
    assert a.lo >= 0 and b.lo > 0
    return I(a.lo/b.hi,a.hi/b.lo)


def exponentielle(t):
    assert 0 <= t <= 1
    p=sum((t**j/factorial(j) for j in range(61)),F(0))
    # Taylor-Lagrange, exp(xi) < 3 pour xi dans [0,1].
    return I(p,p+3*t**61/factorial(61))


class Reference(unittest.TestCase):
    def test_mouvement_libre_dyadique_exact(self):
        x=[I.de(F(1,2)),I.de(F(1,4))]
        for _ in range(128): x,_=pas(x,F(1,128),0)
        self.assertEqual(x,[I.de(F(3,4)),I.de(F(1,4))])

    def test_separatrice_par_exponentielle_validee(self):
        for signe in (-1,1):
            x=[I.de(0),I.de(2*signe)]
            for j in range(1,129):
                x,_=pas(x,F(1,128),1)
                u=exponentielle(F(j,128))
                v=div_positive(4*u,1+u*u)
                if signe < 0: v=-v
                self.assertLessEqual(x[1].lo,v.lo)
                self.assertGreaterEqual(x[1].hi,v.hi)
                self.assertLess(x[1].hi-x[1].lo,F(1,10**10))
                _,c=sincos(x[0])
                energie=x[1]*x[1]/2+1-c
                self.assertLessEqual(energie.lo,2)
                self.assertGreaterEqual(energie.hi,2)


if __name__ == '__main__': unittest.main()
