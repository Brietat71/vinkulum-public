from fractions import Fraction as F
import unittest
from intervalle_temporel import I
from prototype_trajectoire_certifiee import coefficients, pas, verifier_pas


class Taylor(unittest.TestCase):
    def test_coefficients_analytiques_non_lineaires(self):
        a,b=coefficients([I.de(0),I.de(1)],2,6)
        for valeur,exact in [(a[1],F(1)),(a[3],F(-1,3)),(a[5],F(1,20))]:
            self.assertLessEqual(valeur.lo,exact)
            self.assertGreaterEqual(valeur.hi,exact)

    def test_equilibre_inclus_et_pas_non_admis(self):
        r,d=pas([I.de(0),I.de(0)],F(1,128),F(981,200))
        for x in r:
            self.assertLessEqual(x.lo,0)
            self.assertGreaterEqual(x.hi,0)
            self.assertLess(x.mag(),F(1,10**14))
        with self.assertRaises(ValueError): pas([0,0],1,5)
        with self.assertRaises(ValueError): verifier_pas([I.de(0),I.de(1)],[I.de(0),I.de(1)],F(1,128),5)

    def test_reste_necessaire_contre_solution_fermee(self):
        from test_reference_temporelle import exponentielle, div_positive
        x = [I.de(0), I.de(2)]
        h = F(1, 2)
        domaine = [I(F(-1, 100), F(11, 10)), I(1, F(21, 10))]
        y = verifier_pas(x, domaine, h, 1, ordre=2)
        u = exponentielle(h)
        exact = div_positive(4*u, 1+u*u)
        self.assertLessEqual(y[1].lo, exact.lo)
        self.assertGreaterEqual(y[1].hi, exact.hi)
        # Sans le terme de reste le Taylor d'ordre 1 prédit omega=2,
        # incompatible avec la solution analytique 2/cosh(1/2).
        self.assertLess(exact.hi, 2)

    def test_dimensions_et_ordre_refuses(self):
        for x in ([0], [0, 0, 0]):
            with self.assertRaises(ValueError): pas(x, F(1, 128), 1)
        for ordre in (True, 1, 13, 8.0):
            with self.assertRaises(ValueError): pas([0, 0], F(1, 128), 1, ordre)


if __name__ == '__main__': unittest.main()
