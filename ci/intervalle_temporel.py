"""Intervalles dyadiques dirigés, sans opération flottante (prototype temporel)."""
from dataclasses import dataclass
from fractions import Fraction as F
from math import factorial

BITS = 128
GRILLE = 1 << BITS


def bas(q):
    q = F(q)
    return F((q.numerator*GRILLE)//q.denominator, GRILLE)


def haut(q):
    return -bas(-F(q))


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    def __post_init__(self):
        if F(self.lo) > F(self.hi):
            raise ValueError('intervalle vide')
        object.__setattr__(self, 'lo', bas(self.lo))
        object.__setattr__(self, 'hi', haut(self.hi))

    @staticmethod
    def de(q):
        return q if isinstance(q, I) else I(q, q)

    def __add__(self, autre):
        b = I.de(autre)
        return I(self.lo+b.lo, self.hi+b.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, autre):
        return self+-I.de(autre)

    def __rsub__(self, autre):
        return I.de(autre)+-self

    def __mul__(self, autre):
        b = I.de(autre)
        p = [self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi]
        return I(min(p), max(p))

    __rmul__ = __mul__

    def __truediv__(self, q):
        q = F(q)
        if q == 0:
            raise ValueError('division par zéro')
        return I(self.lo/q, self.hi/q) if q > 0 else I(self.hi/q, self.lo/q)

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def interieur_de(self, b):
        return b.lo < self.lo and self.hi < b.hi


def sincos(x):
    x = I.de(x)
    if x.mag() > 2:
        raise ValueError('domaine trigonométrique abs(x)<=2 non établi')
    s, c = I.de(0), I.de(0)
    for k in range(40, -1, -1):
        a = F((-1)**((k-1)//2), factorial(k)) if k % 2 else F(0)
        b = F((-1)**(k//2), factorial(k)) if k % 2 == 0 else F(0)
        s, c = s*x+a, c*x+b
    reste = x.mag()**41/factorial(41)
    erreur = I(-reste, reste)
    return s+erreur, c+erreur
