"""Contre-calcul symbolique du prototype : aucun jet de différentiation."""
from fractions import Fraction as F

if __package__:
    from ._geometrie_certifiee import Intervalle, contraintes
else:
    from _geometrie_certifiee import Intervalle, contraintes


class Polynome:
    def __init__(self, dimension, termes):
        self.dimension = dimension
        self.termes = {tuple(p): F(c) for p, c in termes.items() if c}

    def constante(self, q):
        if isinstance(q, Polynome):
            if q.dimension != self.dimension:
                raise ValueError('dimensions polynomiales différentes')
            return q
        return Polynome(self.dimension, {(0,)*self.dimension: F(q)})

    def __add__(self, other):
        b = self.constante(other)
        out = dict(self.termes)
        for p, c in b.termes.items():
            out[p] = out.get(p, F(0))+c
        return Polynome(self.dimension, out)

    __radd__ = __add__

    def __neg__(self):
        return Polynome(self.dimension, {p: -c for p, c in self.termes.items()})

    def __sub__(self, other):
        return self + -self.constante(other)

    def __rsub__(self, other):
        return self.constante(other) + -self

    def __mul__(self, other):
        b = self.constante(other)
        out = {}
        for p, a in self.termes.items():
            for q, c in b.termes.items():
                exposant = tuple(x+y for x, y in zip(p, q, strict=True))
                out[exposant] = out.get(exposant, F(0))+a*c
        return Polynome(self.dimension, out)

    __rmul__ = __mul__

    def __truediv__(self, c):
        return self*(1/F(c))

    def derivee(self, i):
        out = {}
        for p, c in self.termes.items():
            if p[i]:
                q = list(p)
                q[i] -= 1
                out[tuple(q)] = c*p[i]
        return Polynome(self.dimension, out)

    def evalue(self, boite):
        out = Intervalle.de(0)
        for p, c in self.termes.items():
            terme = Intervalle.de(c)
            for x, degre in zip(boite, p, strict=True):
                for _ in range(degre):
                    terme = terme*Intervalle.de(x)
            out += terme
        return out


def polynomes(modele):
    n = 7*modele['corps']
    if not 1 <= n <= 28:
        raise ValueError('dimension hors domaine')
    variables = [Polynome(n, {tuple(int(i == j) for j in range(n)): F(1)}) for i in range(n)]
    return contraintes(modele, variables)


def verifier(modele, centre, rayons, inverse):
    centre, rayons = list(map(F, centre)), list(map(F, rayons))
    n = 7*modele['corps']
    if not 1 <= n <= 28 or len(centre) != n or len(rayons) != n or any(r <= 0 for r in rayons):
        raise ValueError('boîte hors domaine')
    if len(inverse) != n or any(len(l) != n for l in inverse):
        raise ValueError('matrice hors domaine')
    r = [list(map(F, l)) for l in inverse]
    f, bassins = polynomes(modele)
    boite = [Intervalle(x-d, x+d) for x, d in zip(centre, rayons, strict=True)]
    if any(p.evalue(boite).bas <= 0 for p in bassins):
        raise ValueError('bassin non établi par le contre-calcul polynomial')
    f0 = [p.evalue(centre) for p in f]
    j = [[p.derivee(k).evalue(boite) for k in range(n)] for p in f]
    bornes = []
    for i in range(n):
        z = -sum((f0[k]*r[i][k] for k in range(n)), Intervalle.de(0))
        q = F(0)
        for k in range(n):
            e = Intervalle.de(i == k)-sum((j[l][k]*r[i][l] for l in range(n)), Intervalle.de(0))
            q += e.magnitude()*rayons[k]
        # Recalcul par rayons et magnitudes, et non par image matricielle.
        borne = z.magnitude()+q
        if borne >= rayons[i]:
            raise ValueError('inclusion non établie par le contre-calcul polynomial')
        bornes.append(borne)
    return bornes
