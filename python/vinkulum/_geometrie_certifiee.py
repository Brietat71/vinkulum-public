"""Calculs privés du certificat géométrique : intervalles rationnels et jets.

Les paramètres stockés sont interprétés exactement. Le domaine admis et la
base de confiance sont explicités dans docs/CERTIFICATION_ASSEMBLAGE.md.
"""
from dataclasses import dataclass
from fractions import Fraction as F


@dataclass(frozen=True)
class Intervalle:
    bas: F
    haut: F

    def __post_init__(self):
        object.__setattr__(self, 'bas', F(self.bas))
        object.__setattr__(self, 'haut', F(self.haut))
        if self.bas > self.haut:
            raise ValueError('intervalle vide')

    @staticmethod
    def de(x):
        return x if isinstance(x, Intervalle) else Intervalle(F(x), F(x))

    def __add__(self, other):
        b = self.de(other)
        return Intervalle(self.bas + b.bas, self.haut + b.haut)

    __radd__ = __add__

    def __neg__(self):
        return Intervalle(-self.haut, -self.bas)

    def __sub__(self, other):
        return self + -self.de(other)

    def __rsub__(self, other):
        return self.de(other) + -self

    def __mul__(self, other):
        b = self.de(other)
        p = (self.bas*b.bas, self.bas*b.haut, self.haut*b.bas, self.haut*b.haut)
        return Intervalle(min(p), max(p))

    __rmul__ = __mul__

    def __truediv__(self, q):
        return self * (1/F(q))

    def magnitude(self):
        return max(abs(self.bas), abs(self.haut))


@dataclass(frozen=True)
class Differentielle:
    valeur: Intervalle
    derivees: tuple

    def constante(self, q):
        return q if isinstance(q, Differentielle) else Differentielle(
            Intervalle.de(q), (Intervalle.de(0),)*len(self.derivees))

    def __add__(self, other):
        b = self.constante(other)
        return Differentielle(self.valeur+b.valeur,
                              tuple(x+y for x, y in zip(self.derivees, b.derivees, strict=True)))

    __radd__ = __add__

    def __neg__(self):
        return Differentielle(-self.valeur, tuple(-v for v in self.derivees))

    def __sub__(self, other):
        return self + -self.constante(other)

    def __rsub__(self, other):
        return self.constante(other) + -self

    def __mul__(self, other):
        b = self.constante(other)
        return Differentielle(self.valeur*b.valeur, tuple(
            x*b.valeur+self.valeur*y for x, y in zip(self.derivees, b.derivees, strict=True)))

    __rmul__ = __mul__

    def __truediv__(self, q):
        return self * (1/F(q))


def transpose(a):
    return list(map(list, zip(*a, strict=True)))


def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col, strict=True)) for col in transpose(b)] for row in a]


def mv(a, v):
    return [sum(x*y for x, y in zip(row, v, strict=True)) for row in a]


IDENTITE = [[F(i == j) for j in range(3)] for i in range(3)]


def rotation(q):
    w, x, y, z = q
    return [[w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y)],
            [2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x)],
            [2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z]]


def contraintes(modele, x):
    """Toutes les lignes mécaniques, normes unitaires, puis jauges explicites."""
    if set(modele) != {'corps', 'contraintes', 'jauges'}:
        raise ValueError('paramètre de modèle non qualifié')
    n = modele['corps']
    if type(n) is not int or not 1 <= n <= 4 or len(x) != 7*n:
        raise ValueError('dimension géométrique hors domaine')
    poses = [(x[7*i:7*i+3], rotation(x[7*i+3:7*i+7])) for i in range(n)]

    def pose(i):
        if i is None:
            return [F(0)]*3, IDENTITE
        if type(i) is not int or not 0 <= i < n:
            raise ValueError('référence de corps invalide')
        return poses[i]

    out, bassins = [], []
    for liaison in modele['contraintes']:
        ra, a = pose(liaison['a'])
        rb, b = pose(liaison['b'])
        pa = mv(a, list(map(F, liaison['pa'])))
        pb = mv(b, list(map(F, liaison['pb'])))
        d = [rb[k]+pb[k]-ra[k]-pa[k] for k in range(3)]
        if liaison['type'] == 'distance':
            if set(liaison) != {'type', 'a', 'b', 'pa', 'pb', 'longueur'}:
                raise ValueError('paramètre de distance non qualifié')
            longueur = F(liaison['longueur'])
            if longueur <= 0:
                raise ValueError('distance strictement positive requise')
            out.append(sum(v*v for v in d)-longueur*longueur)
        elif liaison['type'] == 'liaison':
            base = {'type', 'a', 'b', 'pa', 'pb', 'bt', 'br'}
            if set(liaison) not in (base, base | {'ra', 'rb'}):
                raise ValueError('paramètre de liaison non qualifié')
            fa = liaison.get('ra', IDENTITE)
            fb = liaison.get('rb', IDENTITE)
            if len(fa) != 3 or len(fb) != 3 or any(len(l) != 3 for l in list(fa)+list(fb)):
                raise ValueError('repères 3 par 3 requis')
            ua = mm(a, [list(map(F, l)) for l in fa])
            ub = mm(b, [list(map(F, l)) for l in fb])
            for indices in (liaison['bt'], liaison['br']):
                if any(type(k) is not int or not 0 <= k < 3 for k in indices) or len(set(indices)) != len(indices):
                    raise ValueError('lignes de liaison invalides')
            t = mv(transpose(ua), d)
            e = mm(transpose(ua), ub)
            r = [(e[2][1]-e[1][2])/2, (e[0][2]-e[2][0])/2, (e[1][0]-e[0][1])/2]
            out.extend(t[k] for k in liaison['bt'])
            out.extend(r[k] for k in liaison['br'])
            if len(liaison['br']) == 3:
                bassins.append(sum(e[k][k] for k in range(3)))
            elif len(liaison['br']) == 2:
                libre = next(k for k in range(3) if k not in liaison['br'])
                bassins.append(e[libre][libre])
        else:
            raise ValueError('contrainte non qualifiée')
    out.extend(sum(v*v for v in x[7*i+3:7*i+7])-1 for i in range(n))
    indices = [i for i, _ in modele['jauges']]
    if len(set(indices)) != len(indices) or any(type(i) is not int or not 0 <= i < 7*n for i in indices):
        raise ValueError('jauges invalides')
    out.extend(x[i]-F(v) for i, v in modele['jauges'])
    if len(out) != 7*n:
        raise ValueError('toutes les contraintes et jauges doivent former un système carré')
    zero = x[0]*0
    return [v+zero for v in out], [v+zero for v in bassins]


def evaluation(modele, boite):
    n = len(boite)
    variables = [Differentielle(Intervalle.de(v), tuple(Intervalle.de(i == j) for j in range(n)))
                 for i, v in enumerate(boite)]
    f, bassins = contraintes(modele, variables)
    return [v.valeur for v in f], [list(v.derivees) for v in f], [v.valeur for v in bassins]


def inclusion(modele, centre, rayons, inverse):
    """Recalcule l'inclusion stricte, sans se fier à un statut enregistré."""
    centre, rayons = list(map(F, centre)), list(map(F, rayons))
    n = len(centre)
    if not 1 <= n <= 28 or len(rayons) != n or any(r <= 0 for r in rayons):
        raise ValueError('boîte finie de rayon strictement positif requise')
    if len(inverse) != n or any(len(l) != n for l in inverse):
        raise ValueError('inverse de dimension incorrecte')
    inverse = [list(map(F, row)) for row in inverse]
    boite = [Intervalle(x-r, x+r) for x, r in zip(centre, rayons, strict=True)]
    f0, _, _ = evaluation(modele, centre)
    _, j, bassins = evaluation(modele, boite)
    if any(v.bas <= 0 for v in bassins):
        raise ValueError('bassin de rotation non établi')
    z = mv(inverse, [-v for v in f0])
    rj = mm(inverse, j)
    erreur = [[F(i == k)-rj[i][k] for k in range(n)] for i in range(n)]
    image = [z[i]+sum(erreur[i][k]*Intervalle(-rayons[k], rayons[k]) for k in range(n))
             for i in range(n)]
    if any(not -r < v.bas or not v.haut < r for r, v in zip(rayons, image, strict=True)):
        raise ValueError('inclusion stricte de Krawczyk non établie')
    eta = max(sum(erreur[i][k].magnitude()*rayons[k] for k in range(n))/rayons[i]
              for i in range(n))
    if eta >= 1:
        raise ValueError('contraction pondérée non établie')
    return {'rayons': rayons, 'image': image, 'contraction': eta, 'residu': f0}


def proposer_inverse(modele, centre):
    import numpy as np
    _, j, _ = evaluation(modele, centre)
    return np.linalg.inv([[float(v.bas) for v in row] for row in j]).tolist()


def depuis_noyau(noyau, jauges):
    """Extraction privée du prototype ; l'inverse et le quaternion sont proposés."""
    from scipy.spatial.transform import Rotation
    photographie = noyau._geometrie_pour_certificat()
    modele = {k: photographie[k] for k in ('corps', 'contraintes')}
    modele['jauges'] = list(jauges)
    lignes = sum(len(c['bt'])+len(c['br']) if c['type'] == 'liaison' else 1
                 for c in modele['contraintes'])
    if lignes != photographie['lignes_originales']:
        raise ValueError('nombre de contraintes originales incohérent')
    centre = []
    for p in photographie['poses']:
        centre.extend(F(a)+F(b) for a, b in zip(p['position'], p['position_basse'], strict=True))
        x, y, z, w = Rotation.from_matrix(p['rotation']).as_quat()
        centre.extend(map(F, [float(w), float(x), float(y), float(z)]))
    return modele, centre, photographie


def bornes_pose(photographie, centre, rayons):
    """Distance composante à la pose native, conversion quaternion comprise.

    La borne de rotation est le carré de la norme de Frobenius de R-R_native,
    et non une différence de quaternions ni une distance angulaire annoncée.
    """
    centre, rayons = list(map(F, centre)), list(map(F, rayons))
    if len(centre) != 7*len(photographie['poses']) or len(rayons) != len(centre):
        raise ValueError('dimension de pose incohérente')
    out = []
    for i, p in enumerate(photographie['poses']):
        x = centre[7*i:7*i+7]
        r = rayons[7*i:7*i+7]
        if any(v <= 0 for v in r):
            raise ValueError('rayon de pose non positif')
        position = [F(a)+F(b) for a, b in zip(p['position'], p['position_basse'], strict=True)]
        rp = [abs(x[k]-position[k])+r[k] for k in range(3)]
        matrice = rotation([Intervalle(v-d, v+d) for v, d in zip(x[3:], r[3:], strict=True)])
        delta = [[(matrice[j][k]-F(p['rotation'][j][k])).magnitude() for k in range(3)]
                 for j in range(3)]
        out.append({'position': rp,
                    'rotation_frobenius_carree': sum(d*d for row in delta for d in row)})
    return out
