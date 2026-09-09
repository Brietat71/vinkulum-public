"""Propositions numériques pour la certification creuse, distinctes de la preuve."""
from fractions import Fraction as F
import math

import numpy as np

from ._certification_creuse import DIMENSION_MAX, TERMES_MAX, certifier
from ._verification_creuse import ecrire_matrice, permuter, permutation


def scalaire(v):
    if isinstance(v, (bool, np.bool_)) or isinstance(v, (complex, np.complexfloating)):
        raise ValueError('scalaire réel exact ou binary64 requis')
    if type(v) in (int, F):
        return F(v)
    if isinstance(v, (float, np.floating, np.integer)):
        f = float(v)
        if math.isfinite(f):
            return F(f)
    raise ValueError('scalaire fini requis')


def vecteur(v, n):
    if np.shape(v) != (n,):
        raise ValueError('dimension de vecteur invalide')
    return [scalaire(x) for x in v]


def creuse(a, forme):
    from scipy.sparse import issparse
    if not issparse(a) or a.shape != forme or a.dtype.kind not in 'fiu':
        raise ValueError('matrice SciPy creuse réelle de dimension conforme requise')
    if a.nnz > TERMES_MAX:
        raise ValueError('budget de contributions dépassé')
    coo = a.tocoo(copy=True)
    out = [{} for _ in range(forme[0])]
    for i, j, v in zip(coo.row, coo.col, coo.data, strict=True):
        i, j = int(i), int(j)
        if not 0 <= i < forme[0] or not 0 <= j < forme[1]:
            raise ValueError('indice creux hors domaine')
        out[i][j] = out[i].get(j, F(0))+scalaire(v)
    return [{j: v for j, v in row.items() if v} for row in out]


def proposition(a, b):
    from scipy.sparse import csc_matrix
    from scipy.sparse.linalg import splu
    rows, cols, data = [], [], []
    for i, row in enumerate(a):
        for j, v in row.items():
            fv = float(v)
            if not math.isfinite(fv):
                raise ValueError('proposition flottante non représentable')
            rows.append(i); cols.append(j); data.append(fv)
    approx = csc_matrix((data, (rows, cols)), shape=(len(a), len(a)))
    lu = splu(approx)
    bx = np.array([float(v) for v in b])
    if not np.isfinite(bx).all():
        raise ValueError('second membre de proposition non représentable')
    x = lu.solve(bx)
    pr, pc = np.argsort(lu.perm_r).tolist(), np.argsort(lu.perm_c).tolist()
    return vecteur(x, len(a)), creuse(lu.L, approx.shape), creuse(lu.U, approx.shape), pr, pc


def systeme(a, b, x=None, *, facteurs=None, permutations=None,
             delta_a=None, delta_b=None, poids=None):
    from scipy.sparse import issparse
    if not issparse(a) or len(a.shape) != 2 or a.shape[0] != a.shape[1] or not 1 <= a.shape[0] <= DIMENSION_MAX:
        raise ValueError('matrice carrée creuse dans le budget [1,4096] requise')
    n = a.shape[0]
    aa, bb = creuse(a, (n, n)), vecteur(b, n)
    da = [{} for _ in range(n)] if delta_a is None else creuse(delta_a, (n, n))
    db = [F(0)]*n if delta_b is None else vecteur(delta_b, n)
    w = [F(1)]*n if poids is None else vecteur(poids, n)
    return systeme_rationnel(aa, bb, x, facteurs=facteurs, permutations=permutations, da=da, db=db, w=w)


def systeme_rationnel(aa, bb, x, *, facteurs=None, permutations=None, da, db, w):
    n = len(aa)
    if facteurs is None:
        if permutations is not None:
            raise ValueError('permutations réservées aux facteurs fournis')
        propose, l, u, pr, pc = proposition(aa, bb)
        xx = propose if x is None else vecteur(x, n)
    else:
        if x is None or not isinstance(facteurs, (tuple, list)) or len(facteurs) != 2:
            raise ValueError('deux facteurs et un candidat explicite requis')
        l, u = [creuse(v, (n, n)) for v in facteurs]
        xx = vecteur(x, n)
        if permutations is None:
            pr, pc = list(range(n)), list(range(n))
        else:
            if not isinstance(permutations, (tuple, list)) or len(permutations) != 2:
                raise ValueError('deux permutations requises')
            pr, pc = [permutation(p, n) for p in permutations]
    preuve = certifier(permuter(aa, pr, pc), [bb[i] for i in pr], [xx[i] for i in pc], l, u,
                        poids=[w[i] for i in pc], delta_a=permuter(da, pr, pc), delta_b=[db[i] for i in pr])
    bornes = [F(0)]*n
    for i, j in enumerate(pc):
        bornes[j] = preuve['bornes_composantes'][i]
    return {'schema': 'vinkulum.lineaire.creux.1', 'dimension': n, 'a': ecrire_matrice(aa),
            'b': list(map(str, bb)), 'x': list(map(str, xx)), 'l': ecrire_matrice(l), 'u': ecrire_matrice(u),
            'pr': pr, 'pc': pc, 'poids': list(map(str, w)), 'delta_a': ecrire_matrice(da),
            'delta_b': list(map(str, db)), 'bornes_composantes': list(map(str, bornes)),
            'contraction': str(preuve['contraction'])}


def quotient(m, c, t, base, force, d, x=None, *, delta_m=None, delta_c=None,
             delta_t=None, delta_force=None, delta_d=None, poids=None):
    from scipy.sparse import issparse
    from ._certification_creuse import kkt
    from ._verification_structurelle import verifier_document
    if not all(issparse(a) for a in (m, c, t)):
        raise ValueError('matrices SciPy creuses requises')
    n, r, nb = m.shape[0], c.shape[0], t.shape[0]
    if not 1 <= r <= n or n+r > DIMENSION_MAX or not r <= nb <= DIMENSION_MAX:
        raise ValueError('dimensions structurelles hors domaine')
    mm, cc, tt = creuse(m, (n, n)), creuse(c, (r, n)), creuse(t, (nb, r))
    dm = [{} for _ in range(n)] if delta_m is None else creuse(delta_m, (n, n))
    dc = [{} for _ in range(r)] if delta_c is None else creuse(delta_c, (r, n))
    dt = [{} for _ in range(nb)] if delta_t is None else creuse(delta_t, (nb, r))
    f, dd = vecteur(force, n), vecteur(d, r)
    df = [F(0)]*n if delta_force is None else vecteur(delta_force, n)
    ddd = [F(0)]*r if delta_d is None else vecteur(delta_d, r)
    w = [F(1)]*(n+r) if poids is None else vecteur(poids, n+r)
    systeme = systeme_rationnel(kkt(mm, cc), f+dd, x, da=kkt(dm, dc), db=df+ddd, w=w)
    xx, errors = [list(map(F, systeme[nom])) for nom in ('x', 'bornes_composantes')]
    from ._verification_structurelle import reaction
    rp, br = reaction(cc, dc, xx[n:], errors[n:], n)
    doc = {'schema': 'vinkulum.quotient.structurel.1', 'dimension_physique': n,
           'rang': r, 'nombre_contraintes': nb, 'c': ecrire_matrice(cc),
           't': ecrire_matrice(tt), 'base': list(base), 'delta_c': ecrire_matrice(dc),
           'delta_t': ecrire_matrice(dt), 'systeme': systeme,
           'reaction_proposee': list(map(str, rp)), 'bornes_reaction': list(map(str, br))}
    verifier_document(doc)
    return doc
