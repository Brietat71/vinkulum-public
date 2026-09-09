"""Documents creux autonomes : permutations, enveloppes et bornes recalculées."""
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import re
import sys

from prototype_lineaire_creux import DIMENSION_MAX, TERMES_MAX, certifier

SCHEMA = 'vinkulum.lineaire.creux.prototype.1'


def rationnel(s):
    if type(s) is not str or len(s) > 1400 or re.fullmatch(r'-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?', s) is None:
        raise ValueError('rationnel canonique requis')
    q = F(s)
    if str(q) != s or max(abs(q.numerator).bit_length(), q.denominator.bit_length()) > 4096:
        raise ValueError('rationnel hors budget ou non canonique')
    return q


def liste(x, n):
    if type(x) is not list or len(x) != n:
        raise ValueError('dimension de liste invalide')
    return x


def lire_matrice(a, n):
    liste(a, n)
    if sum(len(l) for l in a if type(l) is list) > TERMES_MAX:
        raise ValueError('budget de termes dépassé')
    out = []
    for ligne in a:
        if type(ligne) is not list:
            raise ValueError('ligne creuse invalide')
        row, dernier = {}, -1
        for entree in ligne:
            j, v = liste(entree, 2)
            if type(j) is not int or not dernier < j < n:
                raise ValueError('indices non canoniques')
            row[j] = rationnel(v)
            if not row[j]:
                raise ValueError('zéro explicite non canonique')
            dernier = j
        out.append(row)
    return out


def ecrire_matrice(a):
    return [[[j, str(v)] for j, v in sorted(row.items()) if v] for row in a]


def permutation(p, n):
    liste(p, n)
    if any(type(i) is not int for i in p) or set(p) != set(range(n)):
        raise ValueError('permutation invalide')
    return p


def permuter(a, pr, pc):
    inverse = {j: i for i, j in enumerate(pc)}
    return [{inverse[j]: v for j, v in a[i].items()} for i in pr]


def verifier_document(doc):
    attendus = {'schema', 'dimension', 'a', 'b', 'x', 'l', 'u', 'pr', 'pc',
                'poids', 'delta_a', 'delta_b', 'bornes_composantes', 'contraction'}
    if type(doc) is not dict or set(doc) != attendus or doc['schema'] != SCHEMA:
        raise ValueError('schéma ou champs invalides')
    n = doc['dimension']
    if type(n) is not int or not 1 <= n <= DIMENSION_MAX:
        raise ValueError('dimension hors domaine')
    pr, pc = permutation(doc['pr'], n), permutation(doc['pc'], n)
    a, l, u, da = [lire_matrice(doc[nom], n) for nom in ('a', 'l', 'u', 'delta_a')]
    b, x, w, db, bornes = [list(map(rationnel, liste(doc[nom], n)))
                           for nom in ('b', 'x', 'poids', 'delta_b', 'bornes_composantes')]
    preuve = certifier(permuter(a, pr, pc), [b[i] for i in pr], [x[i] for i in pc], l, u,
                        poids=[w[i] for i in pc], delta_a=permuter(da, pr, pc), delta_b=[db[i] for i in pr])
    q = rationnel(doc['contraction'])
    if not preuve['contraction'] <= q < 1:
        raise ValueError('contraction annoncée non établie')
    originales = [F(0)]*n
    for i, original in enumerate(pc):
        originales[original] = preuve['bornes_composantes'][i]
        if bornes[original] < originales[original]:
            raise ValueError('borne composante sous-estimée')
    return {'dimension': n, 'contraction': str(preuve['contraction']),
            'bornes_composantes': list(map(str, originales)),
            'produits_facteurs': preuve['produits_facteurs'], 'termes_defaut': preuve['termes_defaut']}


def document(e, pr, pc):
    """Défait les permutations d'une proposition ; archive les données originales."""
    if set(e) != {'a', 'b', 'x', 'l', 'u'}:
        raise ValueError('constructeur nominal : arguments supplémentaires non admis')
    n = len(e['a'])
    a, b, x = [{} for _ in range(n)], [F(0)]*n, [F(0)]*n
    for i in range(n):
        a[pr[i]] = {pc[j]: v for j, v in e['a'][i].items()}
        b[pr[i]], x[pc[i]] = e['b'][i], e['x'][i]
    preuve = certifier(**e)
    bornes = [F(0)]*n
    for i, j in enumerate(pc): bornes[j] = preuve['bornes_composantes'][i]
    doc = {'schema': SCHEMA, 'dimension': n, 'a': ecrire_matrice(a), 'b': list(map(str, b)),
           'x': list(map(str, x)), 'l': ecrire_matrice(e['l']), 'u': ecrire_matrice(e['u']),
           'pr': pr, 'pc': pc, 'poids': ['1']*n, 'delta_a': [[] for _ in range(n)],
           'delta_b': ['0']*n, 'bornes_composantes': list(map(str, bornes)),
           'contraction': str(preuve['contraction'])}
    verifier_document(doc)
    return doc


def sans_doublon(paires):
    out = {}
    for k, v in paires:
        if k in out:
            raise ValueError('clé dupliquée')
        out[k] = v
    return out


def lire_document(p):
    p = Path(p)
    ouvre = gzip.open if p.suffix == '.gz' else open
    with ouvre(p, 'rb') as f:
        data = f.read(16_000_001)
    if len(data) > 16_000_000:
        raise ValueError('document trop volumineux')
    return json.loads(data, object_pairs_hook=sans_doublon)


if __name__ == '__main__':
    try:
        if len(sys.argv) != 2:
            raise ValueError('usage : archive_lineaire_creux.py certificat.json[.gz]')
        resultat = verifier_document(lire_document(sys.argv[1]))
        resultat['erreur_max'] = str(max(map(F, resultat.pop('bornes_composantes'))))
        print(json.dumps(resultat, indent=2))
    except (ValueError, OSError, EOFError, RecursionError) as exc:
        print('Certificat creux refusé : '+str(exc), file=sys.stderr)
        sys.exit(1)
