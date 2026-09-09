"""Quotient perturbé à rang préservé par une paramétrisation explicite."""
from fractions import Fraction as F

if __package__:
    from . import _verification_creuse as lecture
    from ._certification_creuse import DIMENSION_MAX, haut, matrice
else:
    import _verification_creuse as lecture
    from _certification_creuse import DIMENSION_MAX, haut, matrice


def reaction(c, dc, nu, erreurs, n):
    rp, br = [F(0)]*n, [F(0)]*n
    for j in range(len(c)):
        for i in c[j].keys() | dc[j].keys():
            v, dv = c[j].get(i, F(0)), dc[j].get(i, F(0))
            rp[i] += v*nu[j]
            br[i] += abs(v)*erreurs[j]+dv*(abs(nu[j])+erreurs[j])
    return rp, list(map(haut, br))


def verifier_document(doc):
    cles = {'schema', 'dimension_physique', 'rang', 'nombre_contraintes', 'c', 't',
            'base', 'delta_c', 'delta_t', 'systeme', 'reaction_proposee', 'bornes_reaction'}
    if type(doc) is not dict or set(doc) != cles or doc['schema'] != 'vinkulum.quotient.structurel.1':
        raise ValueError('schéma structurel invalide')
    n, r, nb = [doc[k] for k in ('dimension_physique', 'rang', 'nombre_contraintes')]
    if any(type(i) is not int for i in (n, r, nb)) or not 1 <= r <= n or n+r > DIMENSION_MAX or not r <= nb <= DIMENSION_MAX:
        raise ValueError('dimensions structurelles hors domaine')
    c, dc = [lecture.lire_matrice(doc[k], r, n) for k in ('c', 'delta_c')]
    t, dt = [lecture.lire_matrice(doc[k], nb, r) for k in ('t', 'delta_t')]
    matrice(dc, r, colonnes=n, positif=True)
    matrice(dt, nb, colonnes=r, positif=True)
    base = lecture.liste(doc['base'], r)
    if any(type(i) is not int or not 0 <= i < nb for i in base) or len(set(base)) != r:
        raise ValueError('base structurelle invalide')
    for j, i in enumerate(base):
        if t[i] != {j: F(1)} or dt[i]:
            raise ValueError('identité des lignes de base non établie')
    s = doc['systeme']
    if type(s) is not dict or s.get('schema') != 'vinkulum.lineaire.creux.1' or s.get('dimension') != n+r:
        raise ValueError('système réduit non conforme')
    preuve = lecture.verifier_document(s)
    a, da = [lecture.lire_matrice(s[k], n+r) for k in ('a', 'delta_a')]
    # Vérifier les deux couplages, le bloc nul et leurs enveloppes. La
    # boîte linéaire peut être plus large que les perturbations corrélées
    # de C ; elle doit ici reproduire exactement leur enveloppe déclarée.
    for j in range(r):
        if a[n+j] != c[j] or da[n+j] != dc[j]:
            raise ValueError('contraintes du KKT différentes de C')
    ct, dct = [{} for _ in range(n)], [{} for _ in range(n)]
    for j in range(r):
        for i, v in c[j].items(): ct[i][j] = v
        for i, v in dc[j].items(): dct[i][j] = v
    for i in range(n):
        if {j-n: v for j, v in a[i].items() if j >= n} != ct[i]:
            raise ValueError('transposée de C non conservée')
        if {j-n: v for j, v in da[i].items() if j >= n} != dct[i]:
            raise ValueError('enveloppe transposée non conservée')
    xx = list(map(lecture.rationnel, lecture.liste(s['x'], n+r)))
    erreurs = list(map(F, preuve['bornes_composantes']))
    rp, br = reaction(c, dc, xx[n:], erreurs[n:], n)
    annonce = list(map(lecture.rationnel, lecture.liste(doc['reaction_proposee'], n)))
    bornes = list(map(lecture.rationnel, lecture.liste(doc['bornes_reaction'], n)))
    if annonce != rp or any(bornes[i] < br[i] for i in range(n)):
        raise ValueError('réaction ou majorants non établis')
    return {'dimension': n+r, 'rang_contraintes': r, 'nombre_contraintes': nb,
            'acceleration_unique': True, 'reaction_generalisee_unique': True,
            'multiplicateurs_uniques': r == nb, 'bornes_composantes': erreurs,
            'borne_erreur_inf': max(erreurs), 'bornes_reaction_generalisee': br,
            'reaction_proposee': rp, 'contrat': 'G_prime=T_prime*C_prime;T_prime[base,:]=I'}
