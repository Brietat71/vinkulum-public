"""Vérification autonome d'une trace de pendule, bibliothèque standard seule.

Le contrat de ce premier format est volontairement un seul banc physique.
Les dates, paramètres, incertitudes et observables sont imposés par le format.
Il prouve des écarts aux nœuds, sans authentifier l'origine de la trace.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import re
import sys

from intervalle_temporel import BITS, I, sincos
from prototype_trajectoire_certifiee import verifier_pas


SCHEMA = 'vinkulum.trajectoire.pendule.prototype.1'
MODELE = 'pivot-y;m=1;L=1;J=identite;gravite=-g*ez;theta_ddot=-g/2*sin(theta)'
G64 = F(5522539043063071, 562949953421312)
INCERTITUDE = F(1, 10**12)
INITIAL = [I(F(1, 2)-INCERTITUDE, F(1, 2)+INCERTITUDE),
           I(-INCERTITUDE, INCERTITUDE)]
COEFFICIENT = I(F(981, 200), G64/2)
NOMS = ('position', 'rotation', 'vitesse', 'omega')
TAILLES = (3, 9, 3, 3)
MAX_OCTETS = 4_000_000


def cles(objet, attendues):
    if type(objet) is not dict or set(objet) != set(attendues):
        raise ValueError('champs du document non conformes')


def rationnel(texte):
    if type(texte) is not str or len(texte) > 1300:
        raise ValueError('rationnel textuel requis')
    if re.fullmatch(r'-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?', texte) is None:
        raise ValueError('syntaxe rationnelle non canonique')
    try:
        q = F(texte)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError('rationnel invalide') from exc
    if str(q) != texte or max(abs(q.numerator).bit_length(), q.denominator.bit_length()) > 2048:
        raise ValueError('rationnel non canonique ou trop grand')
    return q


def liste(valeurs, taille):
    if type(valeurs) is not list or len(valeurs) != taille:
        raise ValueError('dimension non conforme')
    return valeurs


def lire_intervalle(valeurs):
    lo, hi = map(rationnel, liste(valeurs, 2))
    v = I(lo, hi)
    if (v.lo, v.hi) != (lo, hi):
        raise ValueError('intervalle non dyadique à 128 bits')
    return v


def lire_boite(valeurs):
    return [lire_intervalle(v) for v in liste(valeurs, 2)]


def ecrire_intervalle(v):
    return [str(v.lo), str(v.hi)]


def ecrire_boite(x):
    return list(map(ecrire_intervalle, x))


def observables(x):
    s, c = sincos(x[0])
    z, u = I.de(0), I.de(1)
    return {'position': [-s, z, -c],
            'rotation': [c, z, s, z, u, z, -s, z, c],
            'vitesse': [-c*x[1], z, s*x[1]],
            'omega': [z, x[1], z]}


def ecarts(x, trace):
    cles(trace, NOMS)
    reference = observables(x)
    resultat = {}
    for nom, taille in zip(NOMS, TAILLES, strict=True):
        point = list(map(rationnel, liste(trace[nom], taille)))
        # La trace est celle d'observables binary64, stockées exactement.
        for q in point:
            if not est_binary64(q):
                raise ValueError('observable non représentable en binary64')
        resultat[nom] = max(max(abs(v.lo-q), abs(v.hi-q))
                            for v, q in zip(reference[nom], point, strict=True))
    return resultat


def est_binary64(q):
    n, d = abs(q.numerator), q.denominator
    if d & (d-1):
        return False
    if not n:
        return True
    zeros = (n & -n).bit_length()-1
    significande = n >> zeros
    exposant = zeros-(d.bit_length()-1)
    return (significande.bit_length() <= 53 and exposant >= -1074
            and exposant+significande.bit_length()-1 <= 1023)


def verifier_document(doc):
    cles(doc, ('schema', 'modele', 'bits', 'ordre', 'pas', 'initial',
               'coefficient', 'noeuds', 'majorants'))
    if doc['schema'] != SCHEMA or doc['modele'] != MODELE:
        raise ValueError('schéma ou modèle non admis')
    if type(doc['bits']) is not int or doc['bits'] != BITS or type(doc['ordre']) is not int or doc['ordre'] != 8:
        raise ValueError('arithmétique ou ordre non admis')
    h = rationnel(doc['pas'])
    if h not in (F(1, 64), F(1, 128), F(1, 256)):
        raise ValueError('pas hors contrat')
    x, k = lire_boite(doc['initial']), lire_intervalle(doc['coefficient'])
    if x != INITIAL or k != COEFFICIENT:
        raise ValueError('conditions initiales ou coefficient hors contrat')
    noeuds = liste(doc['noeuds'], h.denominator+1)
    maxima = dict.fromkeys(NOMS, F(0))
    largeur = F(0)
    for j, ligne in enumerate(noeuds):
        cles(ligne, ('t', 'domaine', 'boite', 'trace'))
        if rationnel(ligne['t']) != j*h:
            raise ValueError('date manquante ou grille non conforme')
        if j == 0:
            if ligne['domaine'] is not None:
                raise ValueError('pas de domaine avant le nœud initial')
        else:
            x = verifier_pas(x, lire_boite(ligne['domaine']), h, k, doc['ordre'])
        if lire_boite(ligne['boite']) != x:
            raise ValueError('boîte annoncée différente du calcul validé')
        erreurs = ecarts(x, ligne['trace'])
        for nom in NOMS:
            maxima[nom] = max(maxima[nom], erreurs[nom])
        largeur = max(largeur, *(v.hi-v.lo for v in x))
    cles(doc['majorants'], NOMS)
    for nom in NOMS:
        if rationnel(doc['majorants'][nom]) < maxima[nom]:
            raise ValueError('majorant sous-estimé : '+nom)
    return {'noeuds': len(noeuds), 'horizon': '1',
            'majorants_recalcules': {nom: str(maxima[nom]) for nom in NOMS},
            'largeur_reference_max': str(largeur)}


def sans_doublons(paires):
    d = {}
    for k, v in paires:
        if k in d:
            raise ValueError('clé JSON dupliquée')
        d[k] = v
    return d


def lire_document(chemin):
    with Path(chemin).open('rb') as flux:
        contenu = flux.read(MAX_OCTETS+1)
    if len(contenu) > MAX_OCTETS:
        raise ValueError('document trop volumineux')
    try:
        return json.loads(contenu, object_pairs_hook=sans_doublons)
    except (UnicodeDecodeError, RecursionError) as exc:
        raise ValueError('document illisible') from exc


if __name__ == '__main__':
    try:
        if len(sys.argv) != 2:
            raise ValueError('usage : archive_trajectoire_certifiee.py certificat.json')
        print(json.dumps(verifier_document(lire_document(sys.argv[1])), indent=2))
    except (ValueError, OSError) as exc:
        print('Certificat temporel refusé : '+str(exc), file=sys.stderr)
        sys.exit(1)
