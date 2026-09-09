"""Vérificateur autonome des certificats linéaires (bibliothèque standard seule).

Exécution sans importer Vinkulum : python _verification_lineaire.py preuve.json
La preuve concerne les nombres encodés dans le document, pas leur provenance.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys


SCHEMA = 'vinkulum.lineaire.1'
DIMENSION_MAX = 128


def _rationnel(s):
    if not isinstance(s, str):
        raise ValueError('coefficient hexadécimal requis')
    x = float.fromhex(s)
    if not math.isfinite(x) or x.hex() != s:
        raise ValueError('binary64 fini en encodage hexadécimal canonique requis')
    return F(x)


def verifier_certificat(document):
    """Vérifie exactement unicité et bornes ; rend les bornes en rationnels.

    Refuse un certificat invalide par ValueError. Aucun solveur flottant,
    aucun appel à Vinkulum, aucune confiance dans une décision enregistrée.
    """
    if isinstance(document, dict) and document.get('schema') == 'vinkulum.quotient.structurel.1':
        if __package__:
            from ._verification_structurelle import verifier_document
        else:
            from _verification_structurelle import verifier_document
        return verifier_document(document)
    if isinstance(document, dict) and document.get('schema') in ('vinkulum.lineaire.creux.1', 'vinkulum.lineaire.creux.prototype.1'):
        if __package__:
            from ._verification_creuse import verifier_document
        else:
            from _verification_creuse import verifier_document
        resultat = verifier_document(document)
        resultat['bornes_composantes'] = list(map(F, resultat['bornes_composantes']))
        resultat['borne_erreur_inf'] = max(resultat['bornes_composantes'])
        resultat['unicite_uniforme'] = True
        return resultat
    if isinstance(document, dict) and document.get('schema') == 'vinkulum.assemblage.1':
        if __package__:
            from ._verification_assemblage import verifier_public
        else:
            from _verification_assemblage import verifier_public
        return verifier_public(document)
    if isinstance(document, dict) and document.get('schema') == 'vinkulum.quotient.1':
        return _verifier_quotient(document)
    if not isinstance(document, dict) or document.get('schema') != SCHEMA:
        raise ValueError('schéma de certificat inconnu')
    try:
        a0, r0 = document['matrice'], document['inverse_approche']
        n = len(a0)
        if not 1 <= n <= DIMENSION_MAX:
            raise ValueError('dimension hors budget')
        if len(r0) != n or any(len(l) != n for l in a0+r0):
            raise ValueError('matrices carrées de même dimension requises')
        a = [[_rationnel(v) for v in l] for l in a0]
        r = [[_rationnel(v) for v in l] for l in r0]
        b, x, bornes = ([ _rationnel(v) for v in document[k] ]
                       for k in ('second_membre', 'solution', 'bornes_composantes'))
        borne = _rationnel(document['borne_erreur_inf'])
        if any(len(v) != n for v in (b, x, bornes)):
            raise ValueError('dimension des vecteurs incohérente')
        if borne < 0 or any(v < 0 for v in bornes):
            raise ValueError('borne négative')
        # Calcul indépendant du générateur : toutes les opérations sont ici
        # rationnelles ; le générateur emploie le produit scalaire entier Rust.
        residu = [b[i]-sum((a[i][j]*x[j] for j in range(n)), F(0)) for i in range(n)]
        z = [sum((r[i][j]*residu[j] for j in range(n)), F(0)) for i in range(n)]
        sommes = []
        for i in range(n):
            sommes.append(sum((abs(F(i == j)-sum(
                (r[i][k]*a[k][j] for k in range(n)), F(0)))
                for j in range(n)), F(0)))
        eta = max(sommes)
        if eta >= 1:
            raise ValueError('contraction non établie : ||I-RA||_inf >= 1')
        # Inclusion de la boule d'erreur ; pas de division recopiée du générateur.
        if any(abs(z[i])+sommes[i]*borne > borne for i in range(n)):
            raise ValueError('borne globale insuffisante')
        if any(abs(z[i])+sommes[i]*borne > bornes[i] for i in range(n)):
            raise ValueError('borne de composante insuffisante')
        # Pour l'export natif, contrôler aussi l'agrégation des contributions.
        if 'contributions' in document:
            reconstruit = [[F(0) for _ in range(n)] for _ in range(n)]
            for i, j, v in document['contributions']:
                if type(i) is not int or type(j) is not int or not (0 <= i < n and 0 <= j < n):
                    raise ValueError('indice de contribution invalide')
                reconstruit[i][j] += _rationnel(v)
            if reconstruit != a:
                raise ValueError('matrice différente de la somme exacte des contributions')
        return {'dimension': n, 'eta_exact': eta, 'borne_erreur_inf': borne,
                'bornes_composantes': bornes}
    except (KeyError, TypeError, OverflowError, IndexError) as e:
        raise ValueError('document de certificat mal formé') from e


def _fraction(s):
    if not isinstance(s, str):
        raise ValueError('rationnel canonique requis')
    v = F(s)
    if str(v) != s:
        raise ValueError('encodage rationnel non canonique')
    return v


def _verifier_quotient(d):
    try:
        # Un seul niveau, pour ne pas admettre de preuves récursives.
        s = d['systeme']
        if not isinstance(s, dict) or s.get('schema') != SCHEMA:
            raise ValueError('certificat linéaire direct requis')
        q = verifier_certificat(s)
        total, n = q['dimension'], d['composantes_dynamiques']
        if type(n) is not int or not 1 <= n <= total:
            raise ValueError('dimension dynamique invalide')
        m = total-n
        actif = d['actives']
        if len(actif) != m or any(type(v) is not bool for v in actif):
            raise ValueError('masque booléen de contraintes requis')
        base = [i for i, oui in enumerate(actif) if oui]
        r = len(base)
        g = [[_rationnel(v) for v in l] for l in d['contraintes']]
        c = [_rationnel(v) for v in d['second_membre_contraintes']]
        t = [[_fraction(v) for v in l] for l in d['dependances']]
        if len(g) != m or any(len(l) != n for l in g) or len(c) != m:
            raise ValueError('dimensions des contraintes invalides')
        if len(t) != m or any(len(l) != r for l in t):
            raise ValueError('dimensions de la dépendance invalides')
        # Identifier le KKT avec jauge réellement certifié. Son inversibilité
        # implique le plein rang de C, sans seuil numérique de rang.
        a = [[_rationnel(v) for v in l] for l in s['matrice']]
        b = [_rationnel(v) for v in s['second_membre']]
        for i in range(m):
            if b[n+i] != (c[i] if actif[i] else 0):
                raise ValueError('second membre différent du système avec jauge')
            for j in range(n):
                attendu = g[i][j] if actif[i] else 0
                if a[n+i][j] != attendu or a[j][n+i] != attendu:
                    raise ValueError('couplage différent des contraintes actives')
            for j in range(m):
                if a[n+i][n+j] != int(i == j and not actif[i]):
                    raise ValueError('bloc des multiplicateurs différent de la jauge')
        # Vérification multiplicative indépendante de l'élimination du générateur.
        for i in range(m):
            for j in range(n):
                if sum((t[i][k]*g[base[k]][j] for k in range(r)), F(0)) != g[i][j]:
                    raise ValueError('dépendance G=T C non établie exactement')
            if sum((t[i][k]*c[base[k]] for k in range(r)), F(0)) != c[i]:
                raise ValueError('compatibilité c=T c_base non établie exactement')
        for k, i in enumerate(base):
            if t[i] != [F(k == j) for j in range(r)]:
                raise ValueError('la dépendance ne conserve pas la base')
        reconstruit = [[F(0) for _ in range(n)] for _ in range(m)]
        for i, j, v in d['contributions_contraintes']:
            if type(i) is not int or type(j) is not int or not (0 <= i < m and 0 <= j < n):
                raise ValueError('indice de contribution de contrainte invalide')
            reconstruit[i][j] += _rationnel(v)
        if reconstruit != g:
            raise ValueError('contraintes différentes de la somme de leurs contributions')
        h = q['bornes_composantes']
        hr = [_rationnel(v) for v in d['bornes_reaction_generalisee']]
        if len(hr) != n or any(v < 0 for v in hr):
            raise ValueError('bornes de réaction généralisée invalides')
        if any(sum((abs(g[j][i])*h[n+j] for j in range(m)), F(0)) > hr[i]
               for i in range(n)):
            raise ValueError('borne de réaction généralisée insuffisante')
        q.update(rang_contraintes=r, nombre_contraintes=m,
                 acceleration_unique=True, reaction_generalisee_unique=True,
                 multiplicateurs_uniques=(r == m), bornes_reaction_generalisee=hr)
        return q
    except (KeyError, TypeError, OverflowError, IndexError, ZeroDivisionError) as e:
        raise ValueError('document de quotient mal formé') from e


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage : python _verification_lineaire.py certificat.json')
    try:
        resultat = verifier_certificat(json.loads(Path(sys.argv[1]).read_text()))
    except (ValueError, OSError) as erreur:
        print('REFUS :', erreur, file=sys.stderr)
        raise SystemExit(1)
    if 'existence_unique_locale' in resultat:
        print(f"Certificat géométrique vérifié : dimension {resultat['dimension']}, existence et unicité locales.")
    elif 'rang_contraintes' in resultat:
        print(f"Certificat de quotient vérifié : rang {resultat['rang_contraintes']}/"
              f"{resultat['nombre_contraintes']}, borne du représentant "
              f"{resultat['borne_erreur_inf']}.")
    else:
        print(f"Certificat linéaire vérifié : dimension {resultat['dimension']}, "
              f"borne exacte {resultat['borne_erreur_inf']}.")
