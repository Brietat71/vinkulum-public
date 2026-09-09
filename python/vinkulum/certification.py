"""Certificats a posteriori : unicité et erreur en avant des systèmes linéaires.

La preuve concerne les coefficients binary64, après leur calcul. Elle ne
certifie ni leur origine mécanique, ni une trajectoire. Un refus de preuve
ne démontre pas que le système est singulier. Certification dense facultative.
"""
from fractions import Fraction as F
import math

import numpy as np

from vinkulum import _vinkulum
from vinkulum._verification_lineaire import SCHEMA, DIMENSION_MAX
from vinkulum._verification_lineaire import verifier_certificat as _verifier


__all__ = ['CertificationImpossible', 'certifier_systeme', 'certifier_initialisation',
           'verifier_certificat', 'certifier_assemblage', 'certifier_systeme_creux',
           'certifier_quotient_structurel']


def certifier_quotient_structurel(m, c, t, base, force, d, x=None, *, delta_m=None,
                                 delta_c=None, delta_t=None, delta_force=None,
                                 delta_d=None, poids=None):
    """Certifie un KKT perturbé et toutes les contraintes G'=T'C', c_complet=T'd'.

    M, C et T sont des matrices SciPy creuses, respectivement n×n, r×n et
    m×r, avec 1<=r<=n, n+r<=4096 et m<=4096. T[base,:] doit être exactement
    l'identité, avec incertitude nulle sur ces lignes. Cette structure est
    un contrat fourni ; elle n'est pas déduite de contraintes proches.

    Les enveloppes delta_* sont positives composante par composante. La
    preuve couvre toutes leurs réalisations conservant la paramétrisation.
    Le KKT réduit [M' C'^T; C' 0] est certifié uniformément inversible :
    accélération et réaction généralisée uniques. x, si fourni, contient
    l'accélération puis les multiplicateurs réduits. Sinon SuperLU propose
    un candidat. Le représentant des multiplicateurs complets porte les
    valeurs réduites sur base, et zéro ailleurs ; aucune norme minimale.

    La positivité physique de M n'est pas certifiée par cette fonction.
    SciPy requis pour produire le document, pas pour le vérifier. Les
    données sont copiées ; un refus de preuve lève CertificationImpossible.
    """
    from ._api_certification_creuse import quotient
    try:
        return quotient(m, c, t, base, force, d, x, delta_m=delta_m, delta_c=delta_c,
                         delta_t=delta_t, delta_force=delta_force, delta_d=delta_d, poids=poids)
    except (ValueError, TypeError, OverflowError, RuntimeError) as exc:
        raise CertificationImpossible(str(exc)) from exc


def certifier_systeme_creux(a, b, x=None, *, facteurs=None, permutations=None,
                           delta_a=None, delta_b=None, poids=None):
    """Certifie des systèmes creux perturbés sans inverse dense, jusqu'à 4096 inconnues.

    A est une matrice SciPy creuse carrée. Les contributions sont converties
    en binary64 puis sommées exactement, y compris les doublons. delta_a et
    delta_b sont des enveloppes positives composante par composante ; les
    coefficients absents de delta_a ont une incertitude nulle. Les poids
    strictement positifs définissent max_i |erreur_i|/poids_i.

    SuperLU propose x et des facteurs lorsque facteurs=None. Des facteurs
    (L,U) fournis exigent un candidat x explicite ; permutations=(pr,pc)
    indique A_permute[i,j]=A[pr[i],pc[j]], sinon l'ordre est inchangé.
    Les propositions flottantes ne constituent pas une preuve. Le contrôle
    rationnel établit une contraction stricte uniforme et les bornes de
    chaque composante dans l'ordre original. Les comparaisons triangulaires
    peuvent être pessimistes : un refus ne démontre pas une singularité.

    SciPy requis pour la génération, bibliothèque standard seule pour la
    relecture. CertificationImpossible en cas de domaine ou budget refusé.
    Les données sont copiées. Cette garantie concerne les coefficients
    et enveloppes déclarés, pas leur provenance mécanique ni une trajectoire.
    """
    from ._api_certification_creuse import systeme
    try:
        return systeme(a, b, x, facteurs=facteurs, permutations=permutations,
                       delta_a=delta_a, delta_b=delta_b, poids=poids)
    except (ValueError, TypeError, OverflowError, RuntimeError) as exc:
        raise CertificationImpossible(str(exc)) from exc


class CertificationImpossible(ValueError):
    """Les hypothèses ou la précision demandée ne sont pas établies."""


def _budget(dimension_max):
    if type(dimension_max) is not int or not 1 <= dimension_max <= DIMENSION_MAX:
        raise CertificationImpossible('dimension_max doit être un entier dans [1, 128]')


def _tableau(valeur, forme, nom):
    if np.iscomplexobj(valeur):
        raise CertificationImpossible(f'{nom} : coefficients réels requis')
    try:
        a = np.array(valeur, dtype=float, copy=True)
    except (ValueError, TypeError, OverflowError) as e:
        raise CertificationImpossible(f'{nom} : coefficients binary64 requis') from e
    if a.shape != forme or not np.isfinite(a).all():
        raise CertificationImpossible(f'{nom} : dimensions invalides ou donnée non finie')
    return a


def _dot(g, u, dt=0.):
    # Même produit exact que le garde de vitesse ; aucune décision de tolérance
    # de ce garde n'est utilisée ici. Seul le résidu entier est lu.
    _, r, _, _ = _vinkulum._certificat_ligne_vitesse(list(g), list(u), float(dt), 0.)
    return F(int(r, 16), 1 << 2148)


def _haut(q):
    try:
        v = float(q)
    except OverflowError as e:
        raise CertificationImpossible('borne non représentable en binary64 fini') from e
    if not math.isfinite(v):
        raise CertificationImpossible('borne non représentable en binary64 fini')
    if F(v) < q:
        v = math.nextafter(v, math.inf)
    # Vérifier le sens de l'arrondi par comparaison rationnelle, y compris
    # lors d'un sous-débordement ; ne pas faire confiance à une marge fixe.
    if not math.isfinite(v) or F(v) < q:
        raise CertificationImpossible('arrondi majorant non établi')
    return v


def certifier_systeme(a, b, x=None, *, inverse_approche=None,
                      erreur_max=None, dimension_max=64):
    """Certifie l'unicité de Ax=b et borne chaque composante de l'erreur de x.

    A carrée réelle, 1 à dimension_max inconnues (128 au plus). Les données
    sont converties en binary64 puis considérées exactes. x et l'inverse
    approchée peuvent être fournis ; sinon NumPy les propose. Aucun résultat
    de NumPy n'est tenu pour une preuve. La vérification exacte requiert
    ||I-RA||_inf < 1. Coût dense cubique, uniquement à la demande.

    erreur_max, si fourni, impose une borne absolue globale finie >= 0 sur
    les valeurs des coordonnées dans leurs unités fournies. Les coordonnées
    peuvent avoir des unités différentes : ce n'est pas une norme physique
    invariante au changement d'unités. La sortie JSON fournit des bornes
    par composante, et les données nécessaires à une vérification autonome.
    """
    _budget(dimension_max)
    forme = np.shape(a)
    if len(forme) != 2 or forme[0] != forme[1] or not 1 <= forme[0] <= dimension_max:
        raise CertificationImpossible('matrice carrée non vide dans le budget requise')
    n = forme[0]
    a, b = _tableau(a, (n, n), 'A'), _tableau(b, (n,), 'b')
    if erreur_max is not None:
        erreur_max = float(erreur_max)
        if not math.isfinite(erreur_max) or erreur_max < 0:
            raise CertificationImpossible('erreur_max doit être finie et positive ou nulle')
    try:
        x = _tableau(np.linalg.solve(a, b) if x is None else x, (n,), 'x')
        r = _tableau(np.linalg.inv(a) if inverse_approche is None else inverse_approche,
                     (n, n), 'inverse approchée')
    except np.linalg.LinAlgError as e:
        raise CertificationImpossible('construction de la solution ou de l’inverse approchée impossible') from e
    residu = [_dot(a[i], x, -b[i]) for i in range(n)]  # Ax-b
    sommes = [sum((abs(_dot(r[i], a[:, j], -float(i == j))) for j in range(n)), F(0))
              for i in range(n)]
    eta = max(sommes)
    if eta >= 1:
        raise CertificationImpossible('contraction non établie : ||I-RA||_inf >= 1')
    z = [sum((F(float(r[i, j]))*residu[j] for j in range(n)), F(0)) for i in range(n)]
    borne = _haut(max(map(abs, z))/(1-eta))
    if erreur_max is not None and F(borne) > F(erreur_max):
        raise CertificationImpossible(f'borne {borne:.6g} supérieure à erreur_max={erreur_max:.6g}')
    composantes = [_haut(abs(z[i])+sommes[i]*F(borne)) for i in range(n)]
    return {
        'schema': SCHEMA,
        'matrice': [[float(v).hex() for v in ligne] for ligne in a],
        'second_membre': [float(v).hex() for v in b],
        'solution': [float(v).hex() for v in x],
        'inverse_approche': [[float(v).hex() for v in ligne] for ligne in r],
        'borne_erreur_inf': borne.hex(),
        'bornes_composantes': [v.hex() for v in composantes],
    }


def _agrege(trip, lignes, colonnes):
    exact = [[F(0) for _ in range(colonnes)] for _ in range(lignes)]
    for i, j, v in trip:
        if not math.isfinite(v):
            raise CertificationImpossible('contribution non finie')
        exact[i][j] += F(v)
    try:
        a = [[float(v) for v in ligne] for ligne in exact]
    except OverflowError as e:
        raise CertificationImpossible('agrégation des contributions non représentable') from e
    if any(not math.isfinite(a[i][j]) or F(a[i][j]) != exact[i][j]
           for i in range(lignes) for j in range(colonnes)):
        raise CertificationImpossible('agrégation des contributions non représentable exactement')
    return a


def _dependances(g, c, actives):
    """Propose T avec G=T C par élimination exacte ; le vérificateur multiplie."""
    base = [i for i, oui in enumerate(actives) if oui]
    r, n = len(base), len(g[0]) if g else 0
    original = [[F(v) for v in g[i]] for i in base]
    ech = [l[:] for l in original]
    u = [[F(i == j) for j in range(r)] for i in range(r)]
    pivots = []
    for j in range(n):
        k = len(pivots)
        p = next((i for i in range(k, r) if ech[i][j]), None)
        if p is None:
            continue
        ech[k], ech[p] = ech[p], ech[k]
        u[k], u[p] = u[p], u[k]
        d = ech[k][j]
        ech[k], u[k] = [v/d for v in ech[k]], [v/d for v in u[k]]
        for i in range(r):
            if i != k:
                d = ech[i][j]
                ech[i] = [a-d*b for a, b in zip(ech[i], ech[k], strict=True)]
                u[i] = [a-d*b for a, b in zip(u[i], u[k], strict=True)]
        pivots.append(j)
    if len(pivots) != r:
        raise CertificationImpossible('les lignes actives sont exactement dépendantes')
    t = [[sum((F(l[pivots[k]])*u[k][j] for k in range(r)), F(0))
          for j in range(r)] for l in g]
    for i, l in enumerate(g):
        if any(sum((t[i][k]*original[k][j] for k in range(r)), F(0)) != F(l[j])
               for j in range(n)):
            raise CertificationImpossible('contrainte écartée indépendante en arithmétique exacte')
        if sum((t[i][k]*F(c[base[k]]) for k in range(r)), F(0)) != F(c[i]):
            raise CertificationImpossible('seconds membres des contraintes exactement incompatibles')
    return t


def _certifie_quotient(systeme, g, c, actives, contributions):
    t = _dependances(g, c, actives)
    n = len(systeme['solution'])-len(g)
    h = list(map(lambda v: F(float.fromhex(v)), systeme['bornes_composantes']))
    reactions = [_haut(sum((abs(F(g[j][i]))*h[n+j] for j in range(len(g))), F(0)))
                 for i in range(n)]
    return {'schema': 'vinkulum.quotient.1', 'systeme': systeme,
            'composantes_dynamiques': n, 'actives': actives,
            'contraintes': [[v.hex() for v in l] for l in g],
            'second_membre_contraintes': [v.hex() for v in c],
            'dependances': [[str(v) for v in l] for l in t],
            'contributions_contraintes': [[i, j, v.hex()] for i, j, v in contributions],
            'bornes_reaction_generalisee': [v.hex() for v in reactions]}


def certifier_initialisation(noyau, *, t=None, erreur_max=None, dimension_max=64,
                            redondances=False):
    """Borne l'erreur du système linéaire d'accélération réellement résolu par le noyau.

    Lecture seule : export natif sur un clone, avant toute intégration.
    Refuse par défaut les contraintes désactivées/redondantes, les partitions gelées,
    les contacts non lisses et un schéma de redémarrage imposé. Ne certifie
    ni la géométrie, ni les lois de forces, ni la trajectoire qui suivra.
    L'agrégation des contributions doit être représentable exactement en
    binary64 ; une perte à cette conversion est refusée, jamais masquée.

    Les 6 premières composantes par corps sont ax,ay,az en m/s² puis
    alphax,alphay,alphaz en rad/s², dans le repère spatial ; suivent les
    multiplicateurs, dans l'ordre des lignes de contraintes du noyau.

    redondances=True exige une preuve exacte G=T C et c=T c_base pour
    les lignes actives C. Le certificat établit l'accélération unique et
    la réaction généralisée unique ; les multiplicateurs restent définis
    modulo ker(G.T). Leur borne porte sur le représentant dont les lignes
    inactives sont nulles, sans revendiquer une norme minimale. Le budget
    erreur_max porte sur ce représentant et les accélérations. La sortie
    utilise le schéma vinkulum.quotient.1, vérifiable par le même outil.
    """
    _budget(dimension_max)
    if type(redondances) is not bool:
        raise CertificationImpossible('redondances doit être un booléen')
    date = noyau.t() if t is None else float(t)
    try:
        trip, b, x, physiques, actives, g_trip, c = noyau._initialisation_pour_certificat(
            date, dimension_max, redondances)
    except ValueError as e:
        raise CertificationImpossible(str(e)) from e
    n = len(b)
    a = _agrege(trip, n, n)
    certificat = certifier_systeme(a, b, x, erreur_max=erreur_max, dimension_max=dimension_max)
    certificat['contributions'] = [[i, j, v.hex()] for i, j, v in trip]
    certificat['origine'] = {'calcul': 'initialisation_lineaire', 'temps': date.hex(),
                            'composantes_dynamiques': physiques,
                            'portee': 'Système de contributions exporté ; pas une certification de trajectoire.'}
    if redondances:
        g = _agrege(g_trip, n-physiques, physiques)
        c = _tableau(c, (n-physiques,), 'second membre complet').tolist()
        certificat = _certifie_quotient(certificat, g, c, actives, g_trip)
    return certificat


def certifier_assemblage(noyau, *, jauges, rayons):
    """Certifie une pose admissible unique dans une boîte explicitement donnée.

    Lecture seule. Domaine : 1 à 4 corps rigides, liaisons holonomes sans loi imposée
    à repères matériels stockés finis, distances positives ; pas de contact,
    poutre, superélément ni corps gelé. Toutes les contraintes originales sont
    conservées. Les cas hors domaine et les inclusions non prouvées sont refusés.

    Coordonnées par corps : (rx, ry, rz, qw, qx, qy, qz). jauges contient les
    couples (indice global, valeur exacte fixée) nécessaires pour rendre le
    système carré après ajout des normes de quaternion. rayons contient 7*N
    rayons strictement positifs, dans les unités des translations et des
    quaternions. Les valeurs int, Fraction et binary64 finies sont admises.
    Aucun assemblage n'est effectué ; la boîte doit entourer la pose recherchée.

    Le document rationnel contient la pose initiale et les bornes de distance :
    composantes de translation, carré de la norme de Frobenius des rotations.
    La conversion initiale en quaternion est comprise dans ces bornes. Il ne
    certifie ni l'origine des paramètres géométriques stockés ni une trajectoire.
    """
    from ._geometrie_certifiee import depuis_noyau, proposer_inverse
    from ._verification_assemblage import fabriquer_public
    try:
        jauges = list(jauges)
        rayons = list(rayons)
        if len(rayons) > 28 or len(jauges) > 28:
            raise ValueError('dimension hors domaine')
        def nombre(v):
            if isinstance(v, np.float64):
                v = float(v)
            if type(v) not in (int, float, F):
                raise ValueError('entier, Fraction ou binary64 réel requis')
            q = F(v)
            if max(q.numerator.bit_length(), q.denominator.bit_length()) > 8192:
                raise ValueError('rationnel hors budget')
            return q
        jauges = [(i, nombre(v)) for i, v in jauges]
        rayons = [nombre(r) for r in rayons]
        modele, centre, photographie = depuis_noyau(noyau, jauges)
        if len(rayons) != len(centre) or any(r <= 0 for r in rayons):
            raise ValueError('7*N rayons strictement positifs requis')
        inverse = proposer_inverse(modele, centre)
        return fabriquer_public(modele, centre, rayons, inverse, photographie)
    except np.linalg.LinAlgError as e:
        raise CertificationImpossible("construction de l’inverse approché impossible") from e
    except (ValueError, TypeError, IndexError, OverflowError, ZeroDivisionError) as e:
        raise CertificationImpossible(str(e)) from e


def verifier_certificat(certificat):
    """Vérifie un certificat par un second calcul rationnel, sans solveur flottant.

    Retourne les bornes sous forme de Fraction ; ValueError en cas de refus.
    Le fichier _verification_lineaire.py peut aussi être exécuté directement,
    avec la bibliothèque standard seule et un document JSON en argument.
    """
    return _verifier(certificat)
