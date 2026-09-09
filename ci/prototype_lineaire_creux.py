"""Certificat par facteurs triangulaires creux et comparaison dirigée.

Les matrices sont des listes de dictionnaires {colonne: Fraction}.
Le système est déjà permuté : toute permutation doit être appliquée à A,
b, x et aux enveloppes avant cet appel. Aucune factorisation n'est crue.
"""
from fractions import Fraction as F


DIMENSION_MAX = 4096
TERMES_MAX = 200_000
PRODUITS_MAX = 2_000_000


def haut(q):
    """Arrondi dyadique supérieur positif, erreur relative <= 2^-127.

    L'exposant reste variable : pas de plancher absolu lié aux unités.
    Le budget d'exposant est une limite de calcul, pas un seuil de preuve.
    """
    q = F(q)
    if q < 0:
        raise ValueError('majoration positive requise')
    if not q:
        return q
    n, d = q.numerator, q.denominator
    k = n.bit_length()-d.bit_length()
    if abs(k) > 16384:
        raise ValueError('budget d’exposant de majoration dépassé')
    shift = 128-k
    if shift >= 0:
        return F(((n << shift)+d-1)//d, 1 << shift)
    den = d << -shift
    return F(((n+den-1)//den) << -shift)


def matrice(a, n, *, positif=False, triangle=None, colonnes=None):
    if type(a) is not list or len(a) != n:
        raise ValueError('dimension de matrice invalide')
    if sum(len(l) for l in a if type(l) is dict) > TERMES_MAX:
        raise ValueError('budget de termes dépassé')
    out = []
    nc = n if colonnes is None else colonnes
    for i, ligne in enumerate(a):
        if type(ligne) is not dict:
            raise ValueError('ligne creuse requise')
        r = {}
        for j, v in ligne.items():
            if type(j) is not int or not 0 <= j < nc:
                raise ValueError('indice hors domaine')
            if type(v) not in (int, F):
                raise ValueError('coefficient rationnel exact requis')
            v = F(v)
            if max(abs(v.numerator).bit_length(), v.denominator.bit_length()) > 2048:
                raise ValueError('budget de coefficient dépassé')
            if positif and v < 0:
                raise ValueError('enveloppe négative')
            if v:
                if triangle == 'L' and j > i or triangle == 'U' and j < i:
                    raise ValueError('facteur non triangulaire')
                r[j] = v
        if triangle and not r.get(i):
            raise ValueError('diagonale de facteur nulle')
        out.append(r)
    return out


def vecteur(v, n, *, positif=False, strict=False):
    if type(v) is not list or len(v) != n:
        raise ValueError('dimension de vecteur invalide')
    out = []
    for x in v:
        if type(x) not in (int, F):
            raise ValueError('composante rationnelle exacte requise')
        x = F(x)
        if max(abs(x.numerator).bit_length(), x.denominator.bit_length()) > 2048:
            raise ValueError('budget de coefficient dépassé')
        if positif and x < 0 or strict and x <= 0:
            raise ValueError('signe de composante invalide')
        out.append(x)
    return out


def defaut_facteurs(a, l, u, poids):
    """Majore |LU-A|w, une ligne à la fois ; aucun tableau n×n."""
    produits = sum(len(u[k]) for ligne in l for k in ligne)
    if produits > PRODUITS_MAX:
        raise ValueError('budget du produit LU dépassé')
    resultat, termes = [], 0
    for i, ligne in enumerate(l):
        row = {j: -v for j, v in a[i].items()}
        for k, v in ligne.items():
            for j, w in u[k].items():
                row[j] = row.get(j, F(0))+v*w
        termes += len(row)
        if termes > TERMES_MAX:
            raise ValueError('budget de remplissage LU-A dépassé')
        resultat.append(haut(sum((abs(v)*poids[j] for j, v in row.items()), F(0))))
    return resultat, produits, termes


def comparaison(l, u, v):
    """Majore |U^-1| |L^-1| v pour v >= 0, avec arrondis supérieurs."""
    n = len(v)
    y, z = [F(0)]*n, [F(0)]*n
    for i in range(n):
        y[i] = haut((v[i]+sum((abs(a)*y[j] for j, a in l[i].items() if j < i), F(0)))/abs(l[i][i]))
    for i in range(n-1, -1, -1):
        z[i] = haut((y[i]+sum((abs(a)*z[j] for j, a in u[i].items() if j > i), F(0)))/abs(u[i][i]))
    return z


def certifier(a, b, x, l, u, *, poids=None, delta_a=None, delta_b=None):
    """Toutes les A+ΔA, b+Δb de l'enveloppe : unicité et erreur en avant.

    ΔA est bornée composante par composante par delta_a, Δb par delta_b.
    Les coefficients absents ont une incertitude nulle. Les comparaisons
    positives sont conservatrices ; un refus ne prouve pas la singularité.
    """
    if type(a) is not list or not 1 <= len(a) <= DIMENSION_MAX:
        raise ValueError('dimension hors budget')
    n = len(a)
    a, l, u = matrice(a, n), matrice(l, n, triangle='L'), matrice(u, n, triangle='U')
    b, x = vecteur(b, n), vecteur(x, n)
    w = vecteur([1]*n if poids is None else poids, n, strict=True)
    da = matrice([{} for _ in range(n)] if delta_a is None else delta_a, n, positif=True)
    db = vecteur([0]*n if delta_b is None else delta_b, n, positif=True)
    e, produits, termes = defaut_facteurs(a, l, u, w)
    e = [haut(e[i]+sum((v*w[j] for j, v in da[i].items()), F(0))) for i in range(n)]
    z = comparaison(l, u, e)
    q = max(z[i]/w[i] for i in range(n))
    if q >= 1:
        raise ValueError('contraction uniforme non établie')
    r = [haut(abs(sum((v*x[j] for j, v in a[i].items()), F(0))-b[i])
              +sum((v*abs(x[j]) for j, v in da[i].items()), F(0))+db[i]) for i in range(n)]
    s = comparaison(l, u, r)
    beta = haut(max(s[i]/w[i] for i in range(n))/(1-q))
    bornes = [min(haut(w[i]*beta), haut(s[i]+z[i]*beta)) for i in range(n)]
    return {'contraction': q, 'borne_ponderee': beta, 'bornes_composantes': bornes,
            'produits_facteurs': produits, 'termes_defaut': termes}


def kkt(m, c):
    n = len(m)
    a = [row.copy() for row in m]+[row.copy() for row in c]
    for i, row in enumerate(c):
        for j, v in row.items():
            a[j][n+i] = v
    return a


def certifier_dependances_structurelles(m, c, t, base, force, d, x, l, u, *,
                                        delta_m=None, delta_c=None, delta_t=None,
                                        delta_force=None, delta_d=None, poids=None):
    """Famille G'=T'C', second membre T'd', avec T'[base,:]=I exact.

    Ce contrat est une paramétrisation fournie, jamais une déduction à
    partir de lignes proches. Les incertitudes indépendantes sur G ne sont
    pas admises. La preuve couvre toutes les contraintes de la famille.
    x contient l'accélération et les multiplicateurs réduits proposés.
    Les facteurs L,U portent sur le KKT réduit, dans l'ordre fourni ici.
    """
    if type(m) is not list or type(c) is not list or type(t) is not list:
        raise ValueError('matrices structurelles requises')
    n, r, nb = len(m), len(c), len(t)
    if not 1 <= r <= n or n+r > DIMENSION_MAX or not r <= nb <= DIMENSION_MAX:
        raise ValueError('dimensions structurelles hors domaine')
    m, c, t = matrice(m, n), matrice(c, r, colonnes=n), matrice(t, nb, colonnes=r)
    dm = matrice([{} for _ in range(n)] if delta_m is None else delta_m, n, positif=True)
    dc = matrice([{} for _ in range(r)] if delta_c is None else delta_c, r, colonnes=n, positif=True)
    dt = matrice([{} for _ in range(nb)] if delta_t is None else delta_t, nb, colonnes=r, positif=True)
    if (type(base) is not list or len(base) != r or any(type(i) is not int or not 0 <= i < nb for i in base)
            or len(set(base)) != r):
        raise ValueError('indices de base invalides')
    for j, i in enumerate(base):
        if t[i] != {j: F(1)} or dt[i]:
            raise ValueError('identité structurelle des lignes de base non établie')
    f, d = vecteur(force, n), vecteur(d, r)
    df = vecteur([0]*n if delta_force is None else delta_force, n, positif=True)
    dd = vecteur([0]*r if delta_d is None else delta_d, r, positif=True)
    x = vecteur(x, n+r)
    preuve = certifier(kkt(m, c), f+d, x, l, u, poids=poids, delta_a=kkt(dm, dc), delta_b=df+dd)
    errors = preuve['bornes_composantes']
    reaction, borne = [F(0)]*n, [F(0)]*n
    for j in range(r):
        for i in c[j].keys() | dc[j].keys():
            valeur, incertitude = c[j].get(i, F(0)), dc[j].get(i, F(0))
            reaction[i] += valeur*x[n+j]
            borne[i] += abs(valeur)*errors[n+j]+incertitude*(abs(x[n+j])+errors[n+j])
    return {'systeme_reduit': preuve, 'contraintes_originales': nb,
            'rang_uniforme': r, 'reaction_proposee': reaction,
            'bornes_reaction': list(map(haut, borne))}
