"""Oracles de Schur en décimales de haute précision, hors chronométrage.

Deux problèmes de référence sont distingués explicitement :

* ``schur_console`` utilise les paramètres physiques de la fixture, sans
  lire K ni D assemblés. Chaque flottant d'entrée représente exactement sa
  valeur binaire (Decimal.from_float). Les rigidités EA/GA/EI/GJ sont celles
  des métadonnées. La masse est recomposée depuis rho, r, L et la constante
  binaire ``math.pi`` commune à la fixture ; les produits et L/n sont ensuite
  évalués en haute précision, sans les arrondis intermédiaires du noyau.
* ``schur_facteur`` résout exactement le problème défini par les coefficients
  flottants fournis de D et M : D.T D est assemblé en Decimal, puis condensé.
  Ce problème peut différer de celui d'une raideur K assemblée en flottants.

Les deux oracles utilisent O(n) stockage et opérations pour une largeur de
bloc fixée <=6. Aucune matrice globale dense n'est formée. Ils ne fournissent
pas de certificat d'arrondi ; augmenter dps permet de contrôler la convergence
de la référence. L'élimination exige des pivots de blocs non singuliers,
condition satisfaite sous la première fréquence intérieure à ports bloqués.
"""
from decimal import Context, Decimal, ROUND_HALF_EVEN, localcontext
import math
import operator

import numpy as np
from scipy.sparse import csr_matrix, issparse


_ZERO, _UN = Decimal(0), Decimal(1)


def _precision(dps):
    if isinstance(dps, (bool, np.bool_)):
        raise ValueError("dps doit être un entier supérieur ou égal à 30")
    try:
        dps = operator.index(dps)
    except TypeError as exc:
        raise ValueError("dps doit être un entier supérieur ou égal à 30") from exc
    if dps < 30:
        raise ValueError("dps doit être un entier supérieur ou égal à 30")
    return Context(prec=dps, rounding=ROUND_HALF_EVEN)


def _decimal(x, nom, positif=False):
    if isinstance(x, Decimal):
        valeur = x
    elif isinstance(x, (int, np.integer)):
        valeur = Decimal(int(x))
    else:
        valeur = Decimal.from_float(float(x))
    if not valeur.is_finite() or (positif and valeur <= 0):
        raise ValueError(nom+" doit être fini"+(" et strictement positif" if positif else ""))
    return valeur


def _frequence(omega):
    omega = _decimal(omega, "omega")
    if omega < 0:
        raise ValueError("omega doit être positif ou nul")
    return omega*omega


def _zeros(n):
    return [[_ZERO for _ in range(n)] for _ in range(n)]


def _resultat(a, retour_decimal):
    if retour_decimal:
        return tuple(tuple(ligne) for ligne in a)
    result = np.array([[float(x) for x in ligne] for ligne in a])
    if not np.all(np.isfinite(result)):
        raise OverflowError("le Schur de référence ne tient pas dans des flottants doubles")
    return result


def _resout_bloc(a, b):
    """Élimination pivotée d'un bloc de taille au plus 6, plusieurs RHS."""
    n = len(a)
    aug = [list(a[j])+list(b[j]) for j in range(n)]
    for j in range(n):
        pivot = max(range(j, n), key=lambda k: abs(aug[k][j]))
        if aug[pivot][j] == 0:
            raise ValueError("pivot de bloc nul dans l'oracle de Schur")
        aug[j], aug[pivot] = aug[pivot], aug[j]
        for k in range(j+1, n):
            if aug[k][j] == 0:
                continue
            facteur = aug[k][j]/aug[j][j]
            aug[k][j] = _ZERO
            for col in range(j+1, 2*n):
                aug[k][col] -= facteur*aug[j][col]
    result = _zeros(n)
    for j in range(n-1, -1, -1):
        for col in range(n):
            result[j][col] = (aug[j][n+col]
                              -sum((aug[j][k]*result[k][col] for k in range(j+1, n)), _ZERO))/aug[j][j]
    return result


def _retire_bloc(pivot, liaison, suivant):
    """suivant - liaison.T pivot^-1 liaison, en conservant la symétrie."""
    n = len(pivot)
    resolu = _resout_bloc(pivot, liaison)
    result = _zeros(n)
    for i in range(n):
        for j in range(i, n):
            cij = sum((liaison[k][i]*resolu[k][j] for k in range(n)), _ZERO)
            cji = sum((liaison[k][j]*resolu[k][i] for k in range(n)), _ZERO)
            result[i][j] = result[j][i] = suivant[i][j]-(cij+cji)/2
    return result


def _scalaire(n, raideur, masse, z):
    courant = 2*raideur-z*masse
    for j in range(1, n):
        if courant == 0:
            raise ValueError("pivot scalaire nul dans l'oracle de Schur")
        dernier = j == n-1
        diagonal = raideur if dernier else 2*raideur
        poids = Decimal("0.5") if dernier else _UN
        courant = diagonal-z*masse*poids-raideur*raideur/courant
    return courant


def _flexion(n, pas, ga, ei, masse, inertie, z):
    ceff = _UN/(_UN/ga+pas*pas/(12*ei))
    a, b = ceff/pas, ceff/2
    c, d = ei/pas+ceff*pas/4, -ei/pas+ceff*pas/4
    liaison = [[-a, b], [-b, d]]
    courant = [[2*a-z*masse, _ZERO], [_ZERO, 2*c-z*inertie]]
    for j in range(1, n):
        if j == n-1:
            suivant = [[a-z*masse/2, -b], [-b, c-z*inertie/2]]
        else:
            suivant = [[2*a-z*masse, _ZERO], [_ZERO, 2*c-z*inertie]]
        # Formule explicite de l'inverse 2x2, sans matrice de dimension n.
        determinant = courant[0][0]*courant[1][1]-courant[0][1]*courant[1][0]
        if determinant == 0:
            raise ValueError("pivot de flexion nul dans l'oracle de Schur")
        resolu = [[(courant[1][1]*liaison[0][j]-courant[0][1]*liaison[1][j])/determinant
                   for j in range(2)],
                  [(-courant[1][0]*liaison[0][j]+courant[0][0]*liaison[1][j])/determinant
                   for j in range(2)]]
        courant = _zeros(2)
        for i in range(2):
            for j in range(i, 2):
                cij = sum((liaison[k][i]*resolu[k][j] for k in range(2)), _ZERO)
                cji = sum((liaison[k][j]*resolu[k][i] for k in range(2)), _ZERO)
                courant[i][j] = courant[j][i] = suivant[i][j]-(cij+cji)/2
    return courant


def schur_console(metadata, omega, dps=70, retour_decimal=False):
    """Schur physique 6x6 de la console homogène circulaire intégrée.

    L'intérieur contient les n-1 premiers corps libres ; le port est le
    dernier corps (tx,ty,tz,rx,ry,rz). Les masses et inerties de ce dernier
    corps ont poids 1/2. Le corps racine fixé est retiré. Le calcul comprend
    deux récurrences scalaires et une récurrence 2x2, commune aux deux plans
    de flexion. ``retour_decimal=True`` conserve les décimales calculées.
    """
    with localcontext(_precision(dps)):
        if (metadata.get("famille") != "console"
                or metadata.get("formulation_poutre") != "integree"
                or metadata.get("coefficient_cisaillement") != 1.0
                or metadata.get("poids_extremites_masse_et_inertie") != .5
                or metadata.get("ordre_ddl_par_corps") != ["tx", "ty", "tz", "rx", "ry", "rz"]):
            raise ValueError("oracle réservé à la console homogène circulaire intégrée")
        try:
            n = operator.index(metadata["n"])
        except TypeError as exc:
            raise ValueError("n doit être un entier supérieur ou égal à 2") from exc
        if isinstance(metadata["n"], (bool, np.bool_)) or n < 2:
            raise ValueError("n doit être un entier supérieur ou égal à 2")
        if metadata.get("ports") != list(range(6*(n-1), 6*n)):
            raise ValueError("ports différents du dernier corps physique")
        noms = ("longueur_m", "EA_N", "GA_N", "EI_N_m2", "GJ_N_m2",
                "rayon_m", "masse_volumique_kg_par_m3")
        longueur, ea, ga, ei, gj, rayon, rho = [_decimal(metadata[x], x, True) for x in noms]
        pas, pi = longueur/n, Decimal.from_float(math.pi)
        masse = rho*pi*rayon*rayon*pas
        inertie = rho*pi*rayon**4*pas/4
        z = _frequence(omega)
        axial = _scalaire(n, ea/pas, masse, z)
        torsion = _scalaire(n, gj/pas, 2*inertie, z)
        plan = _flexion(n, pas, ga, ei, masse, inertie, z)
        result = _zeros(6)
        result[0][0], result[3][3] = axial, torsion
        for i, j in ((0, 0), (0, 1), (1, 0), (1, 1)):
            result[(1, 5)[i]][(1, 5)[j]] = plan[i][j]
            result[(2, 4)[i]][(2, 4)[j]] = plan[i][j] if i == j else -plan[i][j]
        return _resultat(result, retour_decimal)


def _creux_reel(a, nom):
    # Refuser les doublons avant une conversion COO->CSR susceptible de les
    # sommer en doubles : l'oracle doit conserver les coefficients fournis.
    if issparse(a) and not getattr(a, "has_canonical_format", a.format in ("dia", "dok")):
        raise ValueError(nom+" doit avoir un stockage creux canonique")
    a = csr_matrix(a, copy=True)
    if a.dtype.kind == "c" or not np.all(np.isfinite(a.data)):
        raise ValueError(nom+" doit être réel et fini")
    if not a.has_canonical_format:
        raise ValueError(nom+" doit avoir un stockage creux canonique")
    return a


def schur_facteur(d, m, interieur, interface, omega, dps=70, retour_decimal=False):
    """Schur du problème exact D_float.T D_float-omega² M_float.

    D est canonique et chaque ligne ne touche qu'un bloc physique ou deux
    blocs voisins. La largeur de bloc est le nombre de ports, de 1 à 6.
    L'interface contient le dernier bloc et l'intérieur tous les précédents ;
    l'ordre demandé des ports est respecté. M est diagonale, strictement
    positive. Les blocs de D.T D sont construits en Decimal, pas en double.
    Le coût est O(n p³+nnz(D) p), le stockage O(n p²), pour p<=6.

    Un schéma creux quelconque n'est pas accepté implicitement : une ligne
    reliant des blocs non voisins provoque un refus explicite.
    """
    with localcontext(_precision(dps)):
        d, m = _creux_reel(d, "D"), _creux_reel(m, "M")
        ii, ss = np.asarray(interieur), np.asarray(interface)
        if ii.ndim != 1 or ss.ndim != 1 or ii.dtype.kind not in "iu" or ss.dtype.kind not in "iu":
            raise ValueError("indices d'intérieur et de port entiers requis")
        total, p = d.shape[1], len(ss)
        if (not 1 <= p <= 6 or total % p or total < 2*p
                or m.shape != (total, total)
                or not np.array_equal(np.sort(ii), np.arange(total-p))
                or not np.array_equal(np.sort(ss), np.arange(total-p, total))):
            raise ValueError("partition incompatible avec les blocs physiques contigus et le dernier port")
        n = total//p
        diagonaux = [_zeros(p) for _ in range(n)]
        liaisons = [_zeros(p) for _ in range(n-1)]
        for row in range(d.shape[0]):
            debut, fin = d.indptr[row:row+2]
            indices, data = d.indices[debut:fin], d.data[debut:fin]
            actifs = [(int(j), _decimal(v, "D")) for j, v in zip(indices, data) if v != 0]
            if actifs and actifs[-1][0]//p-actifs[0][0]//p > 1:
                raise ValueError("une ligne de D relie des blocs non voisins")
            for ja, va in actifs:
                bloc_a, ca = divmod(ja, p)
                for jb, vb in actifs:
                    bloc_b, cb = divmod(jb, p)
                    if bloc_a == bloc_b:
                        diagonaux[bloc_a][ca][cb] += va*vb
                    elif bloc_b == bloc_a+1:
                        liaisons[bloc_a][ca][cb] += va*vb
        z = _frequence(omega)
        for row in range(total):
            debut, fin = m.indptr[row:row+2]
            if any(j != row and v != 0 for j, v in zip(m.indices[debut:fin], m.data[debut:fin])):
                raise ValueError("M doit être diagonale dans cet oracle")
            masse = _decimal(m[row, row], "masse diagonale", True)
            bloc, col = divmod(row, p)
            diagonaux[bloc][col][col] -= z*masse
        courant = diagonaux[0]
        for j in range(1, n):
            courant = _retire_bloc(courant, liaisons[j-1], diagonaux[j])
        ordre = [int(j)-(total-p) for j in ss]
        courant = [[courant[i][j] for j in ordre] for i in ordre]
        return _resultat(courant, retour_decimal)
