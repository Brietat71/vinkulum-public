"""Certificat expérimental sur ker(B.T), sans base dense du complément.

D, R, E, M et B désignent leurs valeurs binary64 interprétées exactement.
Le certificat rétablit leur relation énergétique par le calcul dirigé de
la 0.11.0, puis encadre une correction de rang faible en arithmétique
dirigée. Il n'étend pas à lui seul le contrat de ReductionMaterielle.
"""
from decimal import Decimal
import hashlib
import time

import numpy as np

from vinkulum._ports.certificat_spectral import (
    CertificationSpectraleImpossible, _Intervalles, _float_dirige,
    _certifier_spectral_frais)
from vinkulum._ports.inverse_selectionnee import (
    InverseSelectionnee, BudgetInverseSelectionnee, _entier_positif)


def _ldl(h, ar):
    """LDL sans pivot : ses pivots inférieurs positifs démontrent H SPD."""
    n = len(h)
    l = [[ar.un if i == j else ar.zero for j in range(n)] for i in range(n)]
    d = []
    for i in range(n):
        pivot = h[i][i]
        for k in range(i):
            pivot = ar.sub(pivot, ar.mul(ar.mul(l[i][k], l[i][k]), d[k]))
        if pivot[0] <= 0:
            raise CertificationSpectraleImpossible(
                "rang des contraintes non démontré : pivot H non strictement positif",
                dict(phase="rang_complement", pivot=i, intervalle=list(map(str, pivot)),
                     operations_decimal=ar.operations))
        d.append(pivot)
        for j in range(i+1, n):
            x = h[j][i]
            for k in range(i):
                x = ar.sub(x, ar.mul(ar.mul(l[j][k], l[i][k]), d[k]))
            l[j][i] = ar.div(x, pivot)
    return l, d


def _resoudre_ldl(l, d, rhs, ar):
    n, p = len(l), len(rhs[0])
    y = [[ar.zero for _ in range(p)] for _ in range(n)]
    for i in range(n):
        for c in range(p):
            v = rhs[i][c]
            for k in range(i):
                v = ar.sub(v, ar.mul(l[i][k], y[k][c]))
            y[i][c] = v
    x = [[ar.zero for _ in range(p)] for _ in range(n)]
    for i in range(n-1, -1, -1):
        for c in range(p):
            v = ar.div(y[i][c], d[i])
            for k in range(i+1, n):
                v = ar.sub(v, ar.mul(l[k][i], x[k][c]))
            x[i][c] = v
    return x


def _resoudre_qr(inv, b, ar):
    """U = E R^-1 R^-T E B, avec seulement O(n s) coefficients denses."""
    n, s = b.shape
    e = [ar.point(v) for v in inv.echelles]
    diag = [ar.point(v) for v in inv.diagonale]
    voisins = [[(j, ar.point(v)) for j, v in ligne] for ligne in inv.voisins]
    rhs = [[ar.mul(e[i], ar.point(v)) for v in b[i]] for i in range(n)]
    y = [[ar.zero for _ in range(s)] for _ in range(n)]
    for i in range(n):
        for c in range(s):
            y[i][c] = ar.div(rhs[i][c], diag[i])
        for k, rik in voisins[i]:
            for c in range(s):
                rhs[k][c] = ar.sub(rhs[k][c], ar.mul(rik, y[i][c]))
    z = [[ar.zero for _ in range(s)] for _ in range(n)]
    for i in range(n-1, -1, -1):
        for c in range(s):
            v = y[i][c]
            for k, rik in voisins[i]:
                v = ar.sub(v, ar.mul(rik, z[k][c]))
            z[i][c] = ar.div(v, diag[i])
    return [[ar.mul(e[i], v) for v in z[i]] for i in range(n)]


def _gram_contraintes(b, u, masse, ar):
    n, s = b.shape
    mu = [[ar.zero for _ in range(s)] for _ in range(n)]
    for i in range(n):
        a, fin = masse.indptr[i:i+2]
        for k, valeur in zip(masse.indices[a:fin], masse.data[a:fin]):
            mik = ar.point(valeur)
            for c in range(s):
                mu[i][c] = ar.add(mu[i][c], ar.mul(mik, u[k][c]))
    h = [[ar.zero for _ in range(s)] for _ in range(s)]
    j = [[ar.zero for _ in range(s)] for _ in range(s)]
    # Symétrie exacte par K_Q et M ; un seul triangle suffit à l'encadrement.
    for a in range(s):
        for c in range(a, s):
            hh, jj = ar.zero, ar.zero
            for i in range(n):
                hh = ar.add(hh, ar.mul(ar.point(b[i, a]), u[i][c]))
                jj = ar.add(jj, ar.mul(u[i][a], mu[i][c]))
            h[a][c] = h[c][a] = hh
            j[a][c] = j[c][a] = jj
    return h, j


def certifier_complement(d_i, r, echelles, mii, contraintes, *, precision=80,
                         budget_operations=100_000_000,
                         budget_coefficients=2_000_000,
                         budget_rectangulaire=2_000_000,
                         taille_bloc_masse_max=6):
    """Minore lambda_min de (D_I.T D_I, M_II) restreint à ker(B.T).

    Le rang de B est démontré par une LDL dirigée de H=B.T K_Q^-1 B.
    Une base approchée ou oblique est permise ; B n'est pas identifié à un
    produit M Phi calculé dans une autre arithmétique. Aucune petite valeur
    singulière n'est supprimée silencieusement. Pas de cache mutable fourni
    par l'appelant, pas de K dense ni de base dense du complément.

    Le budget Decimal comprend le certificat initial et la correction.
    Le prétraitement flottant possède séparément ce même budget ; la
    limite n*s+s*s compte les coefficients rectangulaires logiques, et
    non des octets (plusieurs tableaux de cette taille sont utilisés).
    """
    debut = time.perf_counter()
    for nom, a in (("D", d_i), ("R", r), ("E", echelles), ("M", mii), ("B", contraintes)):
        if np.iscomplexobj(a):
            raise ValueError(nom+" réel requis")
    precision = _entier_positif(precision, "precision")
    budget_operations = _entier_positif(budget_operations, "budget_operations")
    budget_rectangulaire = _entier_positif(budget_rectangulaire, "budget_rectangulaire")
    b = np.array(contraintes, dtype=float, copy=True)
    if (b.ndim != 2 or b.shape[0] != r.shape[0]
            or not 0 < b.shape[1] < b.shape[0] or not np.all(np.isfinite(b))):
        raise ValueError("B fini (n,s), avec 0 < s < n requis")
    n, s = b.shape
    coefficients = n*s+s*s
    if coefficients > budget_rectangulaire:
        raise BudgetInverseSelectionnee("budget rectangulaire dépassé",
            dict(phase="complement", coefficients_demandes=coefficients,
                 budget_coefficients=budget_rectangulaire))
    inv = InverseSelectionnee(r, mii, echelles,
        budget_operations=budget_operations, budget_coefficients=budget_coefficients,
        taille_bloc_masse_max=taille_bloc_masse_max)
    base = _certifier_spectral_frais(d_i, inv, precision=precision,
        budget_operations=budget_operations, budget_coefficients=budget_coefficients)
    ar = _Intervalles(precision, budget_operations)
    ar.operations = base["operations_decimal"]
    u = _resoudre_qr(inv, b, ar)
    h, j = _gram_contraintes(b, u, inv.masse, ar)
    l, pivots = _ldl(h, ar)
    hj = _resoudre_ldl(l, pivots, j, ar)
    chi = ar.zero
    for i in range(s):
        chi = ar.add(chi, hj[i][i])
    trace = (Decimal(base["trace_inferieure_decimal"]),
             Decimal(base["trace_superieure_decimal"]))
    complement = ar.sub(trace, chi)
    if complement[1] <= 0:
        raise CertificationSpectraleImpossible("trace complémentaire supérieure non positive",
            dict(phase="trace_complement", intervalle=list(map(str, complement))))
    eta = Decimal(base["eta_superieur_decimal"])
    ar.compter(2)
    lam = ar.bas.divide(ar.bas.subtract(Decimal(1), eta), complement[1])
    valeur = _float_dirige(lam, haut=False)
    if valeur <= 0:
        raise CertificationSpectraleImpossible("minorant complémentaire float non positif",
            dict(phase="conversion_complement"))
    return dict(lambda_min=valeur, lambda_inferieur_decimal=str(lam),
        trace_complement_inferieure_decimal=str(complement[0]),
        trace_complement_superieure_decimal=str(complement[1]),
        correction_inferieure_decimal=str(chi[0]), correction_superieure_decimal=str(chi[1]),
        pivots_h_inferieurs_decimal=[str(p[0]) for p in pivots],
        rang_contraintes=s, dimension_interieure=n, dimension_complement=n-s,
        contraintes_sha256=hashlib.sha256(np.asarray(b, dtype="<f8").tobytes()).hexdigest(),
        contraintes_cibles="noyau de B.T ; valeurs binary64 de B interprétées exactement",
        certificat_total=base, precision_decimal=precision,
        operations_decimal=ar.operations, operations_inverse_float=inv.bilan()["operations_total"],
        coefficients_rectangulaires=coefficients, preparation_s=time.perf_counter()-debut,
        certification_machine=True,
        portee="coercivité sur le complément seulement ; ni Schur global ni réponse certifiés")
