"""Minoration spectrale dirigée pour les valeurs binary64 du D fourni.

Ce certificat porte sur K = D.T D en arithmétique exacte, et sur la
masse fournie, pas sur le produit flottant D.T@D ni sur un continuum.
L'inverse sélectionnée et les différences de Gram sont encadrées par
Decimal avec arrondis dirigés, sans modifier le contexte global.
"""
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
import math
import time

import numpy as np
from scipy.sparse import csr_matrix

from .inverse_selectionnee import (BudgetInverseSelectionnee, InverseSelectionnee,
                                  _entier_positif, _verifier_canonique)


class CertificationSpectraleImpossible(RuntimeError):
    def __init__(self, message, bilan):
        super().__init__(message)
        self.bilan = dict(bilan)


class _Intervalles:
    """Opérations élémentaires arrondies vers l'extérieur, contextes privés."""
    def __init__(self, precision, budget):
        self.bas = Context(prec=precision, rounding=ROUND_FLOOR)
        self.haut = Context(prec=precision, rounding=ROUND_CEILING)
        self.budget = budget
        self.operations = 0
        self.zero = (Decimal(0), Decimal(0))
        self.un = (Decimal(1), Decimal(1))

    def compter(self, n=1):
        if self.operations+n > self.budget:
            raise BudgetInverseSelectionnee(
                "budget d'opérations Decimal dépassé",
                dict(phase="certification", operations_demandees=self.operations+n,
                     budget_operations=self.budget))
        self.operations += n

    @staticmethod
    def point(x):
        p = Decimal.from_float(float(x))
        if not p.is_finite():
            raise ValueError("donnée non finie pour le certificat")
        return (p, p)

    def add(self, a, b):
        self.compter(2)
        return self.bas.add(a[0], b[0]), self.haut.add(a[1], b[1])

    def sub(self, a, b):
        self.compter(2)
        return self.bas.subtract(a[0], b[1]), self.haut.subtract(a[1], b[0])

    @staticmethod
    def neg(a):
        # copy_negate, contrairement à l'opérateur unaire, est exact et
        # indépendant de la précision du contexte Decimal ambiant.
        return a[1].copy_negate(), a[0].copy_negate()

    def mul(self, a, b):
        if a[0] == a[1] and b[0] == b[1]:
            self.compter(2)
            return self.bas.multiply(a[0], b[0]), self.haut.multiply(a[1], b[1])
        self.compter(8)
        return (min(self.bas.multiply(x, y) for x in a for y in b),
                max(self.haut.multiply(x, y) for x in a for y in b))

    def div(self, a, b):
        if b[0] <= 0 <= b[1]:
            raise CertificationSpectraleImpossible(
                "division par un intervalle contenant zéro", dict(phase="intervalles"))
        self.compter(2)
        inverse = (self.bas.divide(Decimal(1), b[1]),
                   self.haut.divide(Decimal(1), b[0]))
        return self.mul(a, inverse)

    def racine_superieure(self, x):
        if x < 0:
            raise ArithmeticError("racine d'un majorant négatif")
        self.compter(2)
        # Decimal.sqrt est correctement arrondi au plus proche, même si
        # le contexte porte ROUND_CEILING. Un successeur donne donc un
        # majorant ; ne pas supposer que sqrt respecte ROUND_CEILING.
        return self.haut.next_plus(self.haut.sqrt(x)) if x else Decimal(0)


def _float_dirige(x, haut):
    valeur = float(x)
    if not math.isfinite(valeur):
        raise CertificationSpectraleImpossible(
            "constante certifiée non représentable en float fini", dict(phase="conversion"))
    point = Decimal.from_float(valeur)
    if (haut and point < x) or (not haut and point > x):
        valeur = float(np.nextafter(valeur, np.inf if haut else -np.inf))
    if not math.isfinite(valeur):
        raise CertificationSpectraleImpossible(
            "conversion dirigée hors domaine float", dict(phase="conversion"))
    return valeur


def _verifier_masse(inv, ar):
    """LDL dirigée de chaque petit bloc ; tous les pivots doivent être >0."""
    from scipy.sparse.csgraph import connected_components
    _, labels = connected_components(inv.masse, directed=False)
    groupes = {}
    for i, label in enumerate(labels):
        groupes.setdefault(int(label), []).append(i)
    for ids in groupes.values():
        b = len(ids)
        l = [[ar.zero for _ in range(b)] for _ in range(b)]
        pivots = []
        for i in range(b):
            pivot = ar.point(inv.masse[ids[i], ids[i]])
            for k in range(i):
                pivot = ar.sub(pivot, ar.mul(ar.mul(l[i][k], l[i][k]), pivots[k]))
            if pivot[0] <= 0:
                raise CertificationSpectraleImpossible(
                    "positivité du bloc de masse non démontrée par LDL dirigée",
                    dict(phase="masse", bloc=ids, pivot_inferieur=str(pivot[0])))
            pivots.append(pivot)
            for j in range(i+1, b):
                terme = ar.point(inv.masse[ids[j], ids[i]])
                for k in range(i):
                    terme = ar.sub(terme, ar.mul(ar.mul(l[j][k], l[i][k]), pivots[k]))
                l[j][i] = ar.div(terme, pivot)


def _certifier_spectral_frais(d_interieur, inverse, *, precision=80,
                            budget_operations=100_000_000,
                            budget_coefficients=2_000_000):
    """Certifie lambda_min(D_I.T D_I, M_II) pour les données float exactes.

    `inverse` fournit le R, E et M d'une InverseSelectionnee. La relation
    entre R et D_I n'est jamais présumée : leur différence énergétique
    est encadrée explicitement. Si eta >= 1, si un intervalle ne conclut
    pas ou si le budget est dépassé, la fonction refuse le certificat.

    Les budgets Decimal sont distincts du calcul flottant préalable.
    Le budget de coefficients limite les paires du Gram d'écart ; le
    motif sélectionné est contrôlé par le budget de l'objet `inverse`.
    """
    debut = time.perf_counter()
    if not isinstance(inverse, InverseSelectionnee):
        raise TypeError("un objet InverseSelectionnee est requis")
    precision = _entier_positif(precision, "precision")
    if precision < 16:
        raise ValueError("au moins 16 chiffres Decimal requis")
    budget_operations = _entier_positif(budget_operations, "budget_operations")
    budget_coefficients = _entier_positif(budget_coefficients, "budget_coefficients")
    _verifier_canonique(d_interieur, "D")
    d = csr_matrix(d_interieur, dtype=float, copy=True)
    _verifier_canonique(d, "D converti")
    if d.shape[1] != inverse.n or not np.all(np.isfinite(d.data)):
        raise ValueError("D intérieur fini avec autant de colonnes que R requis")
    ar = _Intervalles(precision, budget_operations)
    _verifier_masse(inverse, ar)
    e = [ar.point(x) for x in inverse.echelles]
    rdiag = [ar.point(x) for x in inverse.diagonale]
    voisins = [[(j, ar.point(v)) for j, v in row] for row in inverse.voisins]
    valeurs = [dict() for _ in range(inverse.n)]

    def lire(i, j):
        if i > j:
            i, j = j, i
        return valeurs[i][j]

    for i in range(inverse.n-1, -1, -1):
        for j in sorted(j for j in inverse.valeurs[i] if j > i):
            somme = ar.zero
            for k, rik in voisins[i]:
                somme = ar.add(somme, ar.mul(rik, lire(k, j)))
            valeurs[i][j] = ar.neg(ar.div(somme, rdiag[i]))
        somme = ar.zero
        for k, rik in voisins[i]:
            somme = ar.add(somme, ar.mul(rik, valeurs[i][k]))
        valeurs[i][i] = ar.sub(ar.div(ar.un, ar.mul(rdiag[i], rdiag[i])),
                              ar.div(somme, rdiag[i]))
    diagonales = []
    for i in range(inverse.n):
        diagonale = ar.mul(ar.mul(e[i], e[i]), valeurs[i][i])
        if diagonale[1] <= 0:
            raise CertificationSpectraleImpossible(
                "diagonale inverse supérieure non positive", dict(phase="inverse"))
        diagonales.append(diagonale)
    trace = ar.zero
    for i in range(inverse.n):
        a, b = inverse.masse.indptr[i:i+2]
        for j, mij in zip(inverse.masse.indices[a:b], inverse.masse.data[a:b]):
            if j < i:
                continue
            terme = ar.mul(ar.point(mij), ar.mul(ar.mul(e[i], e[j]), lire(i, int(j))))
            trace = ar.add(trace, terme if i == j else ar.add(terme, terme))
    if trace[1] <= 0:
        raise CertificationSpectraleImpossible("trace supérieure non positive", dict(phase="trace"))

    delta = {}
    produits_gram = 0

    def ajouter_gram(ids, vals, signe):
        nonlocal produits_gram
        for a, i in enumerate(ids):
            for b in range(a, len(ids)):
                j = ids[b]
                paire = (int(min(i, j)), int(max(i, j)))
                if paire not in delta and len(delta)+1 > budget_coefficients:
                    raise BudgetInverseSelectionnee(
                        "budget de coefficients du Gram d'écart dépassé",
                        dict(phase="gram_ecart", coefficients_demandes=len(delta)+1,
                             budget_coefficients=budget_coefficients))
                terme = ar.mul(vals[a], vals[b])
                ancien = delta.get(paire, ar.zero)
                delta[paire] = ar.add(ancien, terme) if signe > 0 else ar.sub(ancien, terme)
                produits_gram += 1

    # Aucun produit sparse flottant de Gram n'intervient dans le certificat.
    for row in range(d.shape[0]):
        a, b = d.indptr[row:row+2]
        ajouter_gram(d.indices[a:b], [ar.point(x) for x in d.data[a:b]], 1)
    for i in range(inverse.n):
        a, b = inverse.r.indptr[i:i+2]
        ids = inverse.r.indices[a:b]
        vals = [ar.div(ar.point(v), e[j]) for j, v in zip(ids, inverse.r.data[a:b])]
        ajouter_gram(ids, vals, -1)
    eta = Decimal(0)
    for (i, j), intervalle in delta.items():
        absolu = max(intervalle[0].copy_abs(), intervalle[1].copy_abs())
        ar.compter(3)
        racine = ar.racine_superieure(ar.haut.multiply(diagonales[i][1], diagonales[j][1]))
        contribution = ar.haut.multiply(absolu, racine)
        if i != j:
            ar.compter()
            contribution = ar.haut.multiply(Decimal(2), contribution)
        eta = ar.haut.add(eta, contribution)
    bilan = dict(precision_decimal=precision, operations_decimal=ar.operations,
                 coefficients_gram_ecart=len(delta), produits_gram=produits_gram,
                 coefficients_inverse=inverse.coefficients,
                 trace_inferieure_decimal=str(trace[0]), trace_superieure_decimal=str(trace[1]),
                 eta_superieur_decimal=str(eta),
                 modele="D_I.T D_I et M_II : valeurs binary64 d'entrée interprétées exactement",
                 portee="minoration spectrale seule ; pas de certificat du modèle physique continu",
                 preparation_s=time.perf_counter()-debut)
    if eta >= 1:
        raise CertificationSpectraleImpossible(
            "écart relatif énergétique non inférieur à 1 : coercivité non démontrée", bilan)
    ar.compter(2)
    lambda_bas = ar.bas.divide(ar.bas.subtract(Decimal(1), eta), trace[1])
    valeur = _float_dirige(lambda_bas, haut=False)
    if valeur <= 0:
        raise CertificationSpectraleImpossible("minoration spectrale float non positive", bilan)
    bilan.update(lambda_min=valeur, lambda_inferieur_decimal=str(lambda_bas),
                 eta_superieur=_float_dirige(eta, haut=True),
                 trace_superieure=_float_dirige(trace[1], haut=True),
                 operations_decimal=ar.operations, preparation_s=time.perf_counter()-debut,
                 minoration_D_original_etablie=True, certification_machine=True)
    return bilan


def certifier_spectral(d_interieur, inverse, *, precision=80,
                      budget_operations=100_000_000,
                      budget_coefficients=2_000_000):
    """Certifie depuis les R, M et E actuels d'un objet éventuellement mutable.

    La sélection est reconstruite : aucun voisin, diagonale ou résultat
    numérique mémorisé avant une mutation de l'objet ne participe à la
    preuve. Le budget d'opérations s'applique séparément à cette nouvelle
    sélection et au calcul Decimal.
    """
    if not isinstance(inverse, InverseSelectionnee):
        raise TypeError("un objet InverseSelectionnee est requis")
    frais = InverseSelectionnee(
        inverse.r, inverse.masse, inverse.echelles,
        budget_operations=budget_operations, budget_coefficients=budget_coefficients,
        taille_bloc_masse_max=inverse.bloc_masse_max)
    resultat = _certifier_spectral_frais(
        d_interieur, frais, precision=precision, budget_operations=budget_operations,
        budget_coefficients=budget_coefficients)
    resultat["operations_inverse_float"] = frais.bilan()["operations_total"]
    return resultat


def certifier_spectre(d_i, r, echelles, mii, *, precision=80,
                     budget=100_000_000, budget_coefficients=2_000_000,
                     taille_bloc_masse_max=6):
    """Adaptateur autonome depuis D_I, R, E et M_II.

    `budget` est appliqué séparément à la sélection flottante et au
    certificat Decimal ; leurs compteurs sont tous deux renvoyés.
    """
    inverse = InverseSelectionnee(
        r, mii, echelles, budget_operations=budget,
        budget_coefficients=budget_coefficients,
        taille_bloc_masse_max=taille_bloc_masse_max)
    resultat = _certifier_spectral_frais(
        d_i, inverse, precision=precision, budget_operations=budget,
        budget_coefficients=budget_coefficients)
    resultat.update(eta_upper=resultat["eta_superieur"],
                    trace_upper=resultat["trace_superieure"], precision=precision,
                    operations_inverse_float=inverse.bilan()["operations_total"])
    return resultat
