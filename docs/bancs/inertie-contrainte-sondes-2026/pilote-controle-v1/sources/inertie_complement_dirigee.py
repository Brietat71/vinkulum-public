"""Coercivité du complément par inertie creuse à arrondis dirigés.

Le KKT vise exactement D.T D - gamma M et le B binary64 fournis. Les
pivots intervalles constituent des congruences, jamais les signes d'un
facteur flottant pris pour une preuve. Un budget ou signe indécidable
provoque un refus explicite. Prototype de recherche, hors API publique.
"""
from decimal import Decimal
import hashlib
import math
import time

import numpy as np
from scipy.sparse import csr_matrix

from vinkulum._ports.certificat_spectral import _Intervalles
from vinkulum._ports.inverse_selectionnee import _entier_positif, _verifier_canonique


class InertieImpossible(RuntimeError):
    def __init__(self, message, bilan):
        super().__init__(message)
        self.bilan = dict(bilan)
        self.diagnostic = self.bilan


class IntervallesSignes(_Intervalles):
    """Produits de signes séparés : deux arrondis suffisent, sinon huit."""
    def mul(self, a, b):
        if a == self.zero or b == self.zero:
            return self.zero
        if a[0] >= 0:
            if b[0] >= 0:
                lo, hi = (a[0], b[0]), (a[1], b[1])
            elif b[1] <= 0:
                lo, hi = (a[1], b[0]), (a[0], b[1])
            else:
                lo, hi = (a[1], b[0]), (a[1], b[1])
        elif a[1] <= 0:
            if b[0] >= 0:
                lo, hi = (a[0], b[1]), (a[1], b[0])
            elif b[1] <= 0:
                lo, hi = (a[1], b[1]), (a[0], b[0])
            else:
                lo, hi = (a[0], b[1]), (a[0], b[0])
        elif b[0] >= 0:
            lo, hi = (a[0], b[1]), (a[1], b[1])
        elif b[1] <= 0:
            lo, hi = (a[1], b[0]), (a[0], b[0])
        else:
            return super().mul(a, b)
        self.compter(2)
        return self.bas.multiply(*lo), self.haut.multiply(*hi)

    def carre(self, a):
        self.compter(2 if a[0] >= 0 or a[1] <= 0 else 3)
        if a[0] >= 0:
            return self.bas.multiply(a[0], a[0]), self.haut.multiply(a[1], a[1])
        if a[1] <= 0:
            return self.bas.multiply(a[1], a[1]), self.haut.multiply(a[0], a[0])
        return Decimal(0), max(self.haut.multiply(a[0], a[0]), self.haut.multiply(a[1], a[1]))


class MatriceIntervalles:
    """Dictionnaires symétriques ; seuls les zéros exactement prouvés sortent."""
    def __init__(self, n, budget):
        if n > budget:
            raise InertieImpossible("budget de dimension dépassé", dict(dimension=n, budget=budget))
        self.lignes = [dict() for _ in range(n)]
        self.budget = budget
        self.coefficients = self.pic = 0
        self.crees = 0

    def get(self, i, j, zero):
        return self.lignes[i].get(j, zero)

    def set(self, i, j, v):
        ancien = j in self.lignes[i]
        if v[0] == 0 == v[1]:
            if ancien:
                self.lignes[i].pop(j)
                if i != j:
                    self.lignes[j].pop(i)
                self.coefficients -= 1
            return
        if not all(x.is_finite() for x in v) or v[0] > v[1]:
            raise InertieImpossible("intervalle non fini ou inversé", dict(i=i, j=j))
        if not ancien:
            if self.coefficients >= self.budget:
                raise InertieImpossible("budget de coefficients dépassé",
                    dict(coefficients=self.coefficients+1, budget=self.budget))
            self.coefficients += 1
            self.crees += 1
            self.pic = max(self.pic, self.coefficients)
        self.lignes[i][j] = v
        if i != j:
            self.lignes[j][i] = v

    def retirer(self, ids):
        for i in ids:
            for j in list(self.lignes[i]):
                self.set(i, j, (Decimal(0), Decimal(0)))


def _signature_bloc(a, b, c, ar):
    det = ar.sub(ar.mul(a, c), ar.carre(b))
    if det[1] < 0:
        return (1, 1), det
    if det[0] > 0:
        tr = ar.add(a, c)
        if tr[0] > 0:
            return (2, 0), det
        if tr[1] < 0:
            return (0, 2), det
    return None, det


def eliminer(matrice, ar, ordre=None, *, positifs_seulement=False):
    """Congruences intervalles ; aucune signature inférée du centre."""
    n = len(matrice.lignes)
    ordre = list(range(n)) if ordre is None else list(map(int, ordre))
    rang = {v: i for i, v in enumerate(ordre)}
    if len(rang) != n or set(rang) != set(range(n)):
        raise ValueError("permutation complète requise")
    actifs = set(ordre)
    pivots, signature, largeur = [], [0, 0, 0], 0
    for k in ordre:
        if k not in actifs:
            continue
        a = matrice.get(k, k, ar.zero)
        voisins = sorted((j for j in matrice.lignes[k] if j != k), key=rang.__getitem__)
        if a[0] > 0 or a[1] < 0:
            if positifs_seulement and a[1] < 0:
                raise InertieImpossible("masse non définie positive", dict(pivot=k, intervalle=list(map(str, a))))
            bloc, sig = [k], (1, 0) if a[0] > 0 else (0, 1)
            inverse = ar.div(ar.un, a)
            v = {i: matrice.lignes[k][i] for i in voisins}
            u = {i: ar.mul(v[i], inverse) for i in voisins}
            preuve = dict(indices=bloc, signature=list(sig), diagonal=list(map(str, a)))

            def correction(i, j):
                return ar.mul(ar.carre(v[i]), inverse) if i == j else ar.mul(u[i], v[j])
        else:
            choisi = None
            for l in voisins:
                b, c = matrice.get(k, l, ar.zero), matrice.get(l, l, ar.zero)
                sig, det = _signature_bloc(a, b, c, ar)
                if sig is not None:
                    choisi = (l, b, c, sig, det)
                    break
            if choisi is None:
                raise InertieImpossible("pivot non séparé de zéro",
                    dict(pivot=k, intervalle=list(map(str, a)), pivots_effectues=len(pivots),
                         signature_partielle=signature, operations_decimal=ar.operations))
            l, b, c, sig, det = choisi
            if positifs_seulement and sig[1]:
                raise InertieImpossible("masse non définie positive", dict(pivots=[k, l]))
            bloc = [k, l]
            voisins = sorted((set(matrice.lignes[k]) | set(matrice.lignes[l]))-set(bloc), key=rang.__getitem__)
            invdet = ar.div(ar.un, det)
            p, q, t = ar.mul(c, invdet), ar.mul(ar.neg(b), invdet), ar.mul(a, invdet)
            v = {i: (matrice.get(k, i, ar.zero), matrice.get(l, i, ar.zero)) for i in voisins}
            u = {i: (ar.add(ar.mul(p, x), ar.mul(q, y)), ar.add(ar.mul(q, x), ar.mul(t, y)))
                 for i, (x, y) in v.items()}
            preuve = dict(indices=bloc, signature=list(sig), a=list(map(str, a)),
                          b=list(map(str, b)), c=list(map(str, c)), determinant=list(map(str, det)))

            def correction(i, j):
                return ar.add(ar.mul(u[i][0], v[j][0]), ar.mul(u[i][1], v[j][1]))
        # Chaque paire est visitée avant insertion. Le budget empêche la
        # matérialisation d'une clique hors limite, même après pivot 2×2.
        largeur = max(largeur, len(voisins))
        for ii, i in enumerate(voisins):
            for j in voisins[ii:]:
                matrice.set(i, j, ar.sub(matrice.get(i, j, ar.zero), correction(i, j)))
        matrice.retirer(bloc)
        actifs.difference_update(bloc)
        signature[0] += sig[0]
        signature[1] += sig[1]
        pivots.append(preuve)
    return dict(signature=signature, pivots=pivots, largeur_max=largeur,
                pic_coefficients=matrice.pic, coefficients_crees=matrice.crees)


def _matrice_entree(x, nom):
    _verifier_canonique(x, nom)
    a = csr_matrix(x, dtype=float, copy=True)
    _verifier_canonique(a, nom+" converti")
    if not np.all(np.isfinite(a.data)):
        raise ValueError(nom+" fini requis")
    a.eliminate_zeros()
    return a


def _empreinte(a):
    h = hashlib.sha256()
    h.update(np.asarray(a.shape, dtype="<i8").tobytes())
    for x, dtype in ((a.indptr, "<i8"), (a.indices, "<i8"), (a.data, "<f8")):
        h.update(np.asarray(x, dtype=dtype).tobytes())
    return h.hexdigest()


def certifier_inertie_complement(d_i, mii, b, gamma, *, precision=32,
        budget_operations=100_000_000, budget_coefficients=2_000_000, permutation=None):
    debut = time.perf_counter()
    precision = _entier_positif(precision, "precision")
    if precision < 16:
        raise ValueError("au moins 16 chiffres Decimal requis")
    budget_operations = _entier_positif(budget_operations, "budget_operations")
    budget_coefficients = _entier_positif(budget_coefficients, "budget_coefficients")
    if np.iscomplexobj(gamma) or np.ndim(gamma) != 0:
        raise ValueError("gamma scalaire réel requis")
    gamma = float(gamma)
    if not math.isfinite(gamma) or gamma <= 0:
        raise ValueError("gamma fini strictement positif requis")
    d, m = _matrice_entree(d_i, "D_i"), _matrice_entree(mii, "M_ii")
    n = d.shape[1]
    if m.shape != (n, n) or n == 0 or np.any((m-m.T).data != 0):
        raise ValueError("M_ii exactement symétrique et compatible requise")
    if np.iscomplexobj(b):
        raise ValueError("B réel requis")
    b = np.array(b, dtype=float, copy=True)
    if b.ndim != 2 or b.shape[0] != n or not 0 < b.shape[1] < n or not np.all(np.isfinite(b)):
        raise ValueError("B fini de forme (n,s), 0<s<n requis")
    s = b.shape[1]
    if permutation is None:
        perm = np.arange(n)
    else:
        perm = np.asarray(permutation)
        if perm.dtype.kind not in "iu" or perm.shape != (n,) or not np.array_equal(np.sort(perm), np.arange(n)):
            raise ValueError("permutation physique entière complète requise")
    ar = IntervallesSignes(precision, budget_operations)
    if (m.nnz == n and np.array_equal(m.indptr, np.arange(n+1))
            and np.array_equal(m.indices, np.arange(n))):
        if np.any(m.data <= 0):
            raise InertieImpossible("masse diagonale non définie positive", dict(phase="masse"))
        # Les entrées sont des nombres binary64 interprétés exactement :
        # aucun produit, arrondi ou seuil numérique n'intervient ici.
        preuve_masse = dict(methode="diagonale_binary64_positive", signature=[n, 0, 0],
                            minimum_binary64=float(np.min(m.data)))
    else:
        masse = MatriceIntervalles(n, budget_coefficients)
        for i in range(n):
            for j, v in zip(m.indices[m.indptr[i]:m.indptr[i+1]], m.data[m.indptr[i]:m.indptr[i+1]], strict=True):
                if j >= i:
                    masse.set(i, int(j), ar.point(v))
        preuve_masse = eliminer(masse, ar, perm, positifs_seulement=True)
        preuve_masse["methode"] = "congruences_dirigees"
    if preuve_masse["signature"] != [n, 0, 0]:
        raise InertieImpossible("positivité de M_ii non établie", preuve_masse)
    apres_masse = time.perf_counter()
    a = MatriceIntervalles(n+s, budget_coefficients)
    for ligne in range(d.shape[0]):
        indices = d.indices[d.indptr[ligne]:d.indptr[ligne+1]]
        valeurs = [ar.point(v) for v in d.data[d.indptr[ligne]:d.indptr[ligne+1]]]
        for ii, i in enumerate(indices):
            for jj in range(ii, len(indices)):
                j = int(indices[jj])
                a.set(int(i), j, ar.add(a.get(int(i), j, ar.zero), ar.mul(valeurs[ii], valeurs[jj])))
    gp = ar.point(gamma)
    for i in range(n):
        for j, v in zip(m.indices[m.indptr[i]:m.indptr[i+1]], m.data[m.indptr[i]:m.indptr[i+1]], strict=True):
            if j >= i:
                a.set(i, int(j), ar.sub(a.get(i, int(j), ar.zero), ar.mul(gp, ar.point(v))))
        for j in range(s):
            a.set(i, n+j, ar.point(b[i, j]))
    apres_assemblage = time.perf_counter()
    coefficients_initiaux = a.coefficients
    preuve = eliminer(a, ar, list(map(int, perm))+list(range(n, n+s)))
    if preuve["signature"] != [n, s, 0]:
        raise InertieImpossible("complément non coercif au seuil demandé",
            dict(gamma=gamma, signature=preuve["signature"], signature_attendue=[n, s, 0],
                 preuve=preuve, operations_decimal=ar.operations))
    fin = time.perf_counter()
    return dict(lambda_min=gamma, certification_machine=True, precision_decimal=precision,
        dimension_interieure=n, rang_contraintes=s, dimension_complement=n-s,
        contraintes_sha256=hashlib.sha256(np.asarray(b, dtype="<f8").tobytes()).hexdigest(),
        d_sha256=_empreinte(d), masse_sha256=_empreinte(m),
        permutation_physique=perm.tolist(), preuve_masse=preuve_masse, preuve_kkt=preuve,
        coefficients_initiaux=coefficients_initiaux, operations_decimal=ar.operations,
        preparation_s=fin-debut, phases_s=dict(masse=apres_masse-debut,
            assemblage=apres_assemblage-apres_masse, elimination=fin-apres_assemblage),
        portee="coercivité stricte de D_i.T D_i-gamma M_ii sur ker(B.T) ; aucune réponse certifiée")
