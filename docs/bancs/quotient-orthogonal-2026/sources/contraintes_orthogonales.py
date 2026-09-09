"""Expérience de quotient orthogonal local, indépendante de l'API livrée.

Élimination à droite par réflexions de Householder. Aucun G Gᵀ ni base
dense du noyau. L'ordre glouton minimise la taille de la ligne active ;
ce n'est ni un calcul de largeur arborescente ni un QR révélateur de rang.
Voir docs/CONTRAINTES_ORTHOGONALES_PROTOTYPE.md pour le contrat numérique.
"""
import heapq
import math

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve_triangular


class RangAmbigu(RuntimeError):
    """Le résidu de ligne est trop proche du seuil de suppression."""


class BudgetDepasse(RuntimeError):
    """Le travail scalaire demandé dépasse le budget explicite."""


def normalise(g):
    if np.iscomplexobj(g):
        raise ValueError('contraintes réelles requises')
    raw = csr_matrix(g, dtype=float, copy=True)
    raw.sum_duplicates()
    raw.eliminate_zeros()
    if not np.isfinite(raw.data).all():
        raise ValueError('contraintes non finies')
    # Deux divisions évitent le débordement de la norme d'une ligne.
    scales, norms = np.ones(raw.shape[0]), np.ones(raw.shape[0])
    for i in range(raw.shape[0]):
        row = raw.data[raw.indptr[i]:raw.indptr[i+1]]
        if row.size:
            scales[i] = max(abs(row))
            row /= scales[i]
            norms[i] = np.linalg.norm(row)
            row /= norms[i]
    return raw, scales, norms


class QuotientOrthogonal:
    def __init__(self, g, *, seuil=None, budget=10_000_000):
        raw, self.echelles, self.normes = normalise(g)
        self.m, self.n = raw.shape
        self.seuil = float(seuil) if seuil is not None else 64*np.finfo(float).eps*max(1, *raw.shape)
        if not np.isfinite(self.seuil) or not 0 < self.seuil < 1:
            raise ValueError('seuil fini entre 0 et 1 requis')
        if not isinstance(budget, int) or isinstance(budget, bool) or budget < 0:
            raise ValueError('budget entier positif ou nul requis')
        self.g = raw
        rows = [dict(zip(raw.indices[raw.indptr[i]:raw.indptr[i+1]].tolist(),
                         raw.data[raw.indptr[i]:raw.indptr[i+1]].tolist())) for i in range(self.m)]
        adjacent = [set() for _ in range(self.n)]
        for i, row in enumerate(rows):
            for j in row:
                adjacent[j].add(i)
        heap = [(len(row), i) for i, row in enumerate(rows)]
        heapq.heapify(heap)
        coefficients = [{} for _ in rows]
        self.reflexions = []
        pivots, selected = [], []
        dropped = np.zeros(self.m)
        pivot_values = []
        work = 0
        live_nnz = raw.nnz
        peak_nnz = live_nnz
        largest_support = 0
        while heap:
            size, i = heapq.heappop(heap)
            row = rows[i]
            if row is None or size != len(row):
                continue
            norm = math.hypot(*row.values())
            if self.seuil/4 <= norm <= 4*self.seuil:
                raise RangAmbigu(f'ligne {i} : résidu {norm:.6g}, seuil {self.seuil:.6g}')
            for j in row:
                adjacent[j].remove(i)
            rows[i] = None
            live_nnz -= len(row)
            if norm < self.seuil:
                dropped[i] = norm
                continue
            # Le signe sûr évite la soustraction de deux nombres proches.
            ids = np.array(sorted(row), dtype=np.int64)
            p = min(row, key=lambda j: (len(adjacent[j]), j))
            loc = int(np.searchsorted(ids, p))
            u = np.array([row[j] for j in ids])
            alpha = -math.copysign(norm, u[loc])
            u[loc] -= alpha
            u /= np.linalg.norm(u)
            affected = set().union(*(adjacent[j] for j in ids))
            # Inclut la réflexion de la ligne choisie. Unité : coefficient
            # d'un panneau à transformer, pas FLOP ni octet alloué.
            work += len(ids)*(len(affected)+1)
            if work > budget:
                raise BudgetDepasse(f'{work} coefficients à transformer > budget {budget}')
            largest_support = max(largest_support, len(ids))
            self.reflexions.append((ids, u))
            order = len(pivots)
            coefficients[i][order] = alpha
            pivots.append(p)
            selected.append(i)
            pivot_values.append(norm)
            for other in sorted(affected):
                target = rows[other]
                x = np.array([target.get(j, 0.) for j in ids])
                x -= 2*np.dot(x, u)*u
                if not np.isfinite(x).all():
                    raise RuntimeError('élimination non finie')
                for j, value in zip(ids, x):
                    if j in target:
                        del target[j]
                        adjacent[j].remove(other)
                        live_nnz -= 1
                    if j == p:
                        if value != 0.:
                            coefficients[other][order] = float(value)
                    elif value != 0.:
                        # Pas de suppression de petits coefficients : le
                        # seul seuil intervient sur une ligne résiduelle.
                        target[j] = float(value)
                        adjacent[j].add(other)
                        live_nnz += 1
                heapq.heappush(heap, (len(target), other))
            peak_nnz = max(peak_nnz, live_nnz)
        self.pivots = np.array(pivots, dtype=np.int64)
        self.indices = np.array(selected, dtype=np.int64)
        self.rang = len(pivots)
        self.libres = np.setdiff1d(np.arange(self.n), self.pivots)
        ri, ci, values = [], [], []
        for i, row in enumerate(coefficients):
            for j, value in row.items():
                ri.append(i)
                ci.append(j)
                values.append(value)
        self.r = csr_matrix((values, (ri, ci)), shape=(self.m, self.rang))
        self.t = self.r[self.indices].tocsr()
        self.diagnostic = dict(
            lignes=self.m, colonnes=self.n, rang=self.rang, seuil=self.seuil,
            nnz_entree=int(raw.nnz), nnz_actifs_max=int(peak_nnz),
            nnz_facteur=int(self.r.nnz), coefficients_reflexions=sum(len(u) for _, u in self.reflexions),
            support_max=largest_support, coefficients_transformes=work,
            lignes_supprimees=int(self.m-self.rang),
            norme_residus_supprimes=float(np.linalg.norm(dropped)),
            pivot_min=min(pivot_values, default=None),
            orthogonalite_locale_max=max((abs(float(u@u)-1.) for _, u in self.reflexions), default=0.),
            base_dense=False, gram_forme=False,
        )

    def orthogonal(self, v, *, transpose=False):
        """Q v, ou Qᵀ v ; Q=H₁…Hᵣ. Accepte un vecteur ou plusieurs colonnes."""
        out = np.array(v, dtype=float, copy=True)
        if out.ndim not in (1, 2) or out.shape[0] != self.n or not np.isfinite(out).all():
            raise ValueError('vecteur ou matrice fini de hauteur n requis')
        reflectors = self.reflexions if transpose else reversed(self.reflexions)
        for ids, u in reflectors:
            block = out[ids]
            dot = u@block
            out[ids] = block - 2*(u*dot if out.ndim == 1 else u[:, None]*dot)
        if not np.isfinite(out).all():
            raise RuntimeError('application orthogonale non finie')
        return out

    def projette(self, v):
        out = self.orthogonal(v, transpose=True)
        out[self.pivots] = 0.
        return self.orthogonal(out)

    def injecte(self, y):
        y = np.asarray(y, dtype=float)
        if y.ndim not in (1, 2) or y.shape[0] != len(self.libres):
            raise ValueError('hauteur égale à la dimension libre requise')
        out = np.zeros((self.n,)+y.shape[1:])
        out[self.libres] = y
        return self.orthogonal(out)

    def restreint(self, v):
        return self.orthogonal(v, transpose=True)[self.libres]

    def reactions_normalisees(self, force):
        """G_normaliséᵀ η = (I-P) force, sans minimisation de la norme de η."""
        rhs = self.orthogonal(force, transpose=True)[self.pivots]
        out = np.zeros((self.m,)+rhs.shape[1:])
        if self.rang:
            out[self.indices] = spsolve_triangular(self.t.T.tocsc(), rhs, lower=False)
        return out

    def reactions(self, force):
        eta = self.reactions_normalisees(force)
        scales = self.echelles if eta.ndim == 1 else self.echelles[:, None]
        norms = self.normes if eta.ndim == 1 else self.normes[:, None]
        result = (eta/norms)/scales
        if not np.isfinite(result).all():
            raise RuntimeError('réactions originales non finies')
        return result
