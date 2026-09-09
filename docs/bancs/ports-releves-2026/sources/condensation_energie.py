"""Condensation statique par rotations de Givens sur le facteur d'énergie.

D décrit l'énergie 1/2 ||D u||². Les colonnes intérieures sont équilibrées
avant QR ; les ports sont normalisés par W. L'ordre des colonnes est fixé,
sans pivotage. Le remplissage et un budget explicite limitent le prototype.
Les facteurs calculés définissent un modèle QR voisin ; leur écart au
modèle énergétique d'entrée doit être audité, sans certificat machine.
"""
import time

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.sparse import csc_matrix, csr_matrix, diags
from scipy.sparse.linalg import spsolve_triangular


class BudgetQR(RuntimeError):
    pass


class CondensationEnergie:
    def __init__(self, d, interieur, interface, metrique, budget=10_000_000):
        debut = time.perf_counter()
        self.d = csr_matrix(d, dtype=float, copy=True)
        self.d.sum_duplicates()
        self.d.eliminate_zeros()
        n = self.d.shape[1]
        i, s = np.asarray(interieur), np.asarray(interface)
        if (i.ndim != 1 or s.ndim != 1 or i.dtype.kind not in "iu" or s.dtype.kind not in "iu"
                or len(i) == 0 or len(s) == 0 or not np.all(np.isfinite(self.d.data))
                or not np.array_equal(np.sort(np.r_[i, s]), np.arange(n))):
            raise ValueError("D fini et partition physique complète requis")
        metric = np.asarray(metrique, dtype=float)
        if (metric.shape != (len(s), len(s)) or not np.all(np.isfinite(metric))
                or not np.allclose(metric, metric.T, rtol=0., atol=1e-14*np.max(abs(metric)))):
            raise ValueError("métrique symétrique finie requise")
        l = cholesky(metric, lower=True)
        self.w = solve_triangular(l.T, np.eye(len(s)), lower=False)
        self.i, self.s = i, s
        self.di = self.d[:, i].tocsc()
        ds = self.d[:, s]@self.w
        normes = np.sqrt(np.asarray(self.di.power(2).sum(axis=0)).ravel())
        if np.any(normes == 0) or not np.all(np.isfinite(normes)):
            raise ValueError("colonne intérieure nulle ou hors domaine")
        self.echelles = 1/normes
        a = (self.di@diags(self.echelles)).tocsr()
        lignes = [dict() for _ in i]
        port = np.zeros((len(i), len(s)))
        restes = []
        self.rotations = self.operations = 0
        self.budget = int(budget)
        for row in range(a.shape[0]):
            u, v = a.indptr[row:row+2]
            courant = dict(zip(a.indices[u:v].tolist(), a.data[u:v].tolist()))
            bord = ds[row].copy()
            while courant:
                j = min(courant)
                b = courant.pop(j)
                if b == 0.:
                    continue
                ancien = lignes[j]
                aa = ancien.get(j, 0.)
                h = np.hypot(aa, b)
                if not np.isfinite(h):
                    raise ArithmeticError("rotation QR hors domaine flottant")
                c, sn = aa/h, b/h
                colonnes = (ancien.keys() | courant.keys())-{j}
                self.operations += 6*(len(colonnes)+len(s))
                self.rotations += 1
                if self.operations > self.budget:
                    raise BudgetQR("budget de mises à jour Givens dépassé")
                suivant, nouveau = {}, {j: h}
                for k in colonnes:
                    x, y = ancien.get(k, 0.), courant.get(k, 0.)
                    val, reste = c*x+sn*y, -sn*x+c*y
                    if val != 0.:
                        nouveau[k] = val
                    if reste != 0.:
                        suivant[k] = reste
                ancien_port = port[j].copy()
                port[j] = c*ancien_port+sn*bord
                bord = -sn*ancien_port+c*bord
                lignes[j] = nouveau
                courant = suivant
            if np.any(bord):
                restes.append(bord)
        if any(not row or row.get(j, 0.) <= 0 for j, row in enumerate(lignes)):
            raise ValueError("intérieur de rang insuffisant pour le QR sans pivotage")
        rows, cols, vals = [], [], []
        for j, row in enumerate(lignes):
            for k, val in sorted(row.items()):
                rows.append(j); cols.append(k); vals.append(val)
        self.r = csr_matrix((vals, (rows, cols)), shape=(len(i), len(i)))
        self.rt = self.r.T.tocsr()
        self.p = port
        self.e = np.array(restes).reshape(-1, len(s))
        self.nnz_facteurs = self.r.nnz
        self.coefficients_ports = self.p.size+self.e.size
        self.psi = -self.echelles[:, None]*spsolve_triangular(self.r, self.p, lower=False)
        # Défaut du relèvement stocké, conservé dans K0 et le couplage R0.
        self.h = self.r@(self.psi/self.echelles[:, None])+self.p
        self.k0 = self.e.T@self.e+self.h.T@self.h
        self.r0 = (self.rt@self.h)/self.echelles[:, None]
        self.preparation_s = time.perf_counter()-debut

    def coordonnees(self, x):
        return self.r@(x/self.echelles[:, None])

    def produit(self, x):
        return (self.rt@self.coordonnees(x))/self.echelles[:, None]

    def gram(self, x, y):
        return self.coordonnees(x).T@self.coordonnees(y)

    def dual(self, rhs):
        return spsolve_triangular(self.rt, self.echelles[:, None]*rhs, lower=True)

    def solve(self, rhs):
        return self.echelles[:, None]*spsolve_triangular(self.r, self.dual(rhs), lower=False)

    def matrice(self):
        inv = diags(1/self.echelles)
        return (inv@(self.rt@self.r)@inv).tocsc()

    def audit(self):
        """Résidus de modèle sur les ports et Gram ; diagnostics, pas certificats."""
        x = np.zeros((self.d.shape[1], len(self.s)))
        x[self.i], x[self.s] = self.psi, self.w
        dx = self.d@x
        statique_entree = dx.T@dx
        defaut = self.di.T@dx-self.r0
        return dict(ecart_energie_statique=float(np.linalg.norm(statique_entree-self.k0, 2)),
                    defaut_equilibre_entree_dual=float(np.linalg.norm(self.dual(defaut), 2)),
                    residu_lifting=float(np.linalg.norm(self.h, 2)),
                    rotations=self.rotations, operations=self.operations,
                    nnz_R=self.r.nnz, coefficients_ports=self.coefficients_ports,
                    certification_machine=False)
