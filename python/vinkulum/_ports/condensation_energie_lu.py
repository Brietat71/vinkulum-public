"""Préconditionnement LU et raffinement dans le facteur d'énergie original.

Le Gram creux K_a=fl(D_I.T D_I), équilibré diagonalement, sert seulement à
construire une LU. Les produits, Gram, résidus de raffinement et énergies
statiques sont évalués par D_I. Une résolution conserve au plus deux
corrections ; le relèvement statique en conserve au plus trois.

Cette voie doit être choisie explicitement. Elle ne remplace aucun refus QR
et n'offre aucun certificat arithmétique. ``dual(rhs)=D_I solve(rhs)`` n'est
qu'une représentation approchée de la norme K^-1 : l'audit indépendant de
PortsReleves dans l'énergie d'entrée et la masse reste nécessaire.
"""
import operator
import time

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import splu


class BudgetLU(RuntimeError):
    pass


class CondensationEnergieLU:
    """Interface de CondensationEnergie, avec préconditionneur LU explicite.

    ``budget`` limite le stockage L+U constaté après factorisation, pas son
    pic mémoire. Une LU de rang insuffisant, des pivots non positifs ou une
    estimation déjà alarmante de kappa*epsilon provoquent un refus. Les
    contrôles restants ne prouvent pas la précision de l'inverse approché.
    """

    def __init__(self, d, interieur, interface, metrique, budget=10_000_000):
        debut = time.perf_counter()
        if isinstance(budget, (bool, np.bool_)):
            raise ValueError("budget LU entier positif requis")
        try:
            budget = operator.index(budget)
        except TypeError as exc:
            raise ValueError("budget LU entier positif requis") from exc
        if budget < 1:
            raise ValueError("budget LU entier positif requis")
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
        self.i, self.s = i.copy(), s.copy()
        self.di = self.d[:, self.i].tocsc()
        self._ds = self.d[:, self.s]@self.w
        self._k = (self.di.T@self.di).tocsc()
        diag = self._k.diagonal()
        if np.any(diag <= 0) or not np.all(np.isfinite(self._k.data)):
            raise ValueError("Gram intérieur nul ou hors domaine flottant")
        self.echelles = 1/np.sqrt(diag)
        e = diags(self.echelles, format="csc")
        self._ke = (e@self._k@e).tocsc()
        try:
            self._lu = splu(self._ke, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.,
                            options={"SymmetricMode": True})
        except RuntimeError as exc:
            raise ValueError("Gram assemblé inutilisable comme préconditionneur LU") from exc
        pivots = self._lu.U.diagonal()
        if (not np.array_equal(self._lu.perm_r, self._lu.perm_c)
                or not np.all(np.isfinite(pivots)) or np.any(pivots <= 0)):
            raise ValueError("préconditionneur LU sans pivots SPD fiables")
        self.nnz_facteurs = self._lu.L.nnz+self._lu.U.nnz
        if self.nnz_facteurs > budget:
            raise BudgetLU("budget de coefficients du préconditionneur LU dépassé")
        self.budget = budget
        self.compteurs = {
            "factorisations_LU": 1, "resolutions_LU": 0, "colonnes_rhs_LU": 0,
            "resolutions_conditionnement": 0, "colonnes_rhs_conditionnement": 0,
            "appels_solve": 0, "appels_dual": 0,
            "corrections_solve_tentees": 0, "corrections_statiques_tentees": 0,
            "colonnes_corrections_solve_retenues": 0, "colonnes_corrections_statiques_retenues": 0,
        }
        self.condition_estimee = self._condition()
        self.kappa_epsilon_estime = self.condition_estimee*np.finfo(float).eps
        if not np.isfinite(self.kappa_epsilon_estime) or self.kappa_epsilon_estime >= .1:
            raise ArithmeticError("estimation kappa*epsilon du Gram LU supérieure ou égale à 0.1")
        self.statut_lu = "preconditionneur_admis_sous_audit"
        self.residu_resolution_relatif_max = 0.
        rhs = -(self.di.T@self._ds)
        psi = self._resout_lu(rhs)
        h = self.di@psi+self._ds
        residu = self.di.T@h
        scores = self._scores(residu)
        self.historique_statique = [float(np.max(scores, initial=0.))]
        for _ in range(3):
            correction = self._resout_lu(residu)
            self.compteurs["corrections_statiques_tentees"] += 1
            candidat = psi-correction
            h_candidat = self.di@candidat+self._ds
            r_candidat = self.di.T@h_candidat
            suivants = self._scores(r_candidat)
            garder = suivants < scores
            self.compteurs["colonnes_corrections_statiques_retenues"] += int(np.count_nonzero(garder))
            psi[:, garder], h[:, garder], residu[:, garder] = (
                candidat[:, garder], h_candidat[:, garder], r_candidat[:, garder])
            scores[garder] = suivants[garder]
            self.historique_statique.append(float(np.max(scores, initial=0.)))
        self.psi, self.h, self.r0 = psi, h, residu
        self.k0 = self.h.T@self.h
        self.coefficients_ports = self.psi.size+self._ds.size
        self.preparation_s = time.perf_counter()-debut

    def _condition(self):
        """Estimation déterministe de kappa_1 du Gram équilibré, sans borne haute.

        Chaque sonde x de norme 1 donne ||K_e^-1 x||_1. Des directions
        transposées choisissent les sondes suivantes, au plus cinq fois.
        Une dernière sonde alternée complète l'estimation. Une petite valeur
        estimée ne prouve donc jamais un bon conditionnement.
        """
        n = self._ke.shape[0]
        x = np.full(n, 1/n)
        estimation, precedent = 0., -1

        def inverse(v, transpose=False):
            self.compteurs["resolutions_LU"] += 1
            self.compteurs["colonnes_rhs_LU"] += 1
            self.compteurs["resolutions_conditionnement"] += 1
            self.compteurs["colonnes_rhs_conditionnement"] += 1
            result = self._lu.solve(v, trans="T" if transpose else "N")
            if not np.all(np.isfinite(result)):
                raise ArithmeticError("sonde de conditionnement LU hors domaine flottant")
            return result

        for _ in range(5):
            y = inverse(x)
            estimation = max(estimation, float(np.linalg.norm(y, 1)))
            signe = np.where(y >= 0, 1., -1.)
            z = inverse(signe, True)
            j = int(np.argmax(np.abs(z)))
            if j == precedent or abs(z[j]) <= float(np.dot(z, x)):
                break
            x = np.zeros(n)
            x[j], precedent = 1., j
        x = (-1.)**np.arange(n)*(1+np.arange(n)/max(n-1, 1))
        x /= np.linalg.norm(x, 1)
        estimation = max(estimation, float(np.linalg.norm(inverse(x), 1)))
        norme = float(np.max(np.asarray(abs(self._ke).sum(axis=0))))
        return norme*estimation

    def _tableau(self, x, nom):
        x = np.asarray(x, dtype=float)
        if x.ndim not in (1, 2) or x.shape[0] != len(self.i) or not np.all(np.isfinite(x)):
            raise ValueError(nom+" doit être fini et avoir la dimension intérieure")
        return (x[:, None], True) if x.ndim == 1 else (x, False)

    def _scores(self, residu):
        scores = np.max(np.abs(self.echelles[:, None]*residu), axis=0)
        if not np.all(np.isfinite(scores)):
            raise ArithmeticError("résidu énergétique LU hors domaine flottant")
        return scores

    def _resout_lu(self, rhs):
        self.compteurs["resolutions_LU"] += 1
        self.compteurs["colonnes_rhs_LU"] += rhs.shape[1]
        result = self.echelles[:, None]*self._lu.solve(self.echelles[:, None]*rhs)
        if not np.all(np.isfinite(result)):
            raise ArithmeticError("résolution LU hors domaine flottant")
        return result

    def coordonnees(self, x):
        x, vectoriel = self._tableau(x, "x")
        result = self.di@x
        return result[:, 0] if vectoriel else result

    def produit(self, x):
        x, vectoriel = self._tableau(x, "x")
        result = self.di.T@(self.di@x)
        return result[:, 0] if vectoriel else result

    def gram(self, x, y):
        return self.coordonnees(x).T@self.coordonnees(y)

    def solve(self, rhs):
        rhs, vectoriel = self._tableau(rhs, "rhs")
        self.compteurs["appels_solve"] += 1
        if rhs.shape[1] == 0 or not np.any(rhs):
            result = np.zeros_like(rhs)
            return result[:, 0] if vectoriel else result
        x = self._resout_lu(rhs)
        residu = rhs-self.di.T@(self.di@x)
        scores = self._scores(residu)
        for _ in range(2):
            correction = self._resout_lu(residu)
            self.compteurs["corrections_solve_tentees"] += 1
            candidat = x+correction
            r_candidat = rhs-self.di.T@(self.di@candidat)
            suivants = self._scores(r_candidat)
            garder = suivants < scores
            self.compteurs["colonnes_corrections_solve_retenues"] += int(np.count_nonzero(garder))
            x[:, garder], residu[:, garder] = candidat[:, garder], r_candidat[:, garder]
            scores[garder] = suivants[garder]
        denominateurs = np.maximum(self._scores(rhs), np.finfo(float).tiny)
        self.residu_resolution_relatif_max = max(self.residu_resolution_relatif_max,
                                                 float(np.max(scores/denominateurs, initial=0.)))
        return x[:, 0] if vectoriel else x

    def dual(self, rhs):
        """Coordonnées duales approchées ; leur norme n'est pas certifiée."""
        self.compteurs["appels_dual"] += 1
        return self.di@self.solve(rhs)

    def matrice(self):
        """Gram double du préconditionneur ; les produits utilisent D original."""
        return self._k.copy()

    def audit(self):
        """Diagnostics dans D original, sans encadrement des erreurs d'arrondi."""
        x = np.zeros((self.d.shape[1], len(self.s)))
        x[self.i], x[self.s] = self.psi, self.w
        dx = self.d@x
        defaut = self.di.T@dx-self.r0
        # Cette évaluation de norme duale reste approchée, comme tout dual().
        defaut_dual = float(np.linalg.norm(self.dual(defaut), 2))
        lifting_dual = float(np.linalg.norm(self.dual(self.r0), 2))
        return dict(backend="lu_energie", statut_lu=self.statut_lu,
                    ecart_energie_statique=float(np.linalg.norm(dx.T@dx-self.k0, 2)),
                    defaut_equilibre_entree_dual=defaut_dual,
                    residu_lifting=lifting_dual,
                    norme_deformation_lifting=float(np.linalg.norm(self.h, 2)),
                    residu_statique_equilibre=float(np.max(self._scores(self.r0), initial=0.)),
                    condition_estimee_equilibree=self.condition_estimee,
                    kappa_epsilon_estime=self.kappa_epsilon_estime,
                    residu_resolution_relatif_max=self.residu_resolution_relatif_max,
                    nnz_L=self._lu.L.nnz, nnz_U=self._lu.U.nnz,
                    coefficients_ports=self.coefficients_ports,
                    compteurs=dict(self.compteurs), certification_machine=False,
                    norme_duale="approchee_par_D_I_solve")
