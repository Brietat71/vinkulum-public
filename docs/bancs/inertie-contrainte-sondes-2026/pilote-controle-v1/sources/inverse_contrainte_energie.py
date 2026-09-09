"""Inverse contrainte du modèle QR par réflecteurs en coordonnées d'énergie.

K_Q=E^-1 R.T R E^-1 et les valeurs stockées de B définissent les données
cibles. Les résolutions et réflecteurs restent flottants, sans certificat
machine. Aucun inverse dense ni base dense du complément n'est formé.
"""
import numpy as np
from scipy.linalg import get_lapack_funcs, qr, solve_triangular, svdvals
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve_triangular

from vinkulum._ports.inverse_selectionnee import _verifier_canonique


def _fini_reel(a, nom):
    if np.iscomplexobj(a):
        raise ValueError(nom+" réel requis")
    a = np.asarray(a, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(nom+" fini requis")
    return a.copy()


class InverseContrainteEnergie:
    """Application de K_Q^-1 sur ker(B.T), avec stockage supplémentaire O(n r).

    Les colonnes de R^-T E B sont équilibrées puis factorisées par QR avec
    pivotage. Deux applications de réflecteurs et l'annulation des r premières
    coordonnées réalisent le projecteur orthogonal, avant la résolution R.
    Aucun calcul ne soustrait deux grandes réponses physiques. Un rang
    numériquement insuffisant est refusé ; aucune direction n'est supprimée.

    Cette classe copie ses entrées. Modifier ensuite ses propres attributs
    invalide le modèle préparé. Les contrôles de diagnostic sont évalués en
    doubles ; ils ne certifient ni le rang, ni les contraintes, ni le résidu.
    """
    def __init__(self, r, echelles, b):
        if np.iscomplexobj(r):
            raise ValueError("R réel requis")
        _verifier_canonique(r, "R")
        self.r = csr_matrix(r, dtype=float, copy=True)
        _verifier_canonique(self.r, "R converti")
        self.r.eliminate_zeros()
        self.n = self.r.shape[0]
        if (self.r.shape != (self.n, self.n) or self.n < 2
                or not np.all(np.isfinite(self.r.data))
                or np.any(self.r.diagonal() <= 0)):
            raise ValueError("R carré fini à diagonale strictement positive requis")
        for i in range(self.n):
            a, fin = self.r.indptr[i:i+2]
            if np.any(self.r.indices[a:fin] < i):
                raise ValueError("R triangulaire supérieur requis")
        self.rt = self.r.T.tocsr()
        self.echelles = _fini_reel(echelles, "E")
        self.b = _fini_reel(b, "B")
        if self.echelles.shape != (self.n,) or np.any(self.echelles <= 0):
            raise ValueError("E positif de taille n requis")
        if (self.b.ndim != 2 or self.b.shape[0] != self.n
                or not 0 < self.b.shape[1] < self.n):
            raise ValueError("B de forme (n,r), avec 0 < r < n requis")
        self.rang = self.b.shape[1]
        g = self._dual(self.b)
        self.echelles_contraintes = np.max(np.abs(g), axis=0)
        if (not np.all(np.isfinite(self.echelles_contraintes))
                or np.any(self.echelles_contraintes == 0)):
            raise ValueError("contraintes nulles ou non représentables en énergie")
        equilibre = g/self.echelles_contraintes[None, :]
        (self._reflecteurs, self._tau), self._triangle, self._permutation = qr(
            equilibre, mode="raw", pivoting=True, check_finite=False)
        valeurs = svdvals(self._triangle, check_finite=False)
        seuil = np.finfo(float).eps*max(self.n, self.rang)*valeurs[0]
        if valeurs[-1] <= seuil:
            raise ValueError("rang numérique insuffisant des contraintes en énergie")
        self.condition_triangle_equilibre = float(valeurs[0]/valeurs[-1])
        self._ormqr = get_lapack_funcs("ormqr", (self._reflecteurs,))
        self._travail = {}

    def _matrice(self, a, nom):
        a = _fini_reel(a, nom)
        if a.ndim != 2 or a.shape[0] != self.n or not a.shape[1]:
            raise ValueError(nom+" de forme (n,charges), avec au moins une colonne")
        return a

    def _dual(self, rhs):
        y = spsolve_triangular(self.rt, self.echelles[:, None]*rhs, lower=True)
        if not np.all(np.isfinite(y)):
            raise ArithmeticError("résolution duale non finie")
        return y

    def _rotation(self, a, transpose=False):
        a = np.array(a, dtype=float, order="F", copy=True)
        trans = "T" if transpose else "N"
        cle = (a.shape[1], trans)
        lwork = self._travail.get(cle)
        if lwork is None:
            _, travail, info = self._ormqr("L", trans, self._reflecteurs, self._tau,
                                           a, -1, overwrite_c=0)
            if info:
                raise ArithmeticError("dimension de travail Householder invalide")
            lwork = max(1, int(travail[0]))
            self._travail[cle] = lwork
        resultat, _, info = self._ormqr("L", trans, self._reflecteurs, self._tau,
                                        a, lwork, overwrite_c=1)
        if info or not np.all(np.isfinite(resultat)):
            raise ArithmeticError("application des réflecteurs Householder impossible")
        return resultat

    def _projeter(self, a):
        coordonnees = self._rotation(a, transpose=True)
        coordonnees[:self.rang] = 0.
        return self._rotation(coordonnees)

    def _physique(self, a):
        x = self.echelles[:, None]*spsolve_triangular(self.r, a, lower=False)
        if not np.all(np.isfinite(x)):
            raise ArithmeticError("résolution physique non finie")
        return x

    def appliquer(self, rhs):
        """Renvoie S_Q rhs, matrice (n,charges), sans différence de réponses."""
        rhs = self._matrice(rhs, "rhs")
        return self._physique(self._projeter(self._dual(rhs)))

    def dual_contraint(self, rhs):
        """Coordonnées duales (n-r,charges), avec Y.T Y≈rhs.T S_Q rhs.

        Renvoie les dernières lignes de Q.T R^-T E rhs, sans rotation retour
        ni résolution R. Leur norme euclidienne est la norme duale contrainte
        pour le modèle K_Q, à l'erreur flottante près. Pour la norme massique,
        construire cette même classe avec un facteur énergétique de M.
        """
        rhs = self._matrice(rhs, "rhs")
        return self._rotation(self._dual(rhs), transpose=True)[self.rang:].copy()

    def canonique(self):
        """Relèvement C minimal en énergie : B.T C≈I, sans former H=B.T K^-1 B.

        C vise les colonnes originales de B, avant équilibrage et pivotage.
        Le défaut de B.T C doit être contrôlé ; ce n'est pas une identité
        machine, notamment pour des colonnes d'échelles très différentes.
        """
        droite = np.zeros((self.rang, self.rang))
        droite[np.arange(self.rang), self._permutation] = 1/self.echelles_contraintes[self._permutation]
        petites = solve_triangular(self._triangle.T, droite, lower=True, check_finite=False)
        energie = np.zeros((self.n, self.rang))
        energie[:self.rang] = petites
        return self._physique(self._rotation(energie))

    def diagnostic(self, rhs, solution, d_i=None):
        """Contrôle B.T x et le résidu dual sur ker(B.T), dans Q et D éventuel.

        Le projecteur annule les réactions de contrainte : un résidu complet
        non nul peut être normal. Les nombres retournés ne sont pas des
        majorants dirigés. D_I, s'il est fourni, est utilisé par D_I.T@(D_I x),
        sans remplacer ses coefficients par un Gram flottant préassemblé.
        """
        rhs, x = self._matrice(rhs, "rhs"), self._matrice(solution, "solution")
        if rhs.shape != x.shape:
            raise ValueError("rhs et solution de mêmes dimensions requis")
        y = self._dual(rhs)
        defaut = self._projeter(y-self.r@(x/self.echelles[:, None]))
        norme = float(np.linalg.norm(y, 2))
        residu = float(np.linalg.norm(defaut, 2))
        bnorm = self.b/np.max(np.abs(self.b), axis=0)[None, :]
        resultat = dict(defaut_contrainte=float(np.linalg.norm(self.b.T@x, 2)),
            defaut_contrainte_colonnes_normalisees=float(np.linalg.norm(bnorm.T@x, 2)),
            residu_Q_dual_contraint=residu,
            residu_Q_dual_contraint_relatif=residu/norme if norme else (0. if residu == 0 else None),
            norme_charge_duale=norme, condition_triangle_equilibre=self.condition_triangle_equilibre,
            certification_machine=False)
        if d_i is not None:
            if np.iscomplexobj(d_i):
                raise ValueError("D_I réel requis")
            d = csr_matrix(d_i, dtype=float, copy=True)
            if d.shape[1] != self.n or not np.all(np.isfinite(d.data)):
                raise ValueError("D_I fini avec n colonnes requis")
            residu_d = float(np.linalg.norm(self._projeter(self._dual(rhs-d.T@(d@x))), 2))
            resultat.update(residu_D_dual_contraint=residu_d,
                residu_D_dual_contraint_relatif=residu_d/norme if norme else (0. if residu_d == 0 else None))
        return resultat
