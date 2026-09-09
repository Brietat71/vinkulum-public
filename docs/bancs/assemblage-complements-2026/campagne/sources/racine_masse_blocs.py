"""Racine massique creuse par petits blocs physiques, sans certificat.

Les données stockées doivent être exactement symétriques. La positivité
est un contrôle numérique par Cholesky, pas une preuve d'arrondi dirigé.
Aucune conversion dense de la matrice globale ni permutation physique.
"""
import operator

import numpy as np
from scipy.linalg import cholesky
from scipy.sparse import csr_matrix, issparse
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve_triangular


def _budget(valeur):
    if isinstance(valeur, (bool, np.bool_)):
        raise ValueError("taille_bloc_max doit être un entier de 1 à 6")
    try:
        valeur = operator.index(valeur)
    except TypeError as exc:
        raise ValueError("taille_bloc_max doit être un entier de 1 à 6") from exc
    if not 1 <= valeur <= 6:
        raise ValueError("taille_bloc_max doit être un entier de 1 à 6")
    return valeur


def _valider_stockage(m):
    """Ne pas se fier aux drapeaux canoniques éventuellement mis en cache."""
    if not issparse(m) or m.format not in ("csr", "csc"):
        raise ValueError("M au format CSR ou CSC requis")
    if m.dtype.kind != "f" or m.dtype.itemsize != 8:
        raise ValueError("M doit contenir des valeurs réelles binary64")
    if len(m.shape) != 2 or not m.shape[0] or m.shape[0] != m.shape[1]:
        raise ValueError("M carrée non vide requise")
    ptr, ids, donnees = m.indptr, m.indices, m.data
    n = m.shape[0]
    if (ptr.ndim != 1 or ids.ndim != 1 or donnees.ndim != 1
            or ptr.dtype.kind not in "iu" or ids.dtype.kind not in "iu"
            or len(ptr) != n+1 or not len(ptr) or ptr[0] != 0
            or ptr[-1] != len(ids) or len(donnees) != len(ids)
            or np.any(ptr < 0) or np.any(ptr > len(ids))
            or np.any(ptr[1:] < ptr[:-1])
            or np.any(ids < 0) or np.any(ids >= n)):
        raise ValueError("stockage CSR/CSC invalide")
    if not np.all(np.isfinite(donnees)):
        raise ValueError("M finie requise")
    if len(ids) > 1:
        # Exclure les frontières des lignes/colonnes, y compris vides.
        mauvais = ids[1:] <= ids[:-1]
        frontieres = ptr[1:-1]
        frontieres = frontieres[(frontieres > 0) & (frontieres < len(ids))]
        mauvais[frontieres-1] = False
        if np.any(mauvais):
            raise ValueError("M canonique requise : indices triés sans doublons")


class RacineMasseBlocs:
    """Construit un R CSR supérieur tel que M ≈ R.T R numériquement.

    Les blocs sont les composantes connexes du motif non nul de M ;
    leurs indices peuvent être dispersés dans les coordonnées physiques.
    Chaque bloc comporte au plus taille_bloc_max coordonnées (1 à 6).
    Le cas diagonal ne lance ni recherche de composantes ni Cholesky.

    `dual(rhs)` calcule R^-T rhs. `norme(champ)` renvoie ||R champ||_2
    pour un vecteur et les normes par colonne pour une matrice. Les RHS
    sont denses, réels et finis, de forme (n,) ou (n,k), k>0.

    M est copiée. Le facteur public `r` sert à la lecture ; sa mutation
    après construction n'est pas supportée. Aucune certification machine.
    """
    certification_machine = False

    def __init__(self, m, taille_bloc_max=6):
        max_bloc = _budget(taille_bloc_max)
        _valider_stockage(m)
        masse = csr_matrix(m, dtype=np.float64, copy=True)
        masse.eliminate_zeros()
        self.n = masse.shape[0]
        self._diagonale = None
        if (masse.nnz == self.n
                and np.array_equal(masse.indptr, np.arange(self.n+1))
                and np.array_equal(masse.indices, np.arange(self.n))):
            if np.any(masse.data <= 0):
                raise ValueError("M définie positive requise : diagonale non positive")
            self._diagonale = np.sqrt(masse.data)
            self.r = csr_matrix((self._diagonale.copy(), masse.indices.copy(),
                                 masse.indptr.copy()), shape=masse.shape)
            return

        transpose = masse.T.tocsr()
        if (not np.array_equal(masse.indptr, transpose.indptr)
                or not np.array_equal(masse.indices, transpose.indices)
                or not np.array_equal(masse.data, transpose.data)):
            raise ValueError("M exactement symétrique requise")
        nombre, labels = connected_components(masse, directed=False)
        tailles = np.bincount(labels, minlength=nombre)
        if np.any(tailles > max_bloc):
            raise ValueError("taille de bloc de masse supérieure à taille_bloc_max")
        groupes = [[] for _ in range(nombre)]
        for indice, label in enumerate(labels):
            groupes[label].append(indice)  # Déjà triés dans l'ordre physique.
        lignes, colonnes, valeurs = [], [], []
        for groupe in groupes:
            ids = np.asarray(groupe, dtype=np.intp)
            bloc = masse[ids][:, ids].toarray()
            try:
                facteur = cholesky(bloc, lower=False, check_finite=False)
            except np.linalg.LinAlgError as exc:
                raise ValueError("bloc de M non défini positif numériquement") from exc
            if (not np.all(np.isfinite(facteur))
                    or np.any(np.diag(facteur) <= 0)):
                raise ValueError("facteur de masse non fini ou non positif")
            i, j = np.triu_indices(len(ids))
            lignes.extend(ids[i])
            colonnes.extend(ids[j])
            valeurs.extend(facteur[i, j])
        self.r = csr_matrix((valeurs, (lignes, colonnes)), shape=masse.shape,
                            dtype=np.float64)
        self.r.eliminate_zeros()
        self._transpose = self.r.T.tocsr()

    def _rhs(self, rhs):
        if np.iscomplexobj(rhs):
            raise ValueError("second membre réel requis")
        try:
            rhs = np.asarray(rhs, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError("second membre dense réel requis") from exc
        if (rhs.ndim not in (1, 2) or rhs.shape[0] != self.n
                or (rhs.ndim == 2 and rhs.shape[1] == 0)
                or not np.all(np.isfinite(rhs))):
            raise ValueError("second membre fini de forme (n,) ou (n,k), k>0 requis")
        return rhs

    def dual(self, rhs):
        """Applique R^-T aux RHS, sans former d'inverse dense."""
        rhs = self._rhs(rhs)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            if self._diagonale is not None:
                denom = self._diagonale if rhs.ndim == 1 else self._diagonale[:, None]
                resultat = rhs/denom
            else:
                resultat = spsolve_triangular(self._transpose, rhs, lower=True)
        if not np.all(np.isfinite(resultat)):
            raise ArithmeticError("résolution massique non finie")
        return resultat

    def norme(self, champ):
        """Normes massiques par colonne ; calcul flottant non certifié."""
        champ = self._rhs(champ)
        with np.errstate(over="ignore", invalid="ignore"):
            transforme = self.r@champ
            # hypot évite le débordement artificiel de la somme des carrés.
            resultat = np.hypot.reduce(transforme, axis=0, initial=0.0)
        if not np.all(np.isfinite(resultat)):
            raise ArithmeticError("norme massique non finie")
        return float(resultat) if champ.ndim == 1 else resultat
