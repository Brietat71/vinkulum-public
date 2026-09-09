"""Trace massique par inverse sélectionnée classique (Takahashi).

Pour A = E^-1 R.T R E^-1, calcule les seules entrées de (R.T R)^-1
nécessaires à trace(M A^-1), plus leur fermeture de dépendances.
La minoration 1/trace est exacte pour les données exactes sous SPD.
Les calculs flottants et l'écart entre le facteur QR et D d'entrée ne
sont pas certifiés : aucun certificat du modèle D original n'est émis.
"""
import math
import operator

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.csgraph import connected_components


class BudgetInverseSelectionnee(RuntimeError):
    """Refus explicite avant de dépasser le budget de fermeture ou de calcul."""

    def __init__(self, message, bilan):
        super().__init__(message)
        self.bilan = dict(bilan)


def _entier_positif(valeur, nom):
    if isinstance(valeur, (bool, np.bool_)):
        raise ValueError(nom+" doit être un entier strictement positif")
    try:
        entier = operator.index(valeur)
    except TypeError as exc:
        raise ValueError(nom+" doit être un entier strictement positif") from exc
    if entier <= 0:
        raise ValueError(nom+" doit être un entier strictement positif")
    return entier


def _verifier_canonique(a, nom):
    """Contrôle les indices réels ; un drapeau sparse mis en cache ne suffit pas."""
    if np.iscomplexobj(a):
        raise ValueError(nom+" doit être réel")
    format_sparse = getattr(a, "format", None)
    if format_sparse is not None and format_sparse not in ("csr", "csc", "coo", "lil", "dok"):
        raise ValueError(nom+" : convertir explicitement en CSR/CSC avant certification")
    if format_sparse in ("csr", "csc"):
        ptr, ids = a.indptr, a.indices
        majeur, mineur = a.shape if format_sparse == "csr" else a.shape[::-1]
        if (len(ptr) != majeur+1 or ptr[0] != 0 or ptr[-1] != len(ids)
                or len(a.data) != len(ids) or np.any(ptr[1:] < ptr[:-1])
                or np.any(ids < 0) or np.any(ids >= mineur)):
            raise ValueError(nom+" : stockage sparse invalide")
        for i in range(majeur):
            row = ids[ptr[i]:ptr[i+1]]
            if np.any(row[1:] <= row[:-1]):
                raise ValueError(nom+" doit être canonique : doublons ou indices non triés")
    elif format_sparse == "coo":
        row, col = a.row, a.col
        if (np.any(row < 0) or np.any(row >= a.shape[0])
                or np.any(col < 0) or np.any(col >= a.shape[1])
                or np.any((row[1:] < row[:-1])
                          | ((row[1:] == row[:-1]) & (col[1:] <= col[:-1])))):
            raise ValueError(nom+" COO doit être ordonné sans doublons")
    elif hasattr(a, "has_canonical_format") and not a.has_canonical_format:
        raise ValueError(nom+" doit être canonique : pas de somme implicite de doublons")


class InverseSelectionnee:
    """Inverse partielle d'un facteur supérieur, avec trace massique.

    `masse` est soit un vecteur diagonal positif, soit une matrice SPD
    par blocs connexes de taille au plus `taille_bloc_masse_max` (6 par
    défaut). Les blocs peuvent être dispersés par une permutation ; ils
    ne doivent pas être contigus. Un couplage hors motif de R entraîne
    une fermeture supplémentaire, jamais son omission silencieuse.

    Les budgets portent sur le nombre de coefficients sélectionnés et
    sur les visites de dépendances plus les opérations scalaires du
    calcul. Ils ne représentent pas des octets ou des FLOP mesurés.
    La validation d'entrée coûte en outre O(nnz(R)+nnz(M)), et O(b^3)
    par petit bloc de masse de taille b. Aucune inverse dense, aucune
    résolution par colonne et aucun recours automatique à celles-ci.
    """

    def __init__(self, r, masse, echelles=None, *,
                 budget_coefficients=2_000_000,
                 budget_operations=50_000_000,
                 taille_bloc_masse_max=6):
        self.budget_coefficients = _entier_positif(budget_coefficients, "budget_coefficients")
        self.budget_operations = _entier_positif(budget_operations, "budget_operations")
        max_bloc = _entier_positif(taille_bloc_masse_max, "taille_bloc_masse_max")
        # Ce prototype valide seulement de petits blocs physiques.
        if max_bloc > 64:
            raise ValueError("validation dense de masse limitée à des blocs de taille 64")
        _verifier_canonique(r, "R")
        self.r = csr_matrix(r, dtype=float, copy=True)
        _verifier_canonique(self.r, "R converti")
        self.r.sum_duplicates()
        self.r.eliminate_zeros()
        self.r.sort_indices()
        n, m = self.r.shape
        self.n = n
        self.operations_symboliques = self.operations_numeriques = 0
        self.coefficients = 0
        if n == 0 or n != m or not np.all(np.isfinite(self.r.data)):
            raise ValueError("R carré non vide et fini requis")
        self._budget_coefficients(n)
        self.diagonale = self.r.diagonal()
        if np.any(self.diagonale <= 0):
            raise ValueError("diagonale R strictement positive requise")
        for i in range(n):
            a, b = self.r.indptr[i:i+2]
            if np.any(self.r.indices[a:b] < i):
                raise ValueError("R doit être triangulaire supérieur")
        if np.iscomplexobj(echelles):
            raise ValueError("les échelles doivent être réelles")
        self.echelles = (np.ones(n) if echelles is None
                         else np.asarray(echelles, dtype=float).copy())
        if (self.echelles.shape != (n,) or not np.all(np.isfinite(self.echelles))
                or np.any(self.echelles <= 0)):
            raise ValueError("échelles positives finies de taille n requises")
        self.masse, self.bloc_masse_max = self._masse(masse, max_bloc)
        self.voisins = []
        motif = [set() for _ in range(n)]
        self.coefficients = n  # Les diagonales sont toutes nécessaires.

        def ajouter(i, j):
            if i > j:
                i, j = j, i
            if i == j or j in motif[i]:
                return
            self._budget_coefficients(self.coefficients+1)
            motif[i].add(j)
            self.coefficients += 1

        for i in range(n):
            a, b = self.r.indptr[i:i+2]
            voisins = [(int(j), float(v)) for j, v in
                       zip(self.r.indices[a:b], self.r.data[a:b]) if j > i]
            self.voisins.append(voisins)
            for j, _ in voisins:
                ajouter(i, j)
        self.coefficients_motif_r = self.coefficients
        for i in range(n):
            a, b = self.masse.indptr[i:i+2]
            for j in self.masse.indices[a:b]:
                if j > i:
                    ajouter(i, int(j))
        self.coefficients_demandes = self.coefficients
        # Toute dépendance a un indice de ligne strictement supérieur à i.
        # Un seul passage croissant suffit donc à propager la fermeture.
        for i in range(n):
            for j in tuple(motif[i]):
                for k, _ in self.voisins[i]:
                    self._operation(symbolique=1)
                    ajouter(k, j)
        self.valeurs = [dict() for _ in range(n)]
        for i in range(n-1, -1, -1):
            rii = float(self.diagonale[i])
            voisins = self.voisins[i]
            # Off-diagonales d'abord : elles entrent ensuite dans S_ii.
            for j in sorted(motif[i]):
                self._operation(numerique=2*len(voisins)+1)
                somme = math.fsum(v*self._lire(k, j) for k, v in voisins)
                self.valeurs[i][j] = self._fini(-somme/rii)
            self._operation(numerique=2*len(voisins)+3)
            somme = math.fsum(v*self.valeurs[i][k] for k, v in voisins)
            inverse_rii = 1./rii
            sii = self._fini(inverse_rii*inverse_rii-somme/rii)
            if sii <= 0:
                raise ArithmeticError("diagonale inverse non positive : calcul flottant inexploitable")
            self.valeurs[i][i] = sii
        termes = []
        for i in range(n):
            a, b = self.masse.indptr[i:i+2]
            for j, mij in zip(self.masse.indices[a:b], self.masse.data[a:b]):
                if j >= i:
                    self._operation(numerique=5)
                    facteur = 1. if j == i else 2.
                    terme = facteur*mij*self.echelles[i]*self.echelles[j]*self._lire(i, int(j))
                    termes.append(self._fini(float(terme)))
        self.trace_masse = self._fini(math.fsum(termes))
        if self.trace_masse <= 0:
            raise ArithmeticError("trace massique non positive : calcul flottant inexploitable")
        self.lambda_trace = self._fini(1./self.trace_masse)
        if self.lambda_trace <= 0:
            raise ArithmeticError("minoration de trace hors domaine flottant")
        self.facteur_annulation_trace = math.fsum(abs(t) for t in termes)/self.trace_masse
        self.certification_machine = False

    def _masse(self, masse, max_bloc):
        if np.iscomplexobj(masse):
            raise ValueError("la masse doit être réelle")
        if not hasattr(masse, "tocsr") and np.ndim(masse) == 1:
            diagonal = np.asarray(masse, dtype=float)
            if (diagonal.shape != (self.n,) or not np.all(np.isfinite(diagonal))
                    or np.any(diagonal <= 0)):
                raise ValueError("masse diagonale positive finie de taille n requise")
            return diags(diagonal, format="csr"), 1
        _verifier_canonique(masse, "M")
        m = csr_matrix(masse, dtype=float, copy=True)
        _verifier_canonique(m, "M convertie")
        m.sum_duplicates()
        m.eliminate_zeros()
        m.sort_indices()
        if (m.shape != (self.n, self.n) or not np.all(np.isfinite(m.data))
                or (m-m.T).nnz):
            raise ValueError("masse finie exactement symétrique de taille n requise")
        _, labels = connected_components(m, directed=False)
        blocs = {}
        for i, label in enumerate(labels):
            blocs.setdefault(int(label), []).append(i)
        maximum = max(map(len, blocs.values()))
        if maximum > max_bloc:
            raise BudgetInverseSelectionnee(
                "bloc connexe de masse trop grand pour la validation SPD locale",
                dict(phase="validation_masse", taille_bloc=maximum, budget_bloc=max_bloc))
        for ids in blocs.values():
            if len(ids) == 1:
                positif = m[ids[0], ids[0]] > 0
            else:
                positif = np.linalg.eigvalsh(m[ids][:, ids].toarray())[0] > 0
            if not positif:
                raise ValueError("chaque bloc de masse doit être défini positif")
        return m, maximum

    @staticmethod
    def _fini(valeur):
        if not np.isfinite(valeur):
            raise ArithmeticError("inverse sélectionnée hors domaine flottant")
        return float(valeur)

    def _budget_coefficients(self, prochain):
        if prochain > self.budget_coefficients:
            raise BudgetInverseSelectionnee(
                "budget de coefficients de fermeture dépassé",
                dict(phase="fermeture", coefficients_demandes=prochain,
                     budget_coefficients=self.budget_coefficients,
                     operations_symboliques=self.operations_symboliques))

    def _operation(self, *, symbolique=0, numerique=0):
        prochain = self.operations_symboliques+self.operations_numeriques+symbolique+numerique
        if prochain > self.budget_operations:
            raise BudgetInverseSelectionnee(
                "budget d'opérations de l'inverse sélectionnée dépassé",
                dict(phase="fermeture" if symbolique else "calcul",
                     operations_demandees=prochain, budget_operations=self.budget_operations,
                     coefficients=self.coefficients))
        self.operations_symboliques += symbolique
        self.operations_numeriques += numerique

    def _lire(self, i, j):
        if i > j:
            i, j = j, i
        return self.valeurs[i][j]

    def entree(self, i, j):
        """Entrée de (R.T R)^-1 ; refuse une entrée non sélectionnée."""
        i, j = operator.index(i), operator.index(j)
        if not (0 <= i < self.n and 0 <= j < self.n):
            raise IndexError("indice inverse hors domaine")
        try:
            return self._lire(i, j)
        except KeyError as exc:
            raise KeyError("entrée non sélectionnée ; aucune approximation n'est renvoyée") from exc

    def bilan(self):
        return dict(dimension=self.n, nnz_R=int(self.r.nnz),
                    coefficients_selectionnes=self.coefficients,
                    coefficients_motif_R=self.coefficients_motif_r,
                    coefficients_demandes=self.coefficients_demandes,
                    coefficients_fermeture=self.coefficients-self.coefficients_demandes,
                    operations_symboliques=self.operations_symboliques,
                    operations_numeriques=self.operations_numeriques,
                    operations_total=self.operations_symboliques+self.operations_numeriques,
                    taille_bloc_masse_max=self.bloc_masse_max,
                    trace_masse=self.trace_masse, lambda_trace=self.lambda_trace,
                    facteur_annulation_trace=self.facteur_annulation_trace,
                    modele="A_QR = E^-1 R.T R E^-1",
                    minoration_D_original_etablie=False,
                    certification_machine=False)


def borne_trace(r, masse, echelles=None, **budgets):
    """Retourne le bilan de la minoration conditionnelle pour le facteur fourni."""
    return InverseSelectionnee(r, masse, echelles, **budgets).bilan()
