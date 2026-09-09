"""Champs de référence Decimal pour une chaîne de blocs physiques <= 6.

Le modèle visé est D_float.T D_float - omega² M_float, avec les valeurs
binary64 interprétées exactement. Le calcul Decimal a une précision finie :
70/90 chiffres testent sa convergence, sans certificat d'arrondi.
Les blocs de raideur sont assemblés une seule fois. Chaque fréquence
condense les intérieurs et conserve les transferts pour toutes les charges
appliquées au dernier bloc. Aucune matrice globale dense n'est construite.
"""
from decimal import localcontext
import operator
import time

import numpy as np
from scipy.sparse import csr_matrix, issparse

from reference_ports_precision import (_ZERO, _UN, _decimal, _precision,
                                       _resout_bloc, _resultat, _zeros)


def _csr_binary64(a, nom):
    """Inspection des indices réels avant conversion ; CSR/CSC seulement."""
    if not issparse(a) or a.format not in ("csr", "csc"):
        raise ValueError(nom+" doit être fournie explicitement en CSR ou CSC canonique")
    if a.dtype != np.dtype("float64") or not np.all(np.isfinite(a.data)):
        raise ValueError(nom+" doit contenir des valeurs binary64 réelles finies")
    majeur, mineur = a.shape if a.format == "csr" else a.shape[::-1]
    ptr, ids = a.indptr, a.indices
    if (len(ptr) != majeur+1 or ptr[0] != 0 or ptr[-1] != len(ids)
            or len(ids) != len(a.data) or np.any(ptr[1:] < ptr[:-1])
            or np.any(ids < 0) or np.any(ids >= mineur)):
        raise ValueError(nom+" : stockage sparse invalide")
    for j in range(majeur):
        ligne = ids[ptr[j]:ptr[j+1]]
        if np.any(ligne[1:] <= ligne[:-1]):
            raise ValueError(nom+" : doublons ou indices non triés interdits")
    # La conversion entre CSR/CSC canoniques ne somme aucun doublon.
    return csr_matrix(a, copy=True)


class OracleChamps:
    """Oracle pour blocs contigus, masse diagonale et port terminal de taille p.

    Chaque ligne de D touche au plus deux blocs voisins. L'ordre physique
    des coordonnées est conservé ; p appartient à [1,6]. M peut être un
    vecteur diagonal binary64 ou une matrice diagonale CSR/CSC canonique.
    Les forces sont appliquées exclusivement au dernier bloc, et peuvent
    contenir un nombre arbitraire de colonnes. Une seule fréquence est
    mémorisée, pour éviter de retenir O(nombre_frequences*n*p²) données.

    L'élimination exige des pivots de blocs non singuliers. En particulier,
    elle peut refuser un système global inversible dont un pivot intérieur
    est nul ; aucune permutation globale à travers les blocs n'est tentée.
    Une fréquence exactement singulière est refusée, même à force nulle.
    Un pivot très petit n'est pas un certificat de précision : il faut
    comparer plusieurs précisions et, si besoin, refuser cette référence.
    """

    def __init__(self, d, m, p=6, dps=70):
        debut = time.perf_counter()
        if isinstance(p, (bool, np.bool_)):
            raise ValueError("p doit être un entier de 1 à 6")
        try:
            p = operator.index(p)
        except TypeError as exc:
            raise ValueError("p doit être un entier de 1 à 6") from exc
        if not 1 <= p <= 6:
            raise ValueError("p doit être un entier de 1 à 6")
        self._contexte = _precision(dps)
        self.dps, self.p = self._contexte.prec, p
        d = _csr_binary64(d, "D")
        self.total = d.shape[1]
        if self.total == 0 or self.total % p:
            raise ValueError("le nombre de coordonnées doit être un multiple positif de p")
        self.nombre_blocs = self.total//p
        if issparse(m):
            m = _csr_binary64(m, "M")
            if m.shape != (self.total, self.total):
                raise ValueError("M doit avoir la dimension physique de D")
            for i in range(self.total):
                a, b = m.indptr[i:i+2]
                if any(j != i and v != 0 for j, v in zip(m.indices[a:b], m.data[a:b])):
                    raise ValueError("M doit être diagonale dans cet oracle")
            masses = m.diagonal()
        else:
            masses = np.asarray(m)
            if masses.dtype != np.dtype("float64") or masses.shape != (self.total,):
                raise ValueError("masse vectorielle binary64 de dimension physique requise")
        if not np.all(np.isfinite(masses)) or np.any(masses <= 0):
            raise ValueError("masses diagonales strictement positives et finies requises")
        self._masses = [_decimal(x, "masse", True) for x in masses]
        self._diagonaux = [_zeros(p) for _ in range(self.nombre_blocs)]
        self._liaisons = [_zeros(p) for _ in range(self.nombre_blocs-1)]
        self.compteurs = dict(assemblages_raideur=1, produits_gram=0,
                              condensations_frequence=0, resolutions_blocs=0,
                              reponses=0, colonnes_force=0)
        with localcontext(self._contexte):
            for row in range(d.shape[0]):
                a, b = d.indptr[row:row+2]
                actifs = [(int(j), _decimal(x, "D")) for j, x in
                          zip(d.indices[a:b], d.data[a:b]) if x != 0]
                if actifs and actifs[-1][0]//p-actifs[0][0]//p > 1:
                    raise ValueError("une ligne de D relie des blocs non voisins")
                for debut_pair, (ja, va) in enumerate(actifs):
                    bloc_a, ca = divmod(ja, p)
                    for jb, vb in actifs[debut_pair:]:
                        bloc_b, cb = divmod(jb, p)
                        cible = (self._diagonaux[bloc_a] if bloc_a == bloc_b
                                 else self._liaisons[bloc_a])
                        cible[ca][cb] += va*vb
                        self.compteurs["produits_gram"] += 1
            for diagonal in self._diagonaux:
                for i in range(p):
                    for j in range(i+1, p):
                        diagonal[j][i] = diagonal[i][j]
        self._omega = None
        self._transferts = self._schur = self._inverse_schur = None
        self.certification_machine = False
        self.preparation_s = time.perf_counter()-debut

    def _preparer_frequence(self, omega):
        if np.iscomplexobj(omega):
            raise ValueError("fréquence réelle requise")
        omega = _decimal(omega, "omega")
        if omega < 0:
            raise ValueError("omega doit être positif ou nul")
        if omega == self._omega:
            return
        z, p = omega*omega, self.p

        def diagonal(i):
            bloc = [row.copy() for row in self._diagonaux[i]]
            for j in range(p):
                bloc[j][j] -= z*self._masses[i*p+j]
            return bloc

        pivot = diagonal(0)
        transferts = []
        for i, liaison in enumerate(self._liaisons):
            try:
                resolu = _resout_bloc(pivot, liaison)
            except ValueError as exc:
                raise ValueError(f"pivot intérieur de bloc {i} nul à omega={omega}") from exc
            suivant = diagonal(i+1)
            for j in range(p):
                for k in range(j, p):
                    jk = sum((liaison[l][j]*resolu[l][k] for l in range(p)), _ZERO)
                    kj = sum((liaison[l][k]*resolu[l][j] for l in range(p)), _ZERO)
                    suivant[j][k] -= (jk+kj)/2
                    suivant[k][j] = suivant[j][k]
            transferts.append(resolu)
            pivot = suivant
        identite = [[_UN if i == j else _ZERO for j in range(p)] for i in range(p)]
        try:
            inverse = _resout_bloc(pivot, identite)
        except ValueError as exc:
            raise ValueError(f"Schur terminal singulier à omega={omega}") from exc
        self._transferts, self._schur, self._inverse_schur = transferts, pivot, inverse
        self._omega = omega
        self.compteurs["condensations_frequence"] += 1
        self.compteurs["resolutions_blocs"] += self.nombre_blocs

    def schur(self, omega, retour_decimal=False):
        """Schur terminal p×p ; un Schur singulier est explicitement refusé."""
        with localcontext(self._contexte):
            self._preparer_frequence(omega)
            return _resultat(self._schur, retour_decimal)

    def reponse(self, omega, force_ports_matrix, retour_decimal=False):
        """Champ total×q, ou vecteur total si la force est un vecteur p.

        retour_decimal=True conserve les nombres Decimal calculés dans un
        tuple (de tuples si plusieurs forces). Les données sont des copies.
        """
        if np.iscomplexobj(force_ports_matrix):
            raise ValueError("forces réelles requises")
        force = np.asarray(force_ports_matrix)
        vectoriel = force.ndim == 1
        if vectoriel:
            force = force[:, None]
        if force.ndim != 2 or force.shape[0] != self.p or force.shape[1] == 0:
            raise ValueError("force de port p ou p×q avec q positif requise")
        q, p = force.shape[1], self.p
        with localcontext(self._contexte):
            rhs = [[_decimal(force[i, j], "force") for j in range(q)] for i in range(p)]
            self._preparer_frequence(omega)
            # L'inverse terminale est seulement p×p et est calculée une
            # fois par fréquence. Le nombre q de charges n'est pas limité
            # au nombre p de colonnes attendu par _resout_bloc.
            courant = [[sum((self._inverse_schur[i][k]*rhs[k][j] for k in range(p)), _ZERO)
                         for j in range(q)] for i in range(p)]
            blocs = [courant]
            for transfert in reversed(self._transferts):
                courant = [[-sum((transfert[i][k]*courant[k][j] for k in range(p)), _ZERO)
                            for j in range(q)] for i in range(p)]
                blocs.append(courant)
            resultat = [row for bloc in reversed(blocs) for row in bloc]
            self.compteurs["reponses"] += 1
            self.compteurs["colonnes_force"] += q
            rendu = _resultat(resultat, retour_decimal)
            if vectoriel:
                return tuple(row[0] for row in rendu) if retour_decimal else rendu[:, 0]
            return rendu
