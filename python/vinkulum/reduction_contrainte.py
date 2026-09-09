"""Réduction fréquentielle avec directions retenues et complément certifié.

Modèle réel K=DᵀD, masse couplée, charges aux ports, sans amortissement.
Les certificats concernent les coefficients binary64 stockés ; les champs
et leurs majorants restent numériques. SciPy est requis.
"""
from copy import deepcopy
import time

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import LinearOperator, eigsh

from ._ports.assemblage_complements import AssemblageComplements
from ._ports.comparaison_masse import BibliothequeComparaison, certifier_comparaison_masse
from ._ports.condensation_energie import CondensationEnergie
from ._ports.controle_masse_comparee import ControleMasseComparee
from ._ports.inertie_binaire_compile import BibliothequeInertie, certifier_inertie_compilee
from ._ports.inverse_selectionnee import _entier_positif
from ._ports.krylov_contraint import KrylovContraint
from ._ports.matrices_certificat import matrice_entree
from ._ports.ordre_separateurs import ordre_separateurs
from ._ports.validation_certificat import InertieImpossible as CertificatIndisponible, _empreinte
from .reduction_ports import ReductionMaterielle

__all__ = ['ReductionContrainte', 'AssemblageContraint', 'CertificatIndisponible']


def _scalaire(x, nom, minimum=0.):
    if np.iscomplexobj(x) or np.ndim(x) != 0:
        raise ValueError(nom+' scalaire réel requis')
    x = float(x)
    if not np.isfinite(x) or x < minimum:
        raise ValueError(nom+' fini dans le domaine requis')
    return x


def _masse_positive(m, alpha, beta, budgets):
    """Preuve PSD : bloc actif SPD et lignes de masse nulle exactement nulles."""
    n = m.shape[0]
    if not n or m.shape != (n, n) or np.any((m-m.T).data != 0):
        raise ValueError('masse carrée exactement symétrique requise')
    diag = m.diagonal()
    if np.any(diag < 0):
        raise ValueError('diagonale massique négative')
    actifs = np.flatnonzero(diag > 0)
    nuls = np.flatnonzero(diag == 0)
    if m[nuls].nnz:
        raise ValueError('une diagonale massique nulle exige une ligne exactement nulle')
    proof = None
    if len(actifs):
        proof = certifier_comparaison_masse(
            m[actifs][:, actifs].tocsr(), alpha, beta,
            bibliotheque=BibliothequeComparaison(), **budgets)
    return dict(masse_sha256=_empreinte(m), dimension=n,
                indices_actifs=actifs.tolist(), indices_sans_masse=nuls.tolist(),
                comparaison_active=proof, certification_machine=True,
                portee='masse complète PSD ; noyau aligné sur les coordonnées sans masse')


def _selection(qr, mi, rang):
    root = np.sqrt(mi.diagonal())
    sc = diags(1/root)
    a = (sc@mi@sc).tocsr()
    def mult(v):
        return (qr.produit(np.asarray(v).reshape(-1, 1)/root[:, None])/root[:, None]).ravel()
    def inv(v):
        return (root[:, None]*qr.solve(root[:, None]*np.asarray(v).reshape(-1, 1))).ravel()
    k = LinearOperator(mi.shape, matvec=mult, dtype=float)
    ki = LinearOperator(mi.shape, matvec=inv, dtype=float)
    seed = 107+len(root)
    values, vec = eigsh(a, M=k, Minv=ki, k=rang, which='LA', tol=1e-11,
                        ncv=min(len(root), max(20, 2*rang+1)),
                        v0=np.random.default_rng(seed).normal(size=len(root)))
    phi = vec/root[:, None]
    if np.any(values <= 0) or not np.all(np.isfinite(phi)):
        raise ArithmeticError('directions retenues proposées invalides')
    return phi, dict(methode='eigsh_generalise_masse_complete', graine=seed,
                     valeurs_inverses=values.tolist(), certification_machine=False)


def _statut(resultat):
    # Un certificat spectral ne transforme jamais une marge flottante en preuve.
    resultat['statut_controle'] = ('borne_numerique_disponible' if resultat['marge'] > 0
                                   else 'marge_non_positive')
    return resultat


class ReductionContrainte:
    """Conserve les basses directions et certifie le complément au seuil gamma.

    interieur/interface forment une partition complète ; metrique est SPD.
    omega_max et sqrt(gamma) sont en rad/s, avec gamma > omega_max².
    Par défaut gamma=4 omega_max² (1 si la bande se réduit à zéro).
    modes_retenus propose des directions par eigsh ; directions_retenues
    fournit alternativement une matrice (intérieurs, rang). Dans les deux
    cas, le B=M_ii Phi réellement stocké est vérifié par inertie dirigée.

    alpha_masse diag(M) < M < beta_masse diag(M) doit être démontrable sur
    le bloc actif de la masse complète. Des ports sans masse sont permis
    si leurs lignes sont exactement nulles. Les autres masses singulières
    sont refusées. Les budgets s'appliquent à chaque certificat, pas à une
    limite globale de mémoire ou de temps ; le QR a son propre budget.
    """
    def __init__(self, d, m, interieur, interface, metrique, omega_max, *,
                 modes_retenus=1, directions_retenues=None, gamma=None,
                 blocs=4, max_directions=128, profondeur=8,
                 alpha_masse=.25, beta_masse=2., budget_qr=10_000_000,
                 budget_operations=100_000_000, budget_coefficients=2_000_000,
                 budget_rectangulaire=2_000_000):
        debut = time.perf_counter()
        omega_max = _scalaire(omega_max, 'omega_max')
        if omega_max > np.sqrt(np.finfo(float).max/4):
            raise ValueError('bande hors domaine flottant')
        gamma = (4*omega_max**2 if omega_max else 1.) if gamma is None else _scalaire(gamma, 'gamma')
        if gamma <= omega_max**2:
            raise ValueError('gamma strictement supérieur à omega_max² requis')
        alpha_masse = _scalaire(alpha_masse, 'alpha_masse')
        beta_masse = _scalaire(beta_masse, 'beta_masse')
        if not 0 < alpha_masse < beta_masse:
            raise ValueError('0 < alpha_masse < beta_masse requis')
        budgets = {name: _entier_positif(value, name) for name, value in (
            ('budget_operations', budget_operations), ('budget_coefficients', budget_coefficients))}
        budget_rectangulaire = _entier_positif(budget_rectangulaire, 'budget_rectangulaire')
        budget_qr = _entier_positif(budget_qr, 'budget_qr')
        blocs = _entier_positif(blocs, 'blocs')
        max_directions = _entier_positif(max_directions, 'max_directions')
        profondeur = _entier_positif(profondeur, 'profondeur')
        if profondeur < 2:
            raise ValueError('profondeur au moins égale à deux requise')
        modes_retenus = _entier_positif(modes_retenus, 'modes_retenus')
        if max(*budgets.values(), budget_rectangulaire) > np.iinfo(np.int64).max:
            raise ValueError('budgets des certificats représentables en int64 requis')
        d, m = matrice_entree(d, 'D'), matrice_entree(m, 'M')
        if m.shape != (d.shape[1], d.shape[1]):
            raise ValueError('D et M de dimensions compatibles requis')
        if np.iscomplexobj(metrique):
            raise ValueError('métrique réelle requise')
        metric = np.array(metrique, dtype=float, copy=True)
        if metric.ndim != 2 or metric.shape[0] != metric.shape[1] or np.any(metric != metric.T):
            raise ValueError('métrique carrée exactement symétrique requise')
        self._masse = _masse_positive(m, alpha_masse, beta_masse, budgets)
        qr = CondensationEnergie(d, interieur, interface, metric, budget_qr)
        ni = len(qr.i)
        mi = m[qr.i][:, qr.i].tocsr()
        di = d[:, qr.i].tocsr()
        # Une sélection physique valide peut réordonner les colonnes CSR.
        # Trier ces copies ne somme aucun doublon et ne change aucun coefficient.
        mi.sort_indices()
        di.sort_indices()
        if np.any(mi.diagonal() <= 0):
            raise ValueError('masse intérieure définie positive requise')
        if directions_retenues is None:
            rang = modes_retenus
            phi = None
        else:
            if np.iscomplexobj(directions_retenues):
                raise ValueError('directions retenues réelles requises')
            phi = np.array(directions_retenues, dtype=float, copy=True)
            if phi.ndim != 2 or phi.shape[0] != ni or not np.all(np.isfinite(phi)):
                raise ValueError('directions retenues finies de forme (intérieurs, rang) requises')
            rang = phi.shape[1]
        if not 0 < rang < ni:
            raise ValueError('rang retenu strictement entre zéro et le nombre d’intérieurs requis')
        if ni*rang > budget_rectangulaire or ni+rang > budgets['budget_coefficients']:
            raise CertificatIndisponible('budget des contraintes dépassé', dict(dimension=ni, rang=rang))
        if phi is None:
            phi, self._selection = _selection(qr, mi, rang)
        else:
            self._selection = dict(methode='directions_fournies', certification_machine=False)
        b = mi@phi
        self._certificat = certifier_inertie_compilee(
            di, mi, b, gamma, bibliotheque=BibliothequeInertie(),
            permutation=ordre_separateurs(di, mi),
            budget_rectangulaire=budget_rectangulaire, **budgets)
        reduction = KrylovContraint(qr, m, b, phi, gamma, omega_max,
                                    blocs=blocs, max_directions=max_directions)
        self._controle = ControleMasseComparee(reduction, bibliotheque=BibliothequeComparaison(),
            alpha=alpha_masse, beta=beta_masse, profondeur=profondeur, **budgets)
        self.omega_max = omega_max
        self.preparation_s = time.perf_counter()-debut
        self.origine = dict(modele='K=D.T D', certification_champ_machine=False)
        self.coordonnees_physiques = np.arange(m.shape[0])

    @classmethod
    def depuis_noyau(cls, noyau, libres, interface, metrique, omega_max, *,
                     tolerance_contraintes=1e-12, **options):
        """Extrait les facteurs matériels des poutres sur des coordonnées admissibles."""
        return ReductionMaterielle.depuis_noyau.__func__(
            cls, noyau, libres, interface, metrique, omega_max,
            tolerance_contraintes=tolerance_contraintes, **options)

    @property
    def certificats(self):
        """Copies des preuves de masse et du complément ; aucune preuve de réponse."""
        return deepcopy(dict(masse_complete=self._masse, complement=self._certificat,
                             comparaison_interieure=self._controle.comparaison_masse))

    @property
    def diagnostic(self):
        """Dimensions, proposition spectrale et coût de préparation."""
        r = self._controle.reduction
        return deepcopy(dict(selection=self._selection, taille_physique=r.m.shape[0],
            ports=r.p, directions_retenues=r.r, taille_conservee=r.p+r.r,
            taille_complement_reduit=r.taille_complement_reduit,
            preparation_s=self.preparation_s, certification_machine_reponses=False))

    def reponses(self, omega, forces):
        """Dictionnaire groupé : champ (DDL, charges), normes et bornes numériques."""
        return _statut(self._controle.reponses(omega, forces))


class AssemblageContraint:
    """Assemble des réductions : ports physiques partagés et directions privées.

    applications[j] transforme les déplacements globaux en ports de la pièce j.
    Les énergies externes s'ajoutent à celles des pièces, sans masse implicite.
    Une singularité locale n'impose pas le refus ; un bloc global singulier
    l'impose. Une marge globale non positive rend les bornes indisponibles.
    """
    def __init__(self, reductions, applications, metrique, *, facteur_externe=None,
                 masse_externe=None, budget_conserve=2048,
                 alpha_masse=.25, beta_masse=2., budget_operations=100_000_000,
                 budget_coefficients=2_000_000):
        debut = time.perf_counter()
        reductions = tuple(reductions)
        if not reductions or not all(isinstance(r, ReductionContrainte) for r in reductions):
            raise TypeError('réductions contraintes non vides requises')
        alpha_masse = _scalaire(alpha_masse, 'alpha_masse')
        beta_masse = _scalaire(beta_masse, 'beta_masse')
        if not 0 < alpha_masse < beta_masse:
            raise ValueError('0 < alpha_masse < beta_masse requis')
        budgets = {name: _entier_positif(value, name) for name, value in (
            ('budget_operations', budget_operations), ('budget_coefficients', budget_coefficients))}
        if max(budgets.values()) > np.iinfo(np.int64).max:
            raise ValueError('budgets représentables en int64 requis')
        self._masse_externe = None
        if masse_externe is not None:
            me = matrice_entree(masse_externe, 'masse externe')
            self._masse_externe = _masse_positive(me, alpha_masse, beta_masse, budgets)
            masse_externe = me.toarray()
        self._assemblage = AssemblageComplements(
            [r._controle for r in reductions], applications, metrique,
            facteur_externe=facteur_externe, masse_externe=masse_externe,
            budget_conserve=budget_conserve)
        self.omega_max = self._assemblage.omega_max
        self.preparation_s = time.perf_counter()-debut

    @property
    def certificat_masse_externe(self):
        """Copie de la preuve PSD externe ; None si aucune masse n'est ajoutée."""
        return deepcopy(self._masse_externe)

    def reponses(self, omega, forces):
        """Champs locaux, ports globaux et bornes sur tout le bloc conservé."""
        return _statut(self._assemblage.reponses(omega, forces))
