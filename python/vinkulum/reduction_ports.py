"""Réduction matérielle par ports : énergie d'entrée et contrôle du champ.

SciPy est requis. Le modèle analysé est explicitement K=DᵀD, avec masse M,
sans amortissement, charge intérieure ni évolution non linéaire. L'extraction
du noyau fournit la partie matérielle des seules poutres ; elle ne transforme
pas leur facteur en tangente générale sous précontrainte ou en mouvement.
"""
import time

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix

from ._ports.condensation_energie import CondensationEnergie
from ._ports.ports_releves import PortsReleves, AssemblagePorts
from ._ports.champ_interieur import ControleChamp, adapter_champ, observable_locale

__all__ = ["ReductionMaterielle", "AssemblagePorts", "ControleChamp",
           "adapter_champ", "observable_locale", "reduire_poutres"]


def _indices(valeurs, taille, nom, non_vide=True):
    a = np.asarray(valeurs)
    if (a.ndim != 1 or a.dtype.kind not in "iu" or (non_vide and not len(a))
            or np.any(a < 0) or np.any(a >= taille) or len(np.unique(a)) != len(a)):
        raise ValueError(nom+" : indices entiers distincts dans le domaine requis")
    return a.astype(np.intp)


def _csc(donnees):
    nl, nc, ptr, ind, data = donnees
    return csc_matrix((data, ind, ptr), shape=(nl, nc))


class ReductionMaterielle(PortsReleves):
    """Modèle K=DᵀD réduit, en coordonnées physiques choisies par l'appelant.

    La métrique SPD des ports définit leur normalisation et la tolérance du
    Schur. Sans lambda_min fournie, le backend QR calcule un minorant pour D
    d'entrée avec un encadrement en précision étendue. Une constante fournie
    par l'appelant conserve le statut d'hypothèse, quel que soit son nom.
    Les autres majorants restent évalués en doubles et non certifiés.
    """
    def __init__(self, d, m, interieur, interface, metrique, omega_max, *,
                 lambda_min=None, tolerance=1e-8, methode="qr", max_blocs=12,
                 max_directions=128, budget_qr=10_000_000,
                 precision_spectrale=80, budget_spectral=10_000_000):
        debut = time.perf_counter()
        if any(np.iscomplexobj(a) for a in (d, m, metrique)):
            raise ValueError("D, M et la métrique doivent être réels")
        facteur = None
        self.certificat_spectral = None
        if lambda_min is None:
            if methode != "qr":
                raise ValueError("la borne spectrale automatique requiert methode='qr'")
            from ._ports.inverse_selectionnee import _verifier_canonique
            for nom, matrice in (("D", d), ("M", m)):
                _verifier_canonique(matrice, nom)
            # LIL n'expose pas toujours ce drapeau. Sa conversion CSR peut
            # garder des doublons ; les refuser avant les normalisations QR.
            d, m = csr_matrix(d, dtype=float, copy=True), csr_matrix(m, dtype=float, copy=True)
            _verifier_canonique(d, "D")
            _verifier_canonique(m, "M")
            facteur = CondensationEnergie(d, interieur, interface, metrique, budget_qr)
            from ._ports.certificat_spectral import certifier_spectre
            mi = csc_matrix(m)[facteur.i][:, facteur.i]
            self.certificat_spectral = certifier_spectre(
                facteur.di, facteur.r, facteur.echelles, mi,
                precision=precision_spectrale, budget=budget_spectral)
            lambda_min = self.certificat_spectral["lambda_min"]
        super().__init__(d, m, interieur, interface, metrique, lambda_min, omega_max,
                         tolerance=tolerance, max_blocs=max_blocs,
                         max_directions=max_directions, budget_qr=budget_qr,
                         methode=methode, _facteur=facteur)
        self.preparation_s = time.perf_counter()-debut
        self.origine = {"modele": "K=D.T D", "certification_champ_machine": False}
        self.coordonnees_physiques = np.arange(self.m.shape[0])
        self._controle_reponses = None
        self._frequence_reponses = None

    @classmethod
    def depuis_noyau(cls, noyau, libres, interface, metrique, omega_max, *,
                     tolerance_contraintes=1e-12, **options):
        """Réduit la partie matérielle des poutres sur des DDL physiques libres.

        Les indices libres et d'interface visent les coordonnées du noyau
        (tx,ty,tz,rx,ry,rz par corps). Ils doivent éliminer les contraintes
        par suppression de coordonnées ; une contrainte liant plusieurs DDL
        exige une application admissible fournie via le constructeur matriciel.
        tolerance_contraintes borne les équations de position à l'état extrait,
        dans leurs unités natives ; elle ne certifie pas leurs arrondis.
        Les données retournées sont des copies indépendantes de l'état futur.
        """
        debut = time.perf_counter()
        donnees = noyau.facteurs_materiels_poutres()
        for contrainte in donnees["etat"]["contraintes"]:
            if (contrainte["nature"] == "contact_non_lisse" or contrainte["non_holonome"]
                    or contrainte["pilotee"]):
                raise ValueError("extraction réduite : contraintes fixes holonomes requises")
        if not np.isfinite(tolerance_contraintes) or tolerance_contraintes < 0:
            raise ValueError("tolérance de contraintes finie et non négative requise")
        residu_contraintes = float(np.max(np.abs(donnees["phi"]), initial=0.))
        if residu_contraintes > tolerance_contraintes:
            raise ValueError("les contraintes doivent être satisfaites à l'état extrait")
        d, m, g = (_csc(donnees[n]) for n in ("d", "masse", "contraintes"))
        libres = _indices(libres, d.shape[1], "libres")
        ports = _indices(interface, d.shape[1], "interface")
        inverse = {int(v): j for j, v in enumerate(libres)}
        if any(int(v) not in inverse for v in ports):
            raise ValueError("les interfaces doivent appartenir aux coordonnées libres")
        # Une sélection est exactement admissible pour G d'entrée si chaque
        # colonne conservée est nulle. Aucun seuil n'efface une contrainte.
        admissible = g[:, libres]
        if np.any(admissible.data != 0.):
            raise ValueError("coordonnées libres non admissibles : G[:,libres] doit être nul")
        s = np.array([inverse[int(v)] for v in ports], dtype=np.intp)
        i = np.setdiff1d(np.arange(len(libres)), s)
        if not len(i):
            raise ValueError("au moins une coordonnée intérieure requise")
        dl, ml = d[:, libres], m[libres][:, libres]
        # Une permutation d'indices physiques peut désordonner le stockage.
        # Trier ne change aucune valeur et ne somme aucun doublon.
        dl.sort_indices()
        ml.sort_indices()
        reduction = cls(dl, ml, i, s, metrique, omega_max, **options)
        reduction.coordonnees_physiques = libres.copy()
        reduction.origine = dict(modele="partie_materielle_des_poutres",
                                 tangente_globale=False, donnees_natives=donnees,
                                 tolerance_contraintes=tolerance_contraintes,
                                 residu_contraintes_initial=residu_contraintes,
                                 certification_champ_machine=False)
        reduction.preparation_s = time.perf_counter()-debut
        return reduction

    def _force(self, force_physique):
        if np.iscomplexobj(force_physique):
            raise ValueError("force physique de port réelle requise")
        force = np.asarray(force_physique, dtype=float)
        if force.shape != (self.p,) or not np.all(np.isfinite(force)):
            raise ValueError("force physique de port : vecteur fini de taille p requis")
        return force

    def reconstruire(self, omega, port_normalise):
        """Reconstruit le champ réel depuis les coordonnées normalisées des ports."""
        if np.iscomplexobj(port_normalise):
            raise ValueError("coordonnées de port réelles requises")
        return super().reconstruire(omega, port_normalise)

    def _physique(self, omega, resultat, champ=None):
        resultat = dict(resultat)
        y = resultat.pop("deplacement")
        resultat["coordonnees_ports_normalisees"] = y.copy()
        resultat["borne_ports_normalises"] = resultat.pop("borne_ports")
        resultat["normalisation_ports"] = self.qr.w.copy()
        resultat["deplacement_ports_physique"] = self.qr.w@y
        if champ is None:
            # Même ordre d'évaluation que le contrôle groupé. Une seule
            # colonne reste une matrice pour les produits BLAS.
            champ = (self.reconstruire(omega, np.eye(self.p))@y[:, None])[:, 0]
        resultat["champ_physique"] = champ.copy()
        resultat["coordonnees_physiques"] = self.coordonnees_physiques.copy()
        return resultat

    def _preparer_reponses(self, omega):
        if np.iscomplexobj(omega) or np.ndim(omega) != 0:
            raise ValueError("pulsation scalaire réelle requise")
        omega = float(omega)
        if not np.isfinite(omega) or not 0 <= omega <= self.omega_max:
            raise ValueError("fréquence hors bande commune")
        controle = self._controle_reponses
        if controle is None or controle._revisions != (self._revision,):
            controle = ControleChamp(AssemblagePorts([self], [self.qr.w]))
            self._controle_reponses = controle
            self._frequence_reponses = None
        etat = self._frequence_reponses
        if etat is None or etat.omega != omega:
            etat = controle._preparer(omega)
            self._frequence_reponses = etat
        # Un autre lecteur peut publier une autre fréquence entre cette
        # affectation et le retour : conserver l'état propre à cet appel.
        return etat

    def reponse(self, omega, force_physique):
        """Réponse à une force physique de port, avec contrôles masse/déformation.

        Les coordonnées et forces du port suivent l'ordre d'interface fourni.
        La reconstruction suit l'ordre des coordonnées libres. Le contrôle
        est ponctuel à omega ; il peut refuser une marge près d'une résonance.
        """
        force = self._force(force_physique)
        etat = self._preparer_reponses(omega)
        resultats, champs = etat._evaluer((self.qr.w.T@force)[:, None])
        return self._physique(omega, resultats[0], champs[0][:, 0])

    def reponses(self, omega, forces_physiques):
        """Réponses à plusieurs charges physiques indépendantes à la même pulsation.

        forces_physiques est une matrice réelle (nombre_ports, nombre_charges),
        avec au moins une colonne. Chaque dictionnaire renvoyé suit le contrat
        de reponse ; les champs, normes et bornes sont propres à sa colonne.
        Le Schur, sa factorisation, sa marge et les relèvements sont partagés.
        Seule la dernière fréquence est conservée, avec invalidation après
        enrichissement. Les arrondis des champs et majorants restent non certifiés.
        """
        if np.iscomplexobj(forces_physiques):
            raise ValueError("matrice de forces physiques réelle requise")
        forces = np.asarray(forces_physiques, dtype=float)
        if (forces.ndim != 2 or forces.shape[0] != self.p or not forces.shape[1]
                or not np.all(np.isfinite(forces))):
            raise ValueError("forces physiques : matrice finie (ports, charges) requise")
        etat = self._preparer_reponses(omega)
        resultats, champs = etat._evaluer(self.qr.w.T@forces)
        return [self._physique(omega, r, champs[0][:, j]) for j, r in enumerate(resultats)]

    def adapter(self, frequences, force_physique, tolerance_relative=1e-8, max_etapes=8):
        """Enrichit la base pour les deux normes physiques aux fréquences demandées.

        La force est exprimée dans les unités physiques du port. Le facteur
        statique et son certificat spectral sont conservés. L'acceptation
        porte sur les points demandés et les bornes évaluées en doubles.
        """
        force = self._force(force_physique)
        if np.iscomplexobj(frequences):
            raise ValueError("fréquences réelles requises")
        frequences = np.asarray(frequences, dtype=float)
        assemblage = AssemblagePorts([self], [self.qr.w])
        resultat = adapter_champ(assemblage, frequences, self.qr.w.T@force,
                                 tolerance_relative=tolerance_relative, max_etapes=max_etapes)
        resultat.pop("controle")
        resultat["reponses"] = [self._physique(w, rep)
                                  for w, rep in zip(frequences, resultat["reponses"])]
        return resultat


def reduire_poutres(noyau, libres, interface, metrique, omega_max, **options):
    """Construit une réduction de la partie matérielle des poutres du noyau."""
    return ReductionMaterielle.depuis_noyau(noyau, libres, interface, metrique, omega_max, **options)
