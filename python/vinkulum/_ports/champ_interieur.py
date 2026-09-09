"""Contrôle physique d'une réponse assemblée et enrichissement des intérieurs.

Enveloppes locales uniformes sur la bande ; marge globale et tolérance
relative évaluées aux seules fréquences demandées. Les modèles sont linéaires
conservatifs. Constantes spectrales et arrondis ne sont pas certifiés.
"""
import operator
import time
import warnings

import numpy as np
from scipy.linalg import solve, lu_factor, lu_solve, LinAlgWarning
from scipy.sparse.linalg import splu

from .ports_krylov import _sym, _norme_gram


def borne_relative(borne, norme_candidate):
    """Majore l'erreur relative à la solution exacte, si sa norme est séparée de zéro."""
    if borne is None or norme_candidate <= borne:
        return None
    return float(borne/(norme_candidate-borne))


def observable_locale(reduction, omega, port, covecteur, dual_approche):
    """Correction d'une observable linéaire à port fixé, par résidu primal-dual.

    Le coût des produits physiques et du solve massique est explicite ; le
    dual est fourni par l'appelant. L'observable est un covecteur intérieur,
    donc ses unités et sa normalisation ont un sens physique à conserver.
    """
    r = reduction
    x = r.reconstruire(omega, port)
    ell, dual = np.asarray(covecteur, dtype=float), np.asarray(dual_approche, dtype=float)
    ni = len(r.qr.i)
    if (x.ndim != 1 or ell.shape != (ni,) or dual.shape != (ni,)
            or not np.all(np.isfinite(ell)) or not np.all(np.isfinite(dual))):
        raise ValueError("port, observable et dual vectoriels finis requis")
    mii = r.m[r.qr.i][:, r.qr.i]
    z = float(omega)**2
    residu = r.qr.di.T@(r.qr.d@x)-z*(r.m@x)[r.qr.i]
    residu_dual = ell-r.qr.di.T@(r.qr.di@dual)+z*(mii@dual)
    mf = splu(mii)
    def norme(v):
        return np.sqrt(max(float(v@mf.solve(v)), 0.))
    eta = norme(residu)
    alpha = r.lambda_min-z
    valeur = float(ell@x[r.qr.i])
    correction = -float(dual@residu)
    return dict(valeur_candidate=valeur, correction=correction,
                valeur_corrigee=valeur+correction,
                borne_sans_correction=float(norme(ell)*eta/alpha),
                borne_apres_correction=float(norme(residu_dual)*eta/alpha),
                certification_machine=False)


class ControleChamp:
    """Hypothèses : masses locales PSD, intérieurs coercifs et partition conforme.

    Les masses locales et les constantes de coercivité restent des hypothèses
    fournies par le modèle ; ce contrôleur n'en établit pas un certificat.
    Les deux énergies externes font l'objet d'un contrôle numérique de signe.
    """
    def __init__(self, assemblage):
        self.assemblage = assemblage
        self.structures = assemblage.sous_structures
        self.applications = assemblage.applications
        # Les erreurs d'intérieur ont une trace nulle. Leurs normes s'ajoutent
        # donc sans double compter les degrés de liberté des interfaces.
        self.masse, self.deformation = [], []
        self.erreur_schur = np.zeros_like(assemblage.k_externe)
        for r, e in zip(self.structures, self.applications):
            a = r.lambda_min-r.omega_max**2
            cm = np.sqrt(r.audit_bande["majorant_residu"]/a)
            self.masse.append(float(cm))
            self.deformation.append(float(np.sqrt(r.lambda_min)*cm))
            self.erreur_schur += r.borne_uniforme*(e.T@e)
        self.borne_schur = float(np.linalg.norm(self.erreur_schur, 2))
        # K_ext et M_ext décrivent ici des énergies additionnelles positives.
        # Leur positivité est contrôlée numériquement, pas certifiée.
        self._racines_externes = []
        for nom, a in (("masse", assemblage.m_externe), ("raideur", assemblage.k_externe)):
            if not np.any(a):
                racine = np.zeros((0, a.shape[0]))
            else:
                vals, vecteurs = np.linalg.eigh(_sym(a))
                if vals[0] < 0:
                    raise ValueError(nom+" externe non positive : norme physique indéfinie")
                racine = np.sqrt(vals)[:, None]*vecteurs.T
            self._racines_externes.append(racine)
        self._dimensions = tuple(r.taille_interieure for r in self.structures)
        self._revisions = tuple(r._revision for r in self.structures)
        self._identites = tuple(id(r) for r in self.structures)
        self._applications = [e.copy() for e in self.applications]
        self._externes = (assemblage.k_externe.copy(), assemblage.m_externe.copy())
        self._extensions_locales = []
        for constantes in (self.masse, self.deformation):
            empile = np.vstack([c*e for c, e in zip(constantes, self.applications)])
            self._extensions_locales.append(float(np.linalg.norm(empile, 2)))

    def _verifier(self):
        if (self._dimensions != tuple(r.taille_interieure for r in self.structures)
                or self._revisions != tuple(r._revision for r in self.structures)):
            raise ValueError("contrôle périmé après enrichissement : reconstruire ControleChamp")
        a = self.assemblage
        if (self._identites != tuple(id(r) for r in a.sous_structures)
                or len(self._applications) != len(a.applications)
                or any(not np.array_equal(e, ancien)
                       for e, ancien in zip(a.applications, self._applications))
                or not np.array_equal(a.k_externe, self._externes[0])
                or not np.array_equal(a.m_externe, self._externes[1])):
            raise ValueError("contrôle périmé après modification de l'assemblage")

    def _preparer(self, omega):
        self._verifier()
        if np.iscomplexobj(omega) or np.ndim(omega) != 0:
            raise ValueError("pulsation scalaire réelle requise")
        return _FrequenceChamp(self, float(omega))

    def reponse(self, omega, force):
        """Une charge ; les opérations physiques restent évaluées sur le champ."""
        if np.iscomplexobj(force):
            raise ValueError("un vecteur de force réel est requis")
        force = np.asarray(force, dtype=float)
        if force.shape != (self.assemblage.k_externe.shape[0],) or not np.all(np.isfinite(force)):
            raise ValueError("un vecteur de force fini est requis")
        resultats, _ = self._preparer(omega)._evaluer(force[:, None])
        return resultats[0]

    def reponses(self, omega, forces):
        """Colonnes de charges globales ; un dictionnaire indépendant par colonne."""
        if np.iscomplexobj(forces):
            raise ValueError("matrice de forces réelle requise")
        forces = np.asarray(forces, dtype=float)
        if (forces.ndim != 2 or forces.shape[0] != self.assemblage.k_externe.shape[0]
                or not forces.shape[1] or not np.all(np.isfinite(forces))):
            raise ValueError("matrice de forces finie de forme (ports, charges) requise")
        resultats, _ = self._preparer(omega)._evaluer(forces)
        return resultats


class _FrequenceChamp:
    """Données réutilisables à une fréquence et une révision du modèle.

    Une matrice locale nombre_DDL×nombre_ports est conservée par structure,
    sans la dilater à la taille du port global. Les Gram servent seulement
    aux normes d'opérateur maximales, jamais à la norme du champ demandé.
    """
    def __init__(self, controle, omega):
        a = controle.assemblage
        self.controle, self.omega = controle, omega
        self.s = a.schur(omega)
        self.sigma = float(np.linalg.svd(self.s, compute_uv=False)[-1])
        self.marge = self.sigma-controle.borne_schur
        with warnings.catch_warnings():
            warnings.simplefilter("error", LinAlgWarning)
            try:
                self.facteur = lu_factor(self.s, check_finite=False)
            except LinAlgWarning as exc:
                raise np.linalg.LinAlgError("Schur réduit singulier") from exc
        gm, gd = a.m_externe.copy(), a.k_externe.copy()
        self.relevements = []
        for r, e in zip(controle.structures, controle.applications):
            x = r.reconstruire(omega, np.eye(r.p))
            dx = r.qr.d@x
            gm += e.T@(x.T@(r.m@x))@e
            gd += e.T@(dx.T@dx)@e
            self.relevements.append(x)
        self.extensions = (_norme_gram(gm)+controle._extensions_locales[0],
                           _norme_gram(gd)+controle._extensions_locales[1])

    def _evaluer(self, forces):
        c = self.controle
        c._verifier()
        y = lu_solve(self.facteur, forces, check_finite=False)
        if not np.all(np.isfinite(y)):
            raise ArithmeticError("réponse réduite non finie")
        residus = np.linalg.norm(forces-self.s@y, axis=0)
        bports = ((residus+c.borne_schur*np.linalg.norm(y, axis=0))/self.marge
                  if self.marge > 0 else None)
        normes = [np.linalg.norm(racine@y, axis=0) for racine in c._racines_externes]
        locaux = [np.zeros(y.shape[1]), np.zeros(y.shape[1])]
        champs = []
        for r, e, x, cm, cd in zip(c.structures, c.applications, self.relevements,
                                   c.masse, c.deformation):
            ports = e@y
            champ = x@ports
            champs.append(champ)
            # Former d'abord le champ préserve les petites combinaisons :
            # y.T (X.T M X) y peut s'annuler par arrondi alors que X y != 0.
            nm = np.sqrt(np.maximum(np.sum(champ*(r.m@champ), axis=0), 0.))
            nd = np.linalg.norm(r.qr.d@champ, axis=0)
            normes[0] = np.hypot(normes[0], nm)
            normes[1] = np.hypot(normes[1], nd)
            nport = np.linalg.norm(ports, axis=0)
            locaux[0] = np.hypot(locaux[0], cm*nport)
            locaux[1] = np.hypot(locaux[1], cd*nport)
        resultats = []
        for j in range(y.shape[1]):
            by = float(bports[j]) if bports is not None else None
            result = dict(deplacement=y[:, j].copy(), borne_ports=by,
                          residu_ports=float(residus[j]), sigma_min_reduit=self.sigma,
                          marge=self.marge, borne_schur=c.borne_schur,
                          certification_machine=False)
            for k, nom in enumerate(("masse", "deformation")):
                norme, local, extension = float(normes[k][j]), float(locaux[k][j]), self.extensions[k]
                borne = local+extension*by if by is not None else None
                result[nom] = dict(norme_candidate=norme, borne_absolue=borne,
                                   borne_relative=borne_relative(borne, norme),
                                   defaut_interieur=local, amplification_ports=extension)
            resultats.append(result)
        return resultats, champs


def adapter_champ(assemblage, frequences, force, tolerance_relative=1e-8,
                  max_etapes=8, normes=("masse", "deformation")):
    """Enrichit tous les intérieurs, sans nouvelle factorisation statique.

    L'arrêt utilise les bornes physiques, jamais une solution de référence.
    Aucune garantie entre les fréquences d'acceptation n'est déduite de cet
    échantillonnage ; les enveloppes locales, elles, sont uniformes.
    Une ronde complète traite toutes les sous-structures : le marquage
    sélectif et son optimalité ne sont pas revendiqués par ce prototype.
    """
    if (not np.isfinite(tolerance_relative) or tolerance_relative <= 0
            or isinstance(max_etapes, (bool, np.bool_)) or operator.index(max_etapes) < 0):
        raise ValueError("tolérance positive et budget entier non négatif requis")
    frequences = np.asarray(frequences, dtype=float)
    if (frequences.ndim != 1 or not len(frequences) or not np.all(np.isfinite(frequences))
            or np.any(frequences < 0) or np.any(frequences > assemblage.omega_max)):
        raise ValueError("fréquences non vides dans la bande commune requises")
    if not normes or any(n not in ("masse", "deformation") for n in normes):
        raise ValueError("normes physiques inconnues")
    debut = time.perf_counter()
    historique = []
    statut = "budget_etapes"
    for etape in range(max_etapes+1):
        controle = ControleChamp(assemblage)
        reponses = [controle.reponse(w, force) for w in frequences]
        erreurs = [r[n]["borne_relative"] for r in reponses for n in normes]
        atteint = all(e is not None and e <= tolerance_relative for e in erreurs)
        historique.append(dict(etape=etape,
                                directions=sum(r.taille_interieure for r in controle.structures),
                                borne_schur=controle.borne_schur,
                                bornes_relatives={n: [r[n]["borne_relative"] for r in reponses]
                                                  for n in normes}, accepte=atteint))
        if atteint:
            statut = "tolerance_aux_frequences_demandees"
            break
        if etape == max_etapes:
            break
        changements = [r.enrichir() for r in controle.structures]
        if not any(changements):
            statut = "stagnation_ou_budget_directions"
            break
    # Rendre cohérentes les anciennes métriques de l'assemblage mutable.
    assemblage.erreur_matrice = controle.erreur_schur.copy()
    assemblage.borne_uniforme = controle.borne_schur
    return dict(controle=controle, reponses=reponses, historique=historique,
                statut=statut, preparation_s=time.perf_counter()-debut,
                tolerance_relative=tolerance_relative,
                certification_machine=False)
