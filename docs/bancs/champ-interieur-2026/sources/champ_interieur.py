"""Contrôle physique d'une réponse assemblée et enrichissement des intérieurs.

Enveloppes locales uniformes sur la bande ; marge globale et tolérance
relative évaluées aux seules fréquences demandées. Les modèles sont linéaires
conservatifs. Constantes spectrales et arrondis ne sont pas certifiés.
"""
import operator
import time

import numpy as np
from scipy.linalg import solve
from scipy.sparse.linalg import splu

from ports_krylov import _sym, _norme_gram


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
        for nom, a in (("masse", assemblage.m_externe), ("raideur", assemblage.k_externe)):
            vals = np.linalg.eigvalsh(_sym(a))
            if vals[0] < 0:
                raise ValueError(nom+" externe non positive : norme physique indéfinie")
        self._dimensions = tuple(r.taille_interieure for r in self.structures)

    def reponse(self, omega, force):
        if self._dimensions != tuple(r.taille_interieure for r in self.structures):
            raise ValueError("contrôle périmé après enrichissement : reconstruire ControleChamp")
        a = self.assemblage
        force = np.asarray(force, dtype=float)
        if force.shape != (a.k_externe.shape[0],) or not np.all(np.isfinite(force)):
            raise ValueError("un vecteur de force fini est requis")
        s = a.schur(omega)
        y = solve(s, force, assume_a="sym")
        residu = float(np.linalg.norm(force-s@y))
        sigma = float(np.linalg.svd(s, compute_uv=False)[-1])
        marge = sigma-self.borne_schur
        by = ((residu+self.borne_schur*np.linalg.norm(y))/marge) if marge > 0 else None
        gm, gd = a.m_externe.copy(), a.k_externe.copy()
        poids_m, poids_d = np.zeros_like(gm), np.zeros_like(gd)
        for r, e, cm, cd in zip(self.structures, self.applications, self.masse, self.deformation):
            x = r.reconstruire(omega, np.eye(r.p))@e
            dx = r.qr.d@x
            gm += x.T@(r.m@x)
            gd += dx.T@dx
            poids_m += cm**2*(e.T@e)
            poids_d += cd**2*(e.T@e)
        result = dict(deplacement=y, borne_ports=by, residu_ports=residu,
                      sigma_min_reduit=sigma, marge=marge,
                      borne_schur=self.borne_schur, certification_machine=False)
        for nom, gram, poids in (("masse", gm, poids_m), ("deformation", gd, poids_d)):
            norme = float(np.sqrt(max(float(y@gram@y), 0.)))
            local = float(np.sqrt(max(float(y@poids@y), 0.)))
            extension = _norme_gram(gram)+_norme_gram(poids)
            borne = local+extension*by if by is not None else None
            result[nom] = dict(norme_candidate=norme, borne_absolue=borne,
                               borne_relative=borne_relative(borne, norme),
                               defaut_interieur=local, amplification_ports=extension)
        return result


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
