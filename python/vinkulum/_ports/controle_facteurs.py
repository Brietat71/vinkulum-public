"""Contrôle expérimental : grande norme d'extension et réparation factorisée.

Le solve et les normes des champs physiques gardent le chemin historique.
Les majorants et les défauts supplémentaires sont évalués en binary64 :
aucune certification machine des réponses n'est ajoutée.
"""
import time
import warnings

import numpy as np
from scipy.linalg import qr, lu_factor, lu_solve, LinAlgWarning

from .controle_complement import ControleComplement, norme
from .racine_masse_blocs import RacineMasseBlocs


def normes_majorees(a, axis=0):
    """Normes euclidiennes avec échelle et budget de sommes de carrés.

    Le modèle mixte |fl(t)-t|<=u|t|+tau/2 couvre les sous-flux graduels
    de la réduction normalisée ; tau est le plus petit sous-normal.
    Les opérations scalaires finales sont arrondies vers l'extérieur
    par successeur/prédécesseur. Cela ne certifie pas les entrées a.
    """
    a = np.abs(np.asarray(a,dtype=float))
    if not np.all(np.isfinite(a)):
        raise ArithmeticError("norme d'audit non finie")
    n = a.shape[axis]
    if n == 0:
        return np.sum(a,axis=axis)
    if n == 1:
        return np.squeeze(a,axis=axis)
    echelle = np.max(a,axis=axis)
    denom = np.where(echelle>0,echelle,1.)
    normalise = a/np.expand_dims(denom,axis)
    somme = np.sum(normalise*normalise,axis=axis)
    u = np.finfo(float).eps/2
    nu = np.nextafter((n+3)*u,np.inf)
    if nu >= .5:
        raise ArithmeticError("dimension hors budget de norme")
    gamma = np.nextafter(nu/np.nextafter(1-nu,-np.inf),np.inf)
    absolu = np.nextafter((3*n+2)*np.nextafter(0.,1.),np.inf)
    total = np.nextafter(somme+absolu,np.inf)
    quotient = np.nextafter(total/np.nextafter(1-gamma,-np.inf),np.inf)
    racine = np.nextafter(np.sqrt(quotient),np.inf)
    resultat = np.nextafter(echelle*racine,np.inf)
    resultat = np.where(echelle==0,0.,resultat)
    if not np.all(np.isfinite(resultat)):
        raise ArithmeticError("majorant de norme non représentable")
    return resultat


class ImageFaibleRang:
    """Majore l'action de l'image stockée, défauts de facteurs conservés.

    Les facteurs fournis ne sont jamais supposés reproduire exactement
    l'image, même si elle provient d'un produit de faible rang.
    """
    def __init__(self, image, gauche, droite):
        eps = np.finfo(float).eps/2
        def gamma(k):
            if k*eps >= 1:
                raise ArithmeticError("dimension hors modèle d'arrondi")
            return k*eps/(1-k*eps)
        self.image = np.asarray(image, dtype=float)
        q, r = qr(gauche, mode="economic", check_finite=False)
        self.t = r@droite
        aq = np.abs(q)
        gram = q.T@q
        audit = (np.abs(gram-np.eye(q.shape[1]))/(1-eps)
                 +gamma(q.shape[0])*(aq.T@aq))
        self.kappa = np.sqrt(1+float(np.max(np.sum(audit,axis=1),initial=0.)))
        approx = q@self.t
        defaut = np.abs(self.image-approx)
        calcul = gamma(q.shape[1])*(aq@np.abs(self.t))+eps/(1-eps)*defaut
        # Encadre aussi l'écart de l'ancien produit flottant image@v à
        # l'action exacte de l'image stockée, sous le modèle d'arrondi.
        self.gamma = gamma(self.image.shape[1])
        self.ec = normes_majorees(defaut+calcul)+self.gamma*normes_majorees(self.image)
        if not all(np.all(np.isfinite(v)) for v in (self.t,self.ec,self.kappa)):
            raise ArithmeticError("audit de l'image de réparation non fini")
        self.diagnostic = dict(rang_facteur=self.t.shape[0],
            colonnes_image=self.image.shape[1],defaut_colonnes=self.ec.tolist(),
            kappa=float(self.kappa),certification_machine=False)

    def normes(self, v):
        av = np.abs(v)
        action = self.t@v
        calcul = self.gamma*(np.abs(self.t)@av)
        return self.kappa*(normes_majorees(action)+normes_majorees(calcul))+self.ec@av


class ControleFacteurs(ControleComplement):
    """Majore la grande norme via Gram ou QR, puis les images de réparation.

    Les formules d'arrondi supposent l'absence de dépassement et de
    sous-flux non couvert. Leurs évaluations flottantes restent des
    diagnostics ; le contrôleur ne certifie pas ses réponses en machine.
    """
    def __init__(self, reduction, profondeur=8, *, extension="gram"):
        debut = time.perf_counter()
        if extension not in ("gram", "qr"):
            raise ValueError("extension gram ou qr requise")
        self.extension = extension
        super().__init__(reduction, profondeur)
        r = reduction
        if extension == "qr":
            u = np.column_stack((r.dt, self.dw))
            q, t = qr(u, mode="economic", check_finite=False)
            p = r.t.shape[1]
            self._td, self._wd = t[:, :p], t[:, p:]
            eps = np.finfo(float).eps/2
            def gamma(k):
                if k*eps >= 1:
                    raise ArithmeticError("dimension hors modèle d'arrondi")
                return k*eps/(1-k*eps)
            # Les deux défauts sont évalués en doubles : ce sont des audits
            # numériques des majorations du carnet, pas des intervalles.
            aq, ar = np.abs(q), np.abs(t)
            produit = q@t
            defaut = u-produit
            erreur_qr = gamma(q.shape[1])*(aq@ar)+eps/(1-eps)*np.abs(defaut)
            ec = normes_majorees(np.abs(defaut)+erreur_qr)
            self._ec0, self._ecw = ec[:p], ec[p:]
            gram = q.T@q
            erreur_gram = gamma(q.shape[0])*(aq.T@aq)
            orthogonalite = np.abs(gram-np.eye(q.shape[1]))/(1-eps)+erreur_gram
            self._kappa_d = np.sqrt(1+float(np.max(np.sum(orthogonalite, axis=1), initial=0.)))
            # Ecart entre fl(D fl(T-E fl(WY))) et fl(DT)-fl(DiW)Y.
            ad, adi = abs(r.d).tocsr(), abs(r.di).tocsr()
            a = ad@np.abs(r.t)
            ai = adi@np.abs(r.t[r.i])
            bw = adi@np.abs(self.w)
            gd = gamma(int(np.max(np.diff(ad.indptr), initial=0)))
            gi = gamma(int(np.max(np.diff(adi.indptr), initial=0)))
            gw = gamma(self.w.shape[1])
            kappa = gw+eps*(1+gw)
            h0 = 2*gd*a+eps*(1+gd)*ai
            h1 = (gd*(1+kappa)+kappa+gi)*bw
            self._h0, self._hw = normes_majorees(h0), normes_majorees(h1)
            self._gamma_w, self._eps = gw, eps
            if not all(np.all(np.isfinite(v)) for v in
                       (self._ec0,self._ecw,self._h0,self._hw,self._kappa_d)):
                raise ArithmeticError("audit du facteur énergétique non fini")
            self.facteurs = dict(extension="qr", dimension_image=q.shape[1],
                defaut_qr_colonnes=ec.tolist(), kappa=float(self._kappa_d),
                formation_constante=self._h0.tolist(), formation_lineaire=self._hw.tolist(),
                certification_machine=False)
        else:
            u = np.finfo(float).eps/2
            n, p = r.d.shape[0], r.t.shape[1]
            if max(n,p)*u >= 1:
                raise ArithmeticError("dimensions hors modèle d'arrondi")
            self._gram_gamma = n*u/(1-n*u)
            self._somme_gamma = p*u/(1-p*u)
            self.facteurs = dict(extension="gram",gamma_gram=self._gram_gamma,
                gamma_sommes=self._somme_gamma,certification_machine=False)
        j = np.linalg.solve(r.b.T@r.c, r.b.T@self.w)
        rm = RacineMasseBlocs(r.mi).r
        self._images_h = dict(masse=ImageFaibleRang(self.mh,rm@r.c,j),
                              deformation=ImageFaibleRang(self.dh,r.di@r.c,j))
        self.facteurs["reparation"] = {k:v.diagnostic for k,v in self._images_h.items()}
        self.preparation_s = time.perf_counter()-debut

    def _norme_extension_d(self, y, dx=None, gram=None):
        if self.extension == "gram":
            if gram is None:
                if dx is None:
                    r = self.reduction
                    x = r.t.copy()
                    x[r.i] -= self.w@y
                    dx = r.d@x
                gram = dx.T@dx
            diagonale = np.diag(gram)
            if np.any(diagonale < np.finfo(float).tiny):
                # Un Gram nul/subnormal ne prouve pas une image nulle.
                # Ce repli compte une norme physique, sans changer le solve.
                if dx is None:
                    r = self.reduction
                    x = r.t.copy()
                    x[r.i] -= self.w@y
                    dx = r.d@x
                return float(normes_majorees(np.asarray(dx).ravel()))
            trace = float(np.sum(diagonale))
            gf, gs = self._gram_gamma, self._somme_gamma
            frobenius2 = trace/((1-gf)*(1-gs))
            ag = np.abs(gram)
            lignes = float(np.max(np.sum(ag,axis=1),initial=0.))/(1-gs)
            colonnes = float(np.max(np.sum(ag,axis=0),initial=0.))/(1-gs)
            # Les deux majorations concernent seulement la plus grande
            # valeur propre du Gram exact, jamais une petite action q.T Gq.
            return float(np.sqrt(min(frobenius2,max(lignes,colonnes)+gf*frobenius2)))
        ay = np.abs(y)
        produit = self._wd@y
        petit = self._td-produit
        # Défaut du petit produit et de sa soustraction, avant la norme.
        calcul = (self._gamma_w*(np.abs(self._wd)@ay)
                  +self._eps*(np.abs(self._td)+np.abs(produit)))
        eps_qr = float(normes_majorees(self._ec0+self._ecw@ay))
        eps_formation = float(normes_majorees(self._h0+self._hw@ay))
        return float(self._kappa_d*(norme(petit)+float(normes_majorees(calcul.ravel())))
                     +eps_qr+eps_formation)

    def reponses(self, omega, forces):
        """Schur énergétique de taille ports+r, avec défaut de contrainte borné."""
        r = self.reduction
        if r.base is not self._base:
            raise ValueError("contrôle périmé après enrichissement")
        omega = r._pulsation(omega)
        if np.iscomplexobj(forces):
            raise ValueError("forces réelles requises")
        forces = np.asarray(forces, dtype=float)
        if (forces.ndim != 2 or forces.shape[0] != r.p or not forces.shape[1]
                or not np.all(np.isfinite(forces))):
            raise ValueError("forces finies de forme (ports, charges) requises")
        z, mu = omega**2, omega**2/r.lambda_complement
        petit = np.eye(len(self.theta))-mu*self.theta
        second = self.a0+mu*self.a1
        y = np.linalg.solve(petit, second)
        residu_y = second-petit@y
        norme_residu_y = norme(residu_y)
        erreur_y = norme_residu_y/(1-mu*self.norme_theta)
        x = r.t.copy()
        x[r.i] -= self.w@y
        dx, mx = r.d@x, r.m@x
        gram_m = x.T@mx
        extension_m = float(np.sqrt(max(np.linalg.eigvalsh((gram_m+gram_m.T)*.5)[-1], 0.)))
        gram_d = dx.T@dx
        s = gram_d-z*gram_m
        s = (s+s.T)*.5
        # La résolution du petit Y n'est pas supposée exacte. Ce défaut
        # ajoute un champ W deltaY et son écart de fonctionnel énergétique.
        dm, dd = self.norme_w_m*erreur_y, self.norme_w_d*erreur_y
        cm, cd = self.cm+dm, self.cd+dd
        norme_dx = self._norme_extension_d(y, dx, gram_d)
        defaut_coefficients = (2*norme_dx*dd+dd**2
                              +z*(2*extension_m*dm+dm**2))
        borne_schur = self.borne_schur+defaut_coefficients
        sigma = float(np.linalg.svd(s, compute_uv=False)[-1])
        marge = sigma-borne_schur
        rhs = np.vstack((r.qr.w.T@forces, np.zeros((r.r, forces.shape[1]))))
        with warnings.catch_warnings():
            warnings.simplefilter("error", LinAlgWarning)
            try:
                lu = lu_factor(s, check_finite=False)
            except LinAlgWarning as exc:
                raise np.linalg.LinAlgError("bloc conservé singulier") from exc
        q = lu_solve(lu, rhs, check_finite=False)
        champ = x@q
        if not np.all(np.isfinite(champ)):
            raise ArithmeticError("réponse non finie")
        residu = np.linalg.norm(rhs-s@q, axis=0)
        qnorm = np.linalg.norm(q, axis=0)
        eta = np.minimum(self.delta*qnorm, self.delta_colonnes@np.abs(q))
        alpha = r.lambda_complement-r.omega_max**2
        action_schur = (self.delta*eta/alpha
                        +(self.defaut_fonctionnel+defaut_coefficients)*qnorm)
        bq = (residu+action_schur)/marge if marge > 0 else None
        bornes = {}
        # Former les champs avant les normes évite l'annulation des
        # petites combinaisons de charges dans un Gram projeté.
        yq = y@q
        for nom, constante, extension, h_facteur, norme_hw, facteur_residu in (
                ("masse", cm, extension_m+cm, self.mh, self.norme_h_m+self.norme_w_m, 1.),
                ("deformation", cd, norme_dx+cd, self.dh, self.norme_h_d+self.norme_w_d, np.sqrt(r.lambda_complement))):
            valeurs = (np.sqrt(np.maximum(np.sum(champ*(r.m@champ), axis=0), 0.))
                       if nom == "masse" else np.linalg.norm(r.d@champ, axis=0))
            local = (facteur_residu*eta/alpha+self._images_h[nom].normes(yq)
                     +norme_hw*erreur_y*qnorm)
            local = np.minimum(local, constante*qnorm)
            absolues = local+extension*bq if bq is not None else None
            relatives = [None]*forces.shape[1] if absolues is None else [
                float(b/(n-b)) if n > b else None for b, n in zip(absolues, valeurs)]
            bornes[nom] = dict(normes=valeurs, absolues=absolues, relatives=relatives)
        return dict(champ=champ, coordonnees=q, bornes=bornes,
                    marge=marge, sigma_min_bloc_conserve=sigma,
                    borne_schur=borne_schur, borne_coordonnees=bq,
                    residu_coefficients=norme_residu_y, erreur_coefficients=erreur_y,
                    certification_machine=False, taille_conservee=r.p+r.r,
                    taille_complement_reduit=r.taille_complement_reduit)
