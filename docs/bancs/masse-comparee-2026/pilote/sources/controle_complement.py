"""Enveloppe résiduelle contrainte et contrôle du bloc conservé complet.

Majorants démontrés en arithmétique exacte, évalués ici en binary64 sans
encadrement des arrondis. Une minoration spectrale dirigée fournie à la
réduction ne transforme pas ces évaluations en certificats machine.
"""
import math
import time
import warnings

import numpy as np
from scipy.linalg import cholesky, solve_triangular, lu_factor, lu_solve, LinAlgWarning

from inverse_contrainte_energie import InverseContrainteEnergie
from racine_masse_blocs import RacineMasseBlocs


def norme(a):
    return float(np.linalg.norm(a, 2)) if a.size else 0.


def enveloppe_bernstein(e0, c1, f, theta, d, rho, profondeur=8):
    """Majore ||E0+μC1+μ²F(I−μΘ)^−1D|| pour 0<=μ<=rho.

    Les matrices sont déjà exprimées dans la norme duale contrainte.
    La partie polynomiale est bornée par ses coefficients Bernstein ;
    une queue de résolvante borne les puissances restantes.
    """
    tnorm = norme(theta)
    marge = 1-rho*tnorm
    if marge <= 0:
        raise ArithmeticError("résolvante auxiliaire hors marge positive")
    coefficients = [e0.copy(), rho*c1]
    puissance = f.copy()
    termes = []
    for q in range(2, profondeur+1):
        degre = q-1
        beta = []
        beta_colonnes = []
        for j in range(degre+1):
            x = np.zeros_like(e0)
            for k in range(j+1):
                x += (math.comb(j, k)/math.comb(degre, k))*coefficients[k]
            beta.append(norme(x))
            beta_colonnes.append(np.linalg.norm(x, axis=0))
        polynome = max(beta)
        queue = rho**q*norme(puissance)*norme(d)/marge
        colonnes = (np.max(beta_colonnes, axis=0)
                    +rho**q*norme(puissance)*np.linalg.norm(d, axis=0)/marge)
        termes.append(dict(degre=degre, polynome=polynome, queue=queue,
                           majorant=polynome+queue, majorants_colonnes=colonnes.tolist()))
        coefficients.append(rho**q*(puissance@d))
        puissance = puissance@theta
    return min(t["majorant"] for t in termes), termes


class ControleComplement:
    """Préparation uniforme ; réponse et marge évaluées à chaque fréquence.

    C*=C0(B.T C0)^−1 est le right-inverse théorique pour réparer les
    contraintes. H=C*(B.T W) donne Wbar=W−H. La réponse garde W et l'audit
    ajoute explicitement la différence de champ et de fonctionnel.
    Les valeurs flottantes de ces opérations restent non certifiées.
    """
    def __init__(self, reduction, profondeur=8):
        debut = time.perf_counter()
        r = reduction
        self.reduction = r
        self._base = r.base
        if isinstance(profondeur, bool) or not isinstance(profondeur, int) or profondeur < 2:
            raise ValueError("profondeur entière au moins égale à deux requise")
        if r.v.shape[1]:
            k = (r.dv.T@r.dv)
            l = cholesky((k+k.T)*.5, lower=True, check_finite=False)
            self.w = solve_triangular(l, r.v.T, lower=True, check_finite=False).T
        else:
            self.w = r.v.copy()
        massi = RacineMasseBlocs(r.mi)
        dual = InverseContrainteEnergie(massi.r, np.ones(len(r.i)), r.b)
        self.h = r.c@np.linalg.solve(r.b.T@r.c, r.b.T@self.w)
        self.wbar = self.w-self.h
        self.dw, self.dh = r.di@self.w, r.di@self.h
        self.mw, self.mh = massi.r@self.w, massi.r@self.h
        self.norme_w_m, self.norme_w_d = norme(self.mw), norme(self.dw)
        self.norme_h_m, self.norme_h_d = norme(self.mh), norme(self.dh)
        self.rho = r.omega_max**2/r.lambda_complement
        theta = r.lambda_complement*(self.w.T@(r.mi@self.w))
        self.theta = (theta+theta.T)*.5
        self.norme_theta = norme(self.theta)
        self.a0 = self.w.T@r.g0
        self.a1 = -r.lambda_complement*(self.w.T@r.g1)
        self.resolvante_marge = 1-self.rho*self.norme_theta
        if self.resolvante_marge <= 0:
            raise ArithmeticError("résolvante réduite auxiliaire non séparée de zéro")
        self.ymax = (norme(self.a0)+self.rho*norme(self.a1))/self.resolvante_marge
        p = r.di.T@(r.di@self.wbar)
        q = r.lambda_complement*(r.mi@self.wbar)
        e0 = r.g0-p@self.a0
        e1 = -r.lambda_complement*r.g1-p@self.a1
        f = q-p@self.theta
        c1 = e1+f@self.a0
        d = self.theta@self.a0+self.a1
        self.delta, self.termes = enveloppe_bernstein(
            dual.dual_contraint(e0), dual.dual_contraint(c1),
            dual.dual_contraint(f) if f.shape[1] else np.empty((len(r.i)-r.r, 0)),
            self.theta, d, self.rho, profondeur)
        self.delta_colonnes = np.min([v["majorants_colonnes"] for v in self.termes], axis=0)
        alpha = r.lambda_complement-r.omega_max**2
        self.defaut_masse = self.norme_h_m*self.ymax
        self.defaut_deformation = self.norme_h_d*self.ymax
        self.cm = self.delta/alpha+self.defaut_masse
        self.cd = np.sqrt(r.lambda_complement)*self.delta/alpha+self.defaut_deformation
        # Expansions bilinéaires avant la norme : le produit de normes
        # globales perdrait les orthogonalités entre W et le correctif H.
        mih = r.mi@self.h
        mixte_k = r.dt.T@self.dh
        mixte_m = r.mt[r.i].T@self.h
        quadratique_k = self.dh.T@self.dh-self.dw.T@self.dh-self.dh.T@self.dw
        wh = self.w.T@mih
        quadratique_m = self.h.T@mih-wh-wh.T
        self.defaut_fonctionnel = (
            2*(norme(mixte_k)+r.omega_max**2*norme(mixte_m))*self.ymax
            +(norme(quadratique_k)+r.omega_max**2*norme(quadratique_m))*self.ymax**2)
        self.borne_schur = self.delta**2/alpha+self.defaut_fonctionnel
        if not np.all(np.isfinite([self.delta, self.cm, self.cd, self.borne_schur])):
            raise ArithmeticError("enveloppe de contrôle non finie")
        self.preparation_s = time.perf_counter()-debut

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
        erreur_y = norme(residu_y)/(1-mu*self.norme_theta)
        x = r.t.copy()
        x[r.i] -= self.w@y
        dx, mx = r.d@x, r.m@x
        gram_m = x.T@mx
        extension_m = float(np.sqrt(max(np.linalg.eigvalsh((gram_m+gram_m.T)*.5)[-1], 0.)))
        s = dx.T@dx-z*gram_m
        s = (s+s.T)*.5
        # La résolution du petit Y n'est pas supposée exacte. Ce défaut
        # ajoute un champ W deltaY et son écart de fonctionnel énergétique.
        dm, dd = self.norme_w_m*erreur_y, self.norme_w_d*erreur_y
        cm, cd = self.cm+dm, self.cd+dd
        defaut_coefficients = (2*norme(dx)*dd+dd**2
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
                ("deformation", cd, norme(dx)+cd, self.dh, self.norme_h_d+self.norme_w_d, np.sqrt(r.lambda_complement))):
            valeurs = (np.sqrt(np.maximum(np.sum(champ*(r.m@champ), axis=0), 0.))
                       if nom == "masse" else np.linalg.norm(r.d@champ, axis=0))
            local = (facteur_residu*eta/alpha+np.linalg.norm(h_facteur@yq, axis=0)
                     +norme_hw*erreur_y*qnorm)
            local = np.minimum(local, constante*qnorm)
            absolues = local+extension*bq if bq is not None else None
            relatives = [None]*forces.shape[1] if absolues is None else [
                float(b/(n-b)) if n > b else None for b, n in zip(absolues, valeurs)]
            bornes[nom] = dict(normes=valeurs, absolues=absolues, relatives=relatives)
        return dict(champ=champ, coordonnees=q, bornes=bornes,
                    marge=marge, sigma_min_bloc_conserve=sigma,
                    borne_schur=borne_schur, borne_coordonnees=bq,
                    residu_coefficients=norme(residu_y), erreur_coefficients=erreur_y,
                    certification_machine=False, taille_conservee=r.p+r.r,
                    taille_complement_reduit=r.taille_complement_reduit)
