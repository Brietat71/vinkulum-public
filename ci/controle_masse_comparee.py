"""Contrôle de réponse avec masse connectée, par comparaison de Loewner.

Les équations, Schur, champs et normes physiques emploient toujours M.
Seules les normes auxiliaires utilisent L=diag(M), avec les constantes
prouvées alpha et beta. Aucune certification machine des réponses.
"""
import math
import time

import numpy as np
from scipy.linalg import cholesky,solve_triangular
from scipy.sparse import diags

from comparaison_masse import certifier_comparaison_masse
from controle_complement import norme,enveloppe_bernstein
from controle_facteurs import ControleFacteurs,ImageFaibleRang
from inverse_contrainte_energie import InverseContrainteEnergie
from racine_masse_blocs import RacineMasseBlocs


class ControleMasseComparee(ControleFacteurs):
    """Le coût de la preuve de comparaison fait partie de la préparation.

    Les bornes scalaires sont proposées par l'appelant puis vérifiées sur
    M exact ; elles ne sont jamais déduites des signes d'un facteur approché.
    L'admission requiert alpha L < M < beta L. La restriction aux blocs
    connexes de six coordonnées disparaît, sous les budgets déclarés du
    certificat. Une comparaison trop faible peut refuser le contrôle.
    """
    def __init__(self, reduction, *, bibliotheque, alpha=.25, beta=2., profondeur=8):
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
        self.comparaison_masse = certifier_comparaison_masse(r.mi,alpha,beta,bibliotheque=bibliotheque)
        self.alpha_masse,self.beta_masse = alpha,beta
        massi = RacineMasseBlocs(diags(self.comparaison_masse['diagonale'],format='csr'))
        primal,dual_scale = math.sqrt(beta),1/math.sqrt(alpha)
        dual = InverseContrainteEnergie(massi.r, np.ones(len(r.i)), r.b)
        self.h = r.c@np.linalg.solve(r.b.T@r.c, r.b.T@self.w)
        self.wbar = self.w-self.h
        self.dw, self.dh = r.di@self.w, r.di@self.h
        self.mw, self.mh = primal*(massi.r@self.w), primal*(massi.r@self.h)
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
            dual_scale*dual.dual_contraint(e0), dual_scale*dual.dual_contraint(c1),
            dual_scale*dual.dual_contraint(f) if f.shape[1] else np.empty((len(r.i)-r.r, 0)),
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
        self.extension = 'gram'
        u = np.finfo(float).eps/2
        n,p = r.d.shape[0],r.t.shape[1]
        if max(n,p)*u >= 1:
            raise ArithmeticError('dimensions hors modèle d’arrondi')
        self._gram_gamma = n*u/(1-n*u)
        self._somme_gamma = p*u/(1-p*u)
        self.facteurs = dict(extension='gram',gamma_gram=self._gram_gamma,
            gamma_sommes=self._somme_gamma,certification_machine=False)
        j = np.linalg.solve(r.b.T@r.c,r.b.T@self.w)
        self._images_h = dict(masse=ImageFaibleRang(self.mh,primal*(massi.r@r.c),j),
                             deformation=ImageFaibleRang(self.dh,r.di@r.c,j))
        self.facteurs['reparation'] = {k:v.diagnostic for k,v in self._images_h.items()}
        self.preparation_s = time.perf_counter()-debut

