"""Assemblage énergétique de compléments : ports partagés, modes privés.

Les marges portent sur tout le bloc conservé global. Aucun Schur local
n'est inversé. Les majorants restent numériques, sans certificat machine.
"""
import warnings
import numpy as np
from scipy.linalg import cholesky,solve_triangular,lu_factor,lu_solve,LinAlgWarning
from .controle_complement import norme


def reel(a,nom,ndim=2):
    if np.iscomplexobj(a):raise ValueError(nom+' réel requis')
    a=np.array(a,dtype=float,copy=True)
    if a.ndim!=ndim or not np.all(np.isfinite(a)):raise ValueError(nom+' fini de dimension correcte requis')
    return a


class EtatLocal:
    """Préparation locale indépendante des charges et de l'inversibilité du Schur."""
    def __init__(self,c,omega):
        r=c.reduction
        if r.base is not c._base:raise ValueError('contrôle périmé après enrichissement')
        omega=r._pulsation(omega);z=omega**2;mu=z/r.lambda_complement
        petit=np.eye(len(c.theta))-mu*c.theta;second=c.a0+mu*c.a1
        self.y=np.linalg.solve(petit,second);residu=second-petit@self.y
        self.erreur_y=norme(residu)/(1-mu*c.norme_theta)
        self.x=r.t.copy();self.x[r.i]-=c.w@self.y
        self.dx=r.d@self.x;mx=r.m@self.x;gm=self.x.T@mx;gd=self.dx.T@self.dx
        self.extension_m=float(np.sqrt(max(np.linalg.eigvalsh((gm+gm.T)*.5)[-1],0.)))
        self.extension_d=c._norme_extension_d(self.y,self.dx,gd)
        s=gd-z*gm;self.s=(s+s.T)*.5
        dm,dd=c.norme_w_m*self.erreur_y,c.norme_w_d*self.erreur_y
        self.cm,self.cd=c.cm+dm,c.cd+dd
        self.defaut=(2*self.extension_d*dd+dd**2+z*(2*self.extension_m*dm+dm**2))
        self.borne=c.borne_schur+self.defaut
        self.c=c;self.alpha=r.lambda_complement-r.omega_max**2
        if not all(np.all(np.isfinite(v)) for v in (self.s,self.x,self.borne,self.cm,self.cd)):
            raise ArithmeticError('état local non fini')
    def erreurs(self,q):
        c=self.c;r=c.reduction;qn=np.linalg.norm(q,axis=0)
        eta=np.minimum(c.delta*qn,c.delta_colonnes@np.abs(q));yq=self.y@q
        action=c.delta*eta/self.alpha+(c.defaut_fonctionnel+self.defaut)*qn
        local={}
        for name,constante,coeff,nw,nh in (
            ('masse',self.cm,1.,c.norme_w_m,c.norme_h_m),
            ('deformation',self.cd,np.sqrt(r.lambda_complement),c.norme_w_d,c.norme_h_d)):
            local[name]=np.minimum(constante*qn,coeff*eta/self.alpha+c._images_h[name].normes(yq)
                                    +(nw+nh)*self.erreur_y*qn)
        return action,local


class AssemblageComplements:
    """A_j : déplacements physiques globaux -> ports physiques de la pièce j.

    Les coordonnées retenues de chaque pièce sont ajoutées indépendamment.
    metrique normalise les déplacements globaux ; forces et champs des ports
    sont physiques. Facteur_externe et masse_externe ajoutent leurs énergies
    aux pièces (pas de masse automatique aux interfaces). Les intérieurs sont
    privés, sans charge intérieure ni couplage direct entre leurs énergies.
    Construire un nouvel assemblage pour modifier connexions ou énergies.
    """
    def __init__(self,controles,applications,metrique,*,facteur_externe=None,masse_externe=None,budget_conserve=2048):
        self.controles=tuple(controles)
        if not self.controles or len(self.controles)!=len(applications):raise ValueError('contrôles et applications non vides appariés requis')
        metric=reel(metrique,'métrique');g=metric.shape[0]
        if not g or metric.shape!=(g,g) or np.any(metric!=metric.T):raise ValueError('métrique carrée symétrique requise')
        l=cholesky(metric,lower=True);self.w=solve_triangular(l.T,np.eye(g),lower=False)
        self.g=g;self.metrique=metric
        self.applications=tuple(reel(a,'application') for a in applications)
        self.taille_conservee=g+sum(c.reduction.r for c in self.controles)
        if isinstance(budget_conserve,(bool,np.bool_)) or not isinstance(budget_conserve,(int,np.integer)) or budget_conserve<1:
            raise ValueError('budget conservé entier positif requis')
        if self.taille_conservee>budget_conserve:raise ValueError('budget du bloc conservé dépassé')
        self.base_ids=tuple(c.reduction.base for c in self.controles)
        self.e=[];self.normes_e=[];offset=g
        for c,a in zip(self.controles,self.applications,strict=True):
            r=c.reduction
            if a.shape!=(r.p,g):raise ValueError('application de forme ports locaux × ports globaux requise')
            e=np.zeros((r.p+r.r,self.taille_conservee))
            e[:r.p,:g]=np.linalg.solve(r.qr.w,a@self.w)
            e[r.p:,offset:offset+r.r]=np.eye(r.r);offset+=r.r
            self.e.append(e);self.normes_e.append(norme(e))
        self.de=np.zeros((0,g)) if facteur_externe is None else reel(facteur_externe,'facteur externe')
        self.me=np.zeros((g,g)) if masse_externe is None else reel(masse_externe,'masse externe')
        if self.de.shape[1]!=g or self.me.shape!=(g,g) or np.any(self.me!=self.me.T):raise ValueError('énergies externes de dimensions compatibles requises')
        if np.linalg.eigvalsh(self.me)[0]<0:raise ValueError('masse externe positive semi-définie requise')
        dn=self.de@self.w;mn=self.w.T@self.me@self.w
        self.ke=dn.T@dn;self.mn=(mn+mn.T)*.5
        self.extensions_externes=dict(masse=float(np.sqrt(max(np.linalg.eigvalsh(self.mn)[-1],0.))),deformation=norme(dn))
        self.omega_max=min(c.reduction.omega_max for c in self.controles)
    def _valider(self,omega):
        for c,base in zip(self.controles,self.base_ids,strict=True):
            if c.reduction.base is not base:raise ValueError('assemblage périmé après enrichissement')
            c.reduction._pulsation(omega)
    def reponses(self,omega,forces):
        self._valider(omega);forces=reel(forces,'forces')
        if forces.shape[0]!=self.g or not forces.shape[1]:raise ValueError('forces de forme ports globaux × charges requises')
        states=[EtatLocal(c,omega) for c in self.controles];n=self.taille_conservee
        shat=np.zeros((n,n));enveloppe=np.zeros_like(shat)
        shat[:self.g,:self.g]=self.ke-float(omega)**2*self.mn
        for state,e in zip(states,self.e,strict=True):
            shat+=e.T@state.s@e;enveloppe+=state.borne*(e.T@e)
        shat=(shat+shat.T)*.5;enveloppe=(enveloppe+enveloppe.T)*.5
        borne=float(np.linalg.eigvalsh(enveloppe)[-1]);sigma=float(np.linalg.svd(shat,compute_uv=False)[-1]);marge=sigma-borne
        rhs=np.zeros((n,forces.shape[1]));rhs[:self.g]=self.w.T@forces
        with warnings.catch_warnings():
            warnings.simplefilter('error',LinAlgWarning)
            try:lu=lu_factor(shat,check_finite=False)
            except LinAlgWarning as exc:raise np.linalg.LinAlgError('bloc global conservé singulier') from exc
        q=lu_solve(lu,rhs,check_finite=False);u=self.w@q[:self.g]
        if not np.all(np.isfinite(q)):raise ArithmeticError('réponse globale non finie')
        qs=[e@q for e in self.e];fields=[s.x@ql for s,ql in zip(states,qs,strict=True)]
        errors=[s.erreurs(ql) for s,ql in zip(states,qs,strict=True)]
        action=np.minimum(borne*np.linalg.norm(q,axis=0),sum(ne*err[0] for ne,err in zip(self.normes_e,errors,strict=True)))
        residual=np.linalg.norm(rhs-shat@q,axis=0)
        bq=(residual+action)/marge if marge>0 else None
        bounds={};local_bounds=[]
        for state,field,err,ne in zip(states,fields,errors,self.normes_e,strict=True):
            r=state.c.reduction;item={}
            for name,extension in (('masse',state.extension_m+state.cm),('deformation',state.extension_d+state.cd)):
                values=np.sqrt(np.maximum(np.sum(field*(r.m@field),axis=0),0.)) if name=='masse' else np.linalg.norm(r.d@field,axis=0)
                absolute=None if bq is None else err[1][name]+extension*ne*bq
                item[name]=dict(normes=values,absolues=absolute)
            local_bounds.append(item)
        for name in ('masse','deformation'):
            external=np.sqrt(np.maximum(np.sum(u*(self.me@u),axis=0),0.)) if name=='masse' else np.linalg.norm(self.de@u,axis=0)
            values=np.sqrt(sum(b[name]['normes']**2 for b in local_bounds)+external**2)
            absolute=None if bq is None else np.sqrt(sum(b[name]['absolues']**2 for b in local_bounds)+(self.extensions_externes[name]*bq)**2)
            relative=[None]*forces.shape[1] if absolute is None else [float(b/(v-b)) if v>b else None for b,v in zip(absolute,values,strict=True)]
            bounds[name]=dict(normes=values,absolues=absolute,relatives=relative)
        return dict(champs=fields,champ_ports=u,coordonnees=q,bornes=bounds,bornes_locales=local_bounds,
            borne_coordonnees=bq,marge=marge,sigma_min_bloc_conserve=sigma,borne_schur=borne,
            bornes_schur_locales=[s.borne for s in states],schur=shat,enveloppe=enveloppe,
            taille_conservee=n,certification_machine=False)
