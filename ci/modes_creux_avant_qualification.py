"""Modes locaux contraints par inverse décalé, sans base dense admissible.

Le problème est K_s x + Gᵀη = λ M x, Gx=0, K_s=(K+Kᵀ)/2.
Comme `N.modes`, il porte sur la partie symétrique de la raideur locale ;
il ne remplace pas le spectre amorti/gyroscopique ni une étude de stabilité.
`analyse` conserve les valeurs propres signées proches du décalage choisi,
leurs formes normalisées en masse, les réactions et les résidus.
"""
import operator

import numpy as np


def _csc(export):
    from scipy.sparse import csc_matrix
    nr,nc,ptr,idx,data=export
    return csc_matrix((data,idx,ptr),shape=(nr,nc))


def _norme(a):
    return float(np.max(np.asarray(abs(a).sum(axis=1)),initial=0.))


def _normalise_lignes(b):
    b=b.tocsr(copy=True)
    scales=np.asarray(abs(b).max(axis=1).toarray()).ravel() if b.shape[0] else np.zeros(0)
    scales[scales==0]=1.
    b.data/=np.repeat(scales,np.diff(b.indptr))
    norms=np.sqrt(np.asarray(b.multiply(b).sum(axis=1))).ravel()
    norms[norms==0]=1.
    b.data/=np.repeat(norms,np.diff(b.indptr))
    if not np.isfinite(b.data).all():raise ValueError('normalisation des contraintes non finie')
    return b.tocsc(),scales,norms


class _Contraintes:
    """Projection orthogonale : QR mince ou système augmenté, repli de rang."""
    def __init__(self,b):
        from scipy.linalg import qr
        from scipy.sparse import bmat,eye
        from scipy.sparse.linalg import splu
        self.original,self.scales,self.norms=_normalise_lignes(b)
        self.n=b.shape[1];m=b.shape[0]
        self.q=None;self.lu=None;self.repli=False
        self.indices=np.arange(m)
        self.b=self.original
        self.rank=m
        self.seuil=64*np.finfo(float).eps*max(b.shape,default=1)
        if not m:return
        # Les petites familles utilisent une base des normales seulement,
        # jamais le complément Z de dimension n × (n-rang).
        if m>64 and m<self.n:
            try:
                system=bmat([[eye(self.n),self.b.T],[self.b,None]],format='csc')
                lu=splu(system)
                if np.min(abs(lu.U.diagonal()))>np.sqrt(self.seuil):
                    self.lu=lu
                    return
            except RuntimeError:
                pass
            self.repli=True
        elif m>64:
            self.repli=True
        q,r,piv=qr(self.b.T.toarray(),mode='economic',pivoting=True)
        diag=abs(np.diag(r));threshold=self.seuil*max(float(np.max(diag,initial=0.)),1.)
        rank=int(np.count_nonzero(diag>threshold))
        # Refuser une décision de rang non séparée, plutôt que supprimer
        # une contrainte presque dépendante au seuil sans le signaler.
        if np.any((diag>threshold/4)&(diag<4*threshold)):
            raise RuntimeError('rang des contraintes ambigu au seuil numérique')
        self.rank=rank;self.indices=np.asarray(piv[:rank]);self.q=q[:,:rank]
        self.b=self.original[self.indices].tocsc()
        error=np.linalg.norm(self.original@self.q@self.q.T-self.original.toarray(),ord=np.inf) if rank else _norme(self.original)
        if error>8*self.seuil*max(1.,_norme(self.original)):
            raise RuntimeError('réduction des contraintes non vérifiée')

    def projette(self,v):
        if not self.rank:return v.copy()
        if self.q is not None:
            out=v-self.q@(self.q.T@v)
            return out-self.q@(self.q.T@out)
        rhs=np.concatenate((v,np.zeros((self.rank,)+v.shape[1:])))
        out=self.lu.solve(rhs)[:self.n]
        # Corriger la fermeture sans utiliser le déplacement comme second
        # membre de réaction : raffinement du système de projection complet.
        for _ in range(2):
            violation=self.b@out
            if np.max(abs(violation),initial=0.)<=8*np.finfo(float).eps*max(1.,float(np.max(abs(out),initial=0.))):break
            rhs=np.concatenate((np.zeros_like(out),-violation))
            out+=self.lu.solve(rhs)[:self.n]
        return out

    def reactions(self,force):
        if not self.rank:return np.zeros((self.original.shape[0],)+force.shape[1:])
        if self.q is not None:
            small=self.b@self.q
            reduced=np.linalg.solve(small.T,self.q.T@force)
        else:
            rhs=np.concatenate((force,np.zeros((self.rank,)+force.shape[1:])))
            reduced=self.lu.solve(rhs)[self.n:]
        out=np.zeros((self.original.shape[0],)+force.shape[1:])
        out[self.indices]=reduced
        return out


def analyse(noyau,combien=6,t=None,decalage=None,tol=1e-8,maxiter=1000):
    """Valeurs propres signées proches du décalage, formes et contrôles.

    `decalage` est λ en s⁻², pas une fréquence en Hz. None choisit un
    décalage négatif proportionnel à la norme de la raideur en métrique
    de masse, pour admettre les corps libres. Les résultats sont triés
    par distance à ce décalage ; ils ne certifient pas le bas du spectre
    entier en présence de valeurs propres négatives plus éloignées.
    `tol` borne l'erreur arrière du problème symétrique et les contraintes,
    pas l'erreur relative sur une fréquence très petite ou quasi multiple.
    Les dépendances générales peuvent imposer un QR dense des normales.
    SciPy est requis. Le modèle et ses historiques restent inchangés.
    """
    from scipy.sparse import bsr_matrix,bmat,eye
    from scipy.sparse.linalg import LinearOperator,eigsh,splu
    from scipy.linalg import eigh
    if isinstance(combien,(bool,np.bool_)) or isinstance(maxiter,(bool,np.bool_)):
        raise ValueError('combien et maxiter doivent être des entiers')
    try:combien,maxiter=operator.index(combien),operator.index(maxiter)
    except TypeError:raise ValueError('combien et maxiter doivent être des entiers') from None
    if combien<0 or maxiter<=0:raise ValueError('combien ≥ 0 et maxiter > 0 requis')
    if not np.isfinite(tol) or not 0<tol<1:raise ValueError('0 < tol < 1 requis')
    if decalage is not None and not np.isfinite(decalage):raise ValueError('décalage fini requis')
    if t is not None and not np.isfinite(t):raise ValueError('temps fini requis')
    try:noyau._verifie_domaine_adjoint()
    except ValueError as exc:
        raise ValueError('modes creux hors domaine : '+str(exc).removeprefix('PontNoyau hors domaine : ')) from exc
    k,c,m,g=map(_csc,noyau.k_c_m_g_creux(t))
    n=m.shape[0];nb=n//6
    out=dict(modes=[],ddl_physiques=n,ddl_admissibles=n,rang_contraintes=0,decalage=decalage,
             champ='partie_symetrique_de_la_raideur_locale',tol=tol,
             produits_inverse=0,erreur_arriere_inverse_max=0.,repli_qr_dense=False,base_admissible_dense=False)
    if n==0 or combien==0:return out
    idx=6*np.arange(nb)[:,None]+np.arange(6)
    blocks=np.asarray(m[np.repeat(idx,6,axis=1).ravel(),np.tile(idx,(1,6)).ravel()]).reshape(nb,6,6)
    try:
        l=np.linalg.cholesky(blocks)
        wb=np.linalg.solve(l,np.broadcast_to(np.eye(6),l.shape)).transpose(0,2,1)
    except np.linalg.LinAlgError as exc:raise ValueError('masse physique non définie positive') from exc
    if not np.isfinite(wb).all():raise ValueError('normalisation de masse non représentable')
    w=bsr_matrix((wb,np.arange(nb),np.arange(nb+1)),shape=(n,n)).tocsc();w.eliminate_zeros()
    brut=w.T@k@w
    a=(brut+brut.T)*.5
    norm=_norme(a)
    scale=norm if norm else 1.
    a=(a/scale).tocsc()
    if not np.isfinite(a.data).all():raise ValueError('raideur normalisée non finie')
    constraints=_Contraintes(g@w)
    dimension=n-constraints.rank
    out.update(rang_contraintes=constraints.rank,ddl_admissibles=dimension,
               repli_qr_dense=constraints.repli,defaut_symetrie=_norme(brut-brut.T)/max(_norme(brut),np.finfo(float).tiny),
               norme_raideur_massique=norm,seuil_rang=constraints.seuil)
    if dimension==0:return out
    sigma=-np.sqrt(np.finfo(float).eps) if decalage is None else decalage/scale
    if not np.isfinite(sigma):raise ValueError('décalage normalisé non représentable')
    out['decalage']=sigma*scale
    b=constraints.b;r=constraints.rank
    system=bmat([[a-sigma*eye(n),b.T],[b,None]],format='csc') if r else a-sigma*eye(n)
    systemnorm=_norme(system)
    try:lu=splu(system.tocsc())
    except RuntimeError as exc:raise RuntimeError('inverse décalé singulier : choisir un autre décalage') from exc
    def inverse(v):
        out['produits_inverse']+=1
        pv=constraints.projette(v)
        rhs=np.concatenate((pv,np.zeros((r,)+v.shape[1:])))
        solution=lu.solve(rhs)
        for _ in range(2):
            residual=rhs-system@solution
            if np.max(abs(residual),initial=0.)<=16*np.finfo(float).eps*max(1.,float(np.max(abs(rhs),initial=0.))):break
            solution+=lu.solve(residual)
        residual=rhs-system@solution
        ns=float(np.max(abs(solution),initial=0.));nr=float(np.max(abs(rhs),initial=0.));sc=max(ns,nr)
        backward=(float(np.max(abs(residual),initial=0.))/sc/(systemnorm*(ns/sc)+nr/sc)) if sc else 0.
        if not np.isfinite(backward) or backward>1e-11:
            raise RuntimeError('résolution de l’inverse décalé non vérifiée')
        out['erreur_arriere_inverse_max']=max(out['erreur_arriere_inverse_max'],backward)
        value=constraints.projette(solution[:n])
        if not np.isfinite(value).all():raise RuntimeError('produit inverse non fini')
        return value
    wanted=min(combien,dimension)
    if wanted==n:
        eig,vectors=eigh(a.toarray())
        ids=np.argsort(abs(eig-sigma))[:wanted]
        vectors=vectors[:,ids]
        out['spectre_complet_dense']=True
    else:
        inv=LinearOperator((n,n),matvec=inverse,matmat=inverse,dtype=float)
        rng=np.random.default_rng(371)
        initial=constraints.projette(rng.normal(size=n))
        _,vectors=eigsh(inv,k=wanted,which='LM',tol=min(tol*.01,1e-11),
                        v0=initial,ncv=min(n,max(2*wanted+1,32)),maxiter=maxiter)
        out['spectre_complet_dense']=False
    # Projection, orthonormalisation et Rayleigh–Ritz dans le petit espace
    # calculé : pas de soustraction sigma+1/mu près d'une valeur nulle.
    vectors=constraints.projette(vectors)
    vectors,_=np.linalg.qr(vectors,mode='reduced')
    av=a@vectors
    small=vectors.T@av
    values,rotation=eigh((small+small.T)*.5)
    vectors=vectors@rotation
    ids=np.argsort(abs(values-sigma))
    vectors=vectors[:,ids];values=values[ids]
    orth=float(np.max(abs(vectors.T@vectors-np.eye(wanted))))
    out['defaut_orthogonalite_massique']=orth
    if orth>tol:raise RuntimeError('orthogonalité modale non vérifiée')
    for value,vec in zip(values,vectors.T):
        force=value*vec-a@vec
        reactions=constraints.reactions(force)
        residual=force-constraints.original.T@reactions
        error=float(np.linalg.norm(residual))
        closure=float(np.max(abs(constraints.original@vec),initial=0.))
        backward=error/(1.+abs(value))
        if not np.isfinite(error) or backward>tol or closure>tol:
            raise RuntimeError(f'mode non vérifié : erreur arrière={backward:.3g}, contrainte={closure:.3g}')
        lam=float(value*scale)
        uncertainty=float(error*scale)
        status='positif' if lam>uncertainty else ('negatif' if lam < -uncertainty else 'signe_indetermine')
        shape=w@vec
        # Les réactions réduites n'ont pas une norme minimale en cas de
        # redondance ; leur torseur physique est contrôlé sur toutes les lignes.
        eta=(reactions/constraints.norms)/constraints.scales*scale
        physical=(k@shape+k.T@shape)*.5+g.T@eta-lam*(m@shape)
        original_error=float(np.linalg.norm(w.T@physical)/scale/(1.+abs(value)))
        if not np.isfinite(shape).all() or not np.isfinite(eta).all() or not np.isfinite(original_error) or original_error>tol:
            raise RuntimeError('équation modale physique originale non vérifiée')
        out['modes'].append(dict(valeur_propre=lam,frequence_hz=np.sqrt(lam)/(2*np.pi) if lam>0 else None,
            taux_non_amorti_s=np.sqrt(-lam) if lam<0 else None,statut=status,forme=shape.tolist(),reactions=eta.tolist(),
            erreur_arriere=float(backward),residu_contrainte=closure,residu_lambda_s2=uncertainty,
            erreur_arriere_originale=original_error))
    return out
