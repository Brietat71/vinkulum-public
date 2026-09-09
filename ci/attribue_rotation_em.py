"""Attribution exacte encadrée des écarts de rotation et vitesse spatiale.

La différence de Cayley emploie une identité de résolvante exacte. Les
transformations de repères stockées ne sont jamais supposées orthogonales.
Le post-traitement flottant des observations est séparé explicitement.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path
import time

import numpy as np
from algebre_moments_em import sub,eye,mv,mm,ma,ms,skew,inv
from defauts_moments_em import matrice,encoder,cayley,mnorm,vmax,det,ratio,BUDGET
from attribue_moments_em import (attribuer,proposer,point,imv,iadd,ipoint_add,
                                inorme,largeur,BITS)

BUDGET_R=F(1,10**20)
BUDGET_W=F(1,10**18)


def transpose(a):return [list(x) for x in zip(*a)]
def flat(a):return [x for row in a for x in row]
def reshape(a):return [a[3*i:3*i+3] for i in range(3)]
def izmat():return [[point(F(0)) for _ in range(3)] for _ in range(3)]
def ipmat(a):return [[point(x) for x in row] for row in a]
def imright(a,b):return [imv(transpose(b),row) for row in a]
def iapoint(a,b):return [ipoint_add(row,const) for row,const in zip(a,b)]
def iamat(a,b):return [iadd(x,y) for x,y in zip(a,b)]
def imatvec(a,v):return [imv([v],row)[0] for row in a]


def sqrt_bounds(x):
    q=isqrt((x.numerator<<(2*BITS))//x.denominator)
    lo=F(q,1<<BITS)
    return lo,lo if lo*lo==x else F(q+1,1<<BITS)


def norm_bounds(x):
    low=sum((F(0) if l<=0<=h else min(l*l,h*h)) for l,h in x)
    high=sum(max(l*l,h*h) for l,h in x)
    return sqrt_bounds(low)[0],sqrt_bounds(high)[1]


class RotationVitesse:
    def canon(self,row_a,row_b):
        a=matrice(row_a[2][0]);native_b=matrice(row_b[2][0])
        b=mm(mm(self.mt,native_b),self.ct)
        wa=list(map(F,row_a[3][0]));wb=mv(self.mt,list(map(F,row_b[3][0])))
        return a,b,wa,wb

    def observe(self,index,rcomp,wcomp,ra,rb,wa,wb):
        delta_r=ma(rb,ms(-1,ra));delta_w=sub(wb,wa)
        for rotation in (ra,matrice(self.rows_b[index][2][0])):
            self.orthogonality=max(self.orthogonality,mnorm(ma(mm(transpose(rotation),rotation),ms(-1,eye()))))
            self.determinants_positifs=self.determinants_positifs and det(rotation)>0
        ir=iamat(iamat(rcomp[0],rcomp[1]),rcomp[2])
        iw=iadd(iadd(wcomp[0],wcomp[1]),wcomp[2])
        if not all(l<=v<=h for (l,h),v in zip(flat(ir),flat(delta_r))):
            raise ValueError('rotation hors enclosure')
        if not all(l<=v<=h for (l,h),v in zip(iw,delta_w)):
            raise ValueError('vitesse hors enclosure')
        self.width_r=max(self.width_r,largeur(flat(ir)))
        self.width_w=max(self.width_w,largeur(iw))
        # Valeurs effectivement observées après le post-traitement historique.
        # NumPy propose ces valeurs ; leurs erreurs sont ensuite des différences
        # de rationnels exacts, sans supposer quoi que ce soit sur son arrondi.
        xa,xb=self.rows_a[index],self.rows_b[index]
        nr_a=np.asarray(xa[2][0]).reshape(3,3)
        nr_b=self.world.T@np.asarray(xb[2][0]).reshape(3,3)@self.material.T
        nw_a=np.asarray(xa[3][0]);nw_b=self.world.T@np.asarray(xb[3][0])
        reported_r=matrice((nr_b-nr_a).ravel().tolist())
        reported_w=list(map(F,nw_b-nw_a))
        post_r=ma(reported_r,ms(-1,delta_r));post_w=sub(reported_w,delta_w)
        comps_r=rcomp+[ipmat(post_r)];comps_w=wcomp+[[point(x) for x in post_w]]
        sr=iamat(ir,ipmat(post_r));sw=ipoint_add(iw,post_w)
        if not all(l<=v<=h for (l,h),v in zip(flat(sr),flat(reported_r))):
            raise ValueError('rotation publiée hors enclosure')
        if not all(l<=v<=h for (l,h),v in zip(sw,reported_w)):
            raise ValueError('vitesse publiée hors enclosure')
        for key,actual,comps in (('rotation',flat(reported_r),[flat(x) for x in comps_r]),
                                  ('omega',reported_w,comps_w)):
            square=sum(x*x for x in actual)
            self.maxima[key]=[max(old,norm_bounds(x)[1]) for old,x in zip(self.maxima[key],comps)]
            if square>self.peaks[key]['carre']:
                self.peaks[key]=dict(carre=square,pas=index,norme=sqrt_bounds(square),contributions=[norm_bounds(x) for x in comps],vecteurs=comps)
            if index>0 and index%self.stride==0 and square>self.sample_peaks[key]['carre']:
                self.sample_peaks[key]=dict(carre=square,pas=index,norme=sqrt_bounds(square),contributions=[norm_bounds(x) for x in comps],vecteurs=comps)
        self.last_r,self.last_w=comps_r,comps_w
        self.count=index

    def omega_components(self,index,rcomponents,pcomponents,ra,rb,wa,wb):
        pa=list(map(F,self.rows_a[index][1][0]));pb=mv(self.c,list(map(F,self.rows_b[index][1][0])))
        base=mv(self.ai,pa);coef=mm(rb,self.ai)
        result=[iadd(imatvec(r,base),imv(coef,p)) for r,p in zip(rcomponents,pcomponents)]
        model=mv(mm(rb,ma(self.hb,ms(-1,self.ai))),pb)
        native_rb=matrice(self.rows_b[index][2][0])
        native_pb=list(map(F,self.rows_b[index][1][0]))
        native_wb=list(map(F,self.rows_b[index][3][0]))
        defect_a=sub(wa,mv(ra,base))
        defect_b=sub(native_wb,mv(native_rb,mv(self.bi,native_pb)))
        self.local_omega=max(self.local_omega,ratio(vmax(defect_a),mnorm(ra)*mnorm(self.ai)*vmax(pa)),
                             ratio(vmax(defect_b),mnorm(native_rb)*mnorm(self.bi)*vmax(native_pb)))
        delta_defect=sub(mv(self.mt,defect_b),defect_a)
        result[1]=ipoint_add(result[1],model)
        result[2]=ipoint_add(result[2],delta_defect)
        return result

    def debut(self,da,db,ai,bi,c,pcomponents):
        self.rows_a=da['trace']['echantillons'];self.rows_b=db['trace']['echantillons']
        self.ai,self.bi,self.c=ai,bi,c
        self.ct=transpose(c);self.ci=inv(c);self.cit=transpose(self.ci)
        self.mt=transpose([[F(x) for x in row] for row in db['monde']])
        self.world=np.asarray(db['monde']);self.material=np.asarray(db['matiere'])
        self.hb=mm(mm(self.cit,bi),self.ci)
        self.unit=[skew(mv(ai,[F(i==j) for i in range(3)])) for j in range(3)]
        self.q=eye();self.width_r=self.width_w=self.base_error=F(0)
        self.local_rotation=self.local_omega=self.orthogonality=F(0)
        self.determinants_positifs=True
        self.frame_error=max(mnorm(ma(mm(transpose(c),c),ms(-1,eye()))),
                             mnorm(ma(mm(self.mt,transpose(self.mt)),ms(-1,eye()))))
        self.frames_positifs=det(c)>0 and det(self.mt)>0
        self.stride=max(1,round(db['duree']/db['pas_demande'])//200)
        self.peaks={k:dict(carre=F(-1)) for k in ('rotation','omega')}
        self.sample_peaks={k:dict(carre=F(-1)) for k in ('rotation','omega')}
        self.maxima={k:[F(0)]*4 for k in ('rotation','omega')}
        ra,rb,wa,wb=self.canon(self.rows_a[0],self.rows_b[0])
        self.za=[ipmat(ma(rb,ms(-1,ra))),izmat(),izmat()]
        wcomp=self.omega_components(0,self.za,pcomponents,ra,rb,wa,wb)
        self.observe(0,self.za,wcomp,ra,rb,wa,wb)
        self.ra,self.rb=ra,rb
        self.native_rb=matrice(self.rows_b[0][2][0])

    def pas(self,index,h,mx,mp,my,previous_components,pcomponents):
        ca=cayley([h*x for x in mv(self.ai,mx)])
        native_cb=cayley([h*x for x in mv(self.bi,mp)])
        cb=mm(mm(self.cit,native_cb),self.ct)
        ua,ub=ms(F(1,2),ma(ca,eye())),ms(F(1,2),ma(cb,eye()))
        xb=mm(mm(self.cit,skew([h*x for x in mv(self.bi,mp)])),self.ct)
        xa=skew([h*x for x in mv(self.ai,mx)])
        if ma(cb,ms(-1,ca))!=mm(mm(ub,ma(xb,ms(-1,xa))),ua):
            raise ValueError('identité de résolvante violée')
        qnext=proposer(self.q,ca);qinv=inv(qnext)
        normalized=mm(mm(self.q,ca),qinv)
        self.base_error=max(self.base_error,mnorm(ma(normalized,ms(-1,eye()))))
        left,right=mm(self.rb,ub),mm(ua,qinv)
        maps=[ms(h,mm(mm(left,u),right)) for u in self.unit]
        coeff=[[maps[j][i][k] for j in range(3)] for i in range(3) for k in range(3)]
        model_x=ma(xb,ms(-h,skew(mv(self.ai,my))))
        model=mm(mm(left,model_x),right)
        ra,rb,wa,wb=self.canon(self.rows_a[index],self.rows_b[index])
        native_rb=matrice(self.rows_b[index][2][0])
        defect_a=ma(ra,ms(-1,mm(self.ra,ca)))
        defect_b=ma(native_rb,ms(-1,mm(self.native_rb,native_cb)))
        self.local_rotation=max(self.local_rotation,
            ratio(vmax(flat(defect_a)),max(mnorm(ra),mnorm(self.ra))),
            ratio(vmax(flat(defect_b)),max(mnorm(native_rb),mnorm(self.native_rb))))
        defect=ma(mm(mm(self.mt,defect_b),self.ct),ms(-1,defect_a))
        direct=mm(defect,qinv)
        next_z=[]
        for i,(prev,cur) in enumerate(zip(previous_components,pcomponents)):
            mid=[((a+c)/2,(b+d)/2) for (a,b),(c,d) in zip(prev,cur)]
            forcing=reshape(imv(coeff,mid))
            z=iamat(imright(self.za[i],normalized),forcing)
            if i==1:z=iapoint(z,model)
            if i==2:z=iapoint(z,direct)
            next_z.append(z)
        rcomp=[imright(z,qnext) for z in next_z]
        wcomp=self.omega_components(index,rcomp,pcomponents,ra,rb,wa,wb)
        self.observe(index,rcomp,wcomp,ra,rb,wa,wb)
        self.q,self.za,self.ra,self.rb=qnext,next_z,ra,rb
        self.native_rb=native_rb

    def bilan(self):
        return dict(pas_verifies=self.count,largeur_rotation_max=self.width_r,largeur_omega_max=self.width_w,
                    budget_rotation=BUDGET_R,budget_omega=BUDGET_W,erreur_base_rotation_max=self.base_error,
                    maxima_contributions=self.maxima,pics_tous_pas=self.peaks,pics_echantillonnage_historique=self.sample_peaks,
                    contributions_finales_rotation=self.last_r,contributions_finales_omega=self.last_w,
                    noms_contributions=['initial','donnees_transformees','defauts_algebriques','post_traitement'],
                    defaut_rotation_local_max=self.local_rotation,defaut_omega_local_max=self.local_omega,
                    orthogonalite_native_max=self.orthogonality,orthogonalite_reperes_max=self.frame_error,
                    determinants_positifs=self.determinants_positifs and self.frames_positifs,
                    budget_local=BUDGET,borne_erreur_ode_continue=False,
                    programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),numpy=np.__version__,
                    verifie=(self.width_r<=BUDGET_R and self.width_w<=BUDGET_W
                             and max(self.local_rotation,self.local_omega,self.orthogonality,self.frame_error)<=BUDGET
                             and self.determinants_positifs and self.frames_positifs))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('initial');p.add_argument('autre');p.add_argument('--sortie',required=True);p.add_argument('--limite',type=int)
    args=p.parse_args();start=time.perf_counter()
    d=attribuer(args.initial,args.autre,args.limite,RotationVitesse())
    Path(args.sortie).write_text(json.dumps(encoder(d),indent=2)+'\n')
    print('terminé',time.perf_counter()-start,d['attribution_rotation_vitesse_verifiee'],flush=True)
