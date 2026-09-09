"""Défauts exacts du point milieu et de Cayley sur une trace native.

Les binary64 lus deviennent des rationnels exacts. Le budget local est fixé
à64epsilon avant lecture des défauts. Un contrôle local n'est pas une borne
de propagation globale ni une attribution causale complète des écarts.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
import math
from pathlib import Path
import time

import algebre_moments_em as alg
from algebre_moments_em import add,sub,scale,eye,mv,mm,ma,ms,cross,skew,inv

BUDGET=F(64,2**52)


def matrice(x):return [[F(x[3*i+j]) for j in range(3)] for i in range(3)]
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def vmax(x):return max(map(abs,x),default=F(0))
def mnorm(x):return max(sum(map(abs,row)) for row in x)
def det(x):return dot(x[0],cross(x[1],x[2]))


def cayley(v):
    q=dot(v,v);s=skew(v);den=1+q/4
    return [[F(i==j)+(s[i][j]+(v[i]*v[j]-q*F(i==j))/2)/den for j in range(3)] for i in range(3)]


def ratio(x,s):
    if s:return x/s
    if x:raise ValueError('échelle nulle pour un défaut non nul')
    return F(0)


def encoder(x):
    if isinstance(x,F):return {'n':str(x.numerator),'d':str(x.denominator)}
    if isinstance(x,dict):return {k:encoder(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encoder(v) for v in x]
    return x


def valider_trace(data):
    """Refuse les traces tronquées ou ambiguës avant tout calcul de preuve."""
    def number(x):
        if type(x) not in (int,float) or not math.isfinite(x) or F(float(x))!=F(x):
            raise ValueError('valeur binary64 finie requise')
    def vector(x,n):
        if not isinstance(x,list) or len(x)!=n:raise ValueError('dimension de trace incorrecte')
        for v in x:number(v)
    trace=data['trace']
    if trace['schema']!='vinkulum.trace.moments_em.1' or len(trace['inerties'])!=1:
        raise ValueError('trace mon corps attendue')
    vector(trace['inerties'][0],9)
    j=matrice(trace['inerties'][0])
    if any(j[i][k]!=j[k][i] for i in range(3) for k in range(3)):
        raise ValueError('inertie non symétrique')
    if min(j[0][0],j[0][0]*j[1][1]-j[0][1]**2,det(j))<=0:
        raise ValueError('inertie non définie positive')
    for key in ('monde','matiere'):
        if len(data[key])!=3:raise ValueError('repère de dimension incorrecte')
        for row in data[key]:vector(row,3)
    number(data['duree']);number(data['pas_demande'])
    if data['duree']<0 or data['pas_demande']<=0:raise ValueError('durée ou pas invalide')
    rows=trace['echantillons']
    if not rows:raise ValueError('trace vide')
    previous=None
    for row in rows:
        if len(row)!=4:raise ValueError('échantillon de dimension incorrecte')
        number(row[0]);t=F(row[0])
        for block,n in zip(row[1:],(3,9,3),strict=True):
            if len(block)!=1:raise ValueError('trace mon corps attendue')
            vector(block[0],n)
        if previous is not None:
            h=t-previous
            # Contrat natif numerique::fin_pas : le dernier pas peut
            # absorber le reliquat jusqu'à 1,5 h. Les équations auditées
            # utilisent ensuite la différence rationnelle des dates stockées.
            expected=min(float(previous)+data['pas_demande'],data['duree'])
            if data['duree']-expected<=.5*data['pas_demande']:
                expected=data['duree']
            if h<=0 or t!=F(expected):
                raise ValueError('grille non croissante ou pas manquant')
        previous=t
    if rows[0][0]!=0 or F(rows[-1][0])!=F(data['duree']):
        raise ValueError('horizon incomplet')


def verifier(source,limite=None):
    source=Path(source);payload=source.read_bytes()
    data=json.loads(gzip.decompress(payload));trace=data['trace']
    valider_trace(data)
    if limite is not None and (type(limite)!=int or limite<0):raise ValueError('limite invalide')
    if trace['schema']!='vinkulum.trace.moments_em.1' or len(trace['inerties'])!=1:
        raise ValueError('trace mon corps attendue')
    j=matrice(trace['inerties'][0])
    if any(j[i][k]!=j[k][i] for i in range(3) for k in range(3)):
        raise ValueError('inertie non symétrique')
    if min(j[0][0],j[0][0]*j[1][1]-j[0][1]**2,det(j))<=0:
        raise ValueError('inertie non définie positive')
    ji=inv(j);jin=mnorm(ji)
    rows=trace['echantillons']
    total=len(rows)-1
    if limite is not None:rows=rows[:limite+1]
    if not rows:raise ValueError('trace vide')
    t0=F(rows[0][0]);p0=list(map(F,rows[0][1][0]));r0=matrice(rows[0][2][0])
    energy0=dot(p0,mv(ji,p0))/2;moment20=dot(p0,p0)
    energy=energy0;moment2=moment20
    sums=[F(0),F(0)]
    maxima={k:dict(valeur=F(0),pas=0) for k in ('moment','rotation','omega','energie_relative','moment_carre_relatif')}
    for index,row in enumerate(rows[1:],1):
        t1=F(row[0]);h=t1-t0
        if h<=0:raise ValueError('temps non croissant')
        p1=list(map(F,row[1][0]));r1=matrice(row[2][0]);w1=list(map(F,row[3][0]))
        pm=scale(F(1,2),add(p0,p1));wm=mv(ji,pm)
        defect=sub(sub(p1,p0),scale(h,cross(pm,wm)))
        cm=cayley(scale(h,wm))
        rotated=mm(r0,cm)
        rd=[sub(a,b) for a,b in zip(r1,rotated)]
        wd=sub(w1,mv(r1,mv(ji,p1)))
        energy1=dot(p1,mv(ji,p1))/2;moment21=dot(p1,p1)
        de=dot(wm,defect);dl=2*dot(pm,defect)
        if energy1-energy!=de or moment21-moment2!=dl:
            raise ValueError('identités quadratiques non satisfaites')
        sums[0]+=de;sums[1]+=dl
        values=dict(moment=ratio(vmax(defect),max(vmax(p0),vmax(p1))),
                    rotation=ratio(max(vmax(row) for row in rd),max(mnorm(r0),mnorm(r1))),
                    omega=ratio(vmax(wd),mnorm(r1)*jin*vmax(p1)),
                    energie_relative=ratio(abs(energy1-energy0),energy0),
                    moment_carre_relatif=ratio(abs(moment21-moment20),moment20))
        for key,value in values.items():
            if value>maxima[key]['valeur']:maxima[key]=dict(valeur=value,pas=index)
        t0,p0,r0,energy,moment2=t1,p1,r1,energy1,moment21
        if index%2000==0:print(index,{k:float(v['valeur']) for k,v in maxima.items()},flush=True)
    if energy-energy0!=sums[0] or moment2-moment20!=sums[1]:
        raise ValueError('sommation des défauts incorrecte')
    accepted=all(maxima[k]['valeur']<=BUDGET for k in ('moment','rotation','omega'))
    return dict(schema='vinkulum.diagnostic.defauts_em.1',source_sha256=hashlib.sha256(payload).hexdigest(),
                programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                algebre_sha256=hashlib.sha256(Path(alg.__file__).read_bytes()).hexdigest(),
                cadre=data['cadre'],pas_demande=data['pas_demande'],pas_verifies=len(rows)-1,pas_total=total,
                complet=len(rows)-1==total,budget_local=BUDGET,budget_local_respecte=accepted,
                maxima=maxima,derive_energie=energy-energy0,derive_moment_carre=moment2-moment20,
                sommes_defauts=sums,propagation_globale_certifiee=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source');p.add_argument('--sortie',required=True);p.add_argument('--limite',type=int)
    a=p.parse_args();start=time.perf_counter();d=verifier(a.source,a.limite)
    Path(a.sortie).write_text(json.dumps(encoder(d),indent=2)+'\n')
    print('terminé',time.perf_counter()-start,d['pas_verifies'],d['complet'],d['budget_local_respecte'],flush=True)
