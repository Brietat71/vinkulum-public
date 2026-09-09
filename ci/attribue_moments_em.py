"""Attribution validée des écarts de moments entre deux traces.

Champ quadratique : identité de différence exacte, sans linéarisation
approchée. Une base fondamentale proposée en Decimal est vérifiée en
rationnels à chaque pas. Les enclosures sont arrondies vers l'extérieur.
Ce module ne couvre pas encore la rotation ni la vitesse spatiale.
"""
import argparse
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import time

from algebre_moments_em import add,sub,scale,eye,mv,mm,ma,ms,cross,inv
from defauts_moments_em import matrice,dot,vmax,mnorm,encoder,BUDGET,ratio,valider_trace

BITS=160
LARGEUR_MAX=F(1,10**20)


def plancher(x):return F((x.numerator<<BITS)//x.denominator,1<<BITS)
def plafond(x):return -plancher(-x)
def point(x):return (x,x)
def izero():return [point(F(0)) for _ in range(3)]


def imv(a,x):
    out=[]
    for row in a:
        low=high=F(0)
        for c,(l,h) in zip(row,x):
            low+=c*(l if c>=0 else h)
            high+=c*(h if c>=0 else l)
        out.append((plancher(low),plafond(high)))
    return out


def iadd(x,y):return [(plancher(a+c),plafond(b+d)) for (a,b),(c,d) in zip(x,y)]
def ipoint_add(x,y):return iadd(x,[point(v) for v in y])
def inorme(x):return max(max(abs(l),abs(h)) for l,h in x)
def largeur(x):return max(h-l for l,h in x)


def proposer(t,p):
    # Propositions seulement : ni Decimal ni son arrondi ne sont utilisés
    # comme preuve. Les matrices dyadiques résultantes sont vérifiées ensuite.
    with localcontext() as ctx:
        ctx.prec=70
        td=[[D(v.numerator)/D(v.denominator) for v in row] for row in t]
        pd=[[D(v.numerator)/D(v.denominator) for v in row] for row in p]
        return [[plancher(F(sum(td[i][k]*pd[k][j] for k in range(3)))) for j in range(3)] for i in range(3)]


def lire(path):
    raw=Path(path).read_bytes()
    return json.loads(gzip.decompress(raw)),hashlib.sha256(raw).hexdigest()


def attribuer(source_a,source_b,limite=None,observateur=None):
    da,ha=lire(source_a);db,hb=lire(source_b)
    valider_trace(da);valider_trace(db)
    if limite is not None and (type(limite)!=int or limite<0):raise ValueError('limite invalide')
    if any(da[k]!=db[k] for k in ('duree','pas_demande')):raise ValueError('horizons différents')
    if any(da[k]!=eye() for k in ('monde','matiere')):raise ValueError('repère initial non canonique')
    if da['cadre']!='initial':raise ValueError('repère canonique initial requis')
    a,b=da['trace'],db['trace']
    if any(t['schema']!='vinkulum.trace.moments_em.1' or len(t['inerties'])!=1 for t in (a,b)):
        raise ValueError('deux traces mon corps attendues')
    ja,jb=matrice(a['inerties'][0]),matrice(b['inerties'][0])
    if any(ja[i][j] for i in range(3) for j in range(3) if i!=j):
        raise ValueError('inertie canonique diagonale requise')
    ai,bi=inv(ja),inv(jb)
    c=[[F(v) for v in row] for row in db['matiere']]
    if mm(c,inv(c))!=eye():raise ValueError('transformation non inversible')
    fa=lambda x:cross(x,mv(ai,x))
    fb=lambda x:cross(x,mv(bi,x))
    alpha=[ai[2][2]-ai[1][1],ai[0][0]-ai[2][2],ai[1][1]-ai[0][0]]
    def jac(x):
        return [[F(0),alpha[0]*x[2],alpha[0]*x[1]],
                [alpha[1]*x[2],F(0),alpha[1]*x[0]],
                [alpha[2]*x[1],alpha[2]*x[0],F(0)]]
    ra,rb=a['echantillons'],b['echantillons']
    if len(ra)!=len(rb):raise ValueError('grilles de tailles différentes')
    total=len(ra)-1
    if limite is not None:ra,rb=ra[:limite+1],rb[:limite+1]
    if ra[0][0]!=rb[0][0]:raise ValueError('temps initiaux différents')
    t0=F(ra[0][0]);x0=list(map(F,ra[0][1][0]));p0=list(map(F,rb[0][1][0]));y0=mv(c,p0)
    initial=sub(y0,x0)
    coords=[[point(v) for v in initial],izero(),izero()]
    components=[coords[0],izero(),izero()]
    basis=eye();peaks=[inorme(v) for v in components]
    if observateur is not None:observateur.debut(da,db,ai,bi,c,components)
    observed_peak=vmax(initial);peak_step=0;peak_components=components
    width=F(0);base_error=F(0);residuals=[F(0),F(0)]
    for index,(rowa,rowb) in enumerate(zip(ra[1:],rb[1:],strict=True),1):
        if rowa[0]!=rowb[0]:raise ValueError('temps différents')
        t1=F(rowa[0]);h=t1-t0
        if h<=0:raise ValueError('temps non croissant')
        x1=list(map(F,rowa[1][0]));p1=list(map(F,rowb[1][0]));y1=mv(c,p1)
        mx=scale(F(1,2),add(x0,x1));mp=scale(F(1,2),add(p0,p1));my=mv(c,mp)
        defect_a=sub(sub(x1,x0),scale(h,fa(mx)))
        defect_b=sub(sub(p1,p0),scale(h,fb(mp)))
        residuals[0]=max(residuals[0],ratio(vmax(defect_a),max(vmax(x0),vmax(x1))))
        residuals[1]=max(residuals[1],ratio(vmax(defect_b),max(vmax(p0),vmax(p1))))
        j=jac(scale(F(1,2),add(mx,my)))
        left=ma(eye(),ms(-h/2,j));right=ma(eye(),ms(h/2,j))
        model=scale(h,sub(mv(c,fb(mp)),fa(my)))
        rounding=sub(mv(c,defect_b),defect_a)
        delta0,delta1=sub(y0,x0),sub(y1,x1)
        if mv(left,delta1)!=add(mv(right,delta0),add(model,rounding)):
            raise ValueError('identité discrète violée')
        li=inv(left);transition=mm(li,right)
        next_basis=proposer(transition,basis)
        ni=inv(next_basis)
        normalized=mm(mm(ni,transition),basis)
        err=ma(normalized,ms(-1,eye()))
        base_error=max(base_error,mnorm(err))
        forcing=mm(ni,li)
        sources=[[F(0)]*3,mv(forcing,model),mv(forcing,rounding)]
        coords=[ipoint_add(imv(normalized,z),f) for z,f in zip(coords,sources)]
        previous_components=components
        components=[imv(next_basis,z) for z in coords]
        summed=iadd(iadd(components[0],components[1]),components[2])
        if not all(l<=v<=u for (l,u),v in zip(summed,delta1)):
            raise ValueError('recomposition hors enclosure')
        if observateur is not None:
            observateur.pas(index,h,mx,mp,my,previous_components,components)
        width=max(width,largeur(summed))
        for i,component in enumerate(components):peaks[i]=max(peaks[i],inorme(component))
        observed=vmax(delta1)
        if observed>observed_peak:observed_peak=observed;peak_step=index;peak_components=components
        basis=next_basis;t0,x0,p0,y0=t1,x1,p1,y1
        if index%1000==0:
            print(index,'ecart',float(observed_peak),'largeur',float(width),'base',float(base_error),flush=True)
    result=dict(schema='vinkulum.diagnostic.attribution_moment_em.1',
                sources_sha256=[ha,hb],programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                cadre=db['cadre'],pas_demande=db['pas_demande'],pas_verifies=len(ra)-1,pas_total=total,
                complet=len(ra)-1==total,budget_local=BUDGET,defauts_locaux_max=residuals,
                precision_bits=BITS,budget_largeur=LARGEUR_MAX,largeur_max=width,
                erreur_base_max=base_error,ecart_moment_max=observed_peak,pas_ecart_max=peak_step,
                noms_contributions=['initial','donnees_transformees','defauts_algebriques'],
                maxima_contributions=peaks,contributions_au_pic=peak_components,contributions_finales=components,
                attribution_moment_verifiee=len(ra)-1==total and width<=LARGEUR_MAX and max(residuals)<=BUDGET,
                attribution_rotation_vitesse_verifiee=False)
    if observateur is not None:
        result['rotation_vitesse']=observateur.bilan()
        result['attribution_rotation_vitesse_verifiee']=(result['attribution_moment_verifiee'] and result['rotation_vitesse']['verifie'])
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('initial');p.add_argument('autre');p.add_argument('--sortie',required=True);p.add_argument('--limite',type=int)
    args=p.parse_args();start=time.perf_counter();d=attribuer(args.initial,args.autre,args.limite)
    Path(args.sortie).write_text(json.dumps(encoder(d),indent=2)+'\n')
    print('terminé',time.perf_counter()-start,d['attribution_moment_verifiee'],float(d['largeur_max']),flush=True)
