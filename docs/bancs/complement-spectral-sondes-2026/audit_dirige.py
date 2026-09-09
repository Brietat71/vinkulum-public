"""Audit indépendant, références Fraction des données binary64 exactes."""
from decimal import Decimal
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, "/tmp/vinkulum-energie-native-travail/ci")
from trace_complement_dirigee import certifier_complement
from test_trace_complement import inverse_contrainte, gram, determinant, difference, restreindre, echelle


def rat(a):
    return [[F(float(v)) for v in row] for row in a]


def exact_k(r, e):
    return gram([[F(float(v))/F(float(e[j])) for j, v in enumerate(row)] for row in r])


def noyau(b):
    a = [list(row) for row in zip(*rat(b))]
    n, pivots, row = b.shape[0], [], 0
    for c in range(n):
        k = next((k for k in range(row, len(a)) if a[k][c]), None)
        if k is None:
            continue
        a[k], a[row] = a[row], a[k]
        pivot = a[row][c]
        a[row] = [v/pivot for v in a[row]]
        for k in range(len(a)):
            if k != row:
                f = a[k][c]
                a[k] = [v-f*w for v,w in zip(a[k], a[row])]
        pivots.append(c)
        row += 1
        if row == len(a):
            break
    libres = [i for i in range(n) if i not in pivots]
    z = [[F(0) for _ in libres] for _ in range(n)]
    for j,c in enumerate(libres):
        z[c][j] = F(1)
        for i,p in enumerate(pivots):
            z[p][j] = -a[i][c]
    return z


def verifier(nom, r, e, m, b, precision, d=None):
    if d is None:
        d = r/e
    k = exact_k(r,e)
    _,tc,tau,chi = inverse_contrainte(k, rat(m), rat(b))
    result = dict(nom=nom,precision=precision,tc_exact=str(tc))
    try:
        c = certifier_complement(d, r, e, m, b, precision=precision)
    except Exception as exc:
        result.update(statut="refus", type=type(exc).__name__, message=str(exc), bilan=getattr(exc,"bilan",None))
        return result
    def contient(c, lo, hi, exact):
        assert F(Decimal(c[lo])) <= exact <= F(Decimal(c[hi])), (nom,lo,exact,c[lo],c[hi])
    contient(c, "trace_complement_inferieure_decimal", "trace_complement_superieure_decimal", tc)
    contient(c, "correction_inferieure_decimal", "correction_superieure_decimal", chi)
    base=c["certificat_total"]
    contient(base, "trace_inferieure_decimal", "trace_superieure_decimal", tau)
    eta=F(Decimal(base["eta_superieur_decimal"]))
    assert F(c["lambda_min"]) <= (1-eta)/tc
    z=noyau(b)
    q=difference(restreindre(gram(rat(d)),z),echelle(restreindre(rat(m),z),F(c["lambda_min"])))
    for j in range(1,len(q)+1):
        assert determinant([row[:j] for row in q[:j]]) >= 0, (nom,j)
    result.update(statut="valide_exactement",lambda_min=c["lambda_min"],
        facteur_perte_trace=float(F(Decimal(c["trace_complement_superieure_decimal"]))/tc),
        trace_borne_inferieure_negative=Decimal(c["trace_complement_inferieure_decimal"])<0,
        operations=c["operations_decimal"])
    return result


def main():
    reports=[]
    r0=np.array([[2.,1.,0.,1.],[0.,3.,1.,0.],[0.,0.,2.,1.],[0.,0.,0.,4.]])
    e0=np.array([.5,2.,.25,4.])
    l=np.array([[1.,2.,0.,0.],[0.,2.,1.,0.],[0.,0.,3.,1.],[0.,0.,0.,2.]])
    m=l.T@l
    b=np.array([[1.,0.],[0.,1.],[2.,-1.],[1.,3.]])
    for p in (16,40,80):
        for exposant in (50,200,400):
            be=b*np.array([2.**-exposant,2.**exposant])
            reports.append(verifier(f"colonnes_B_2_pm{exposant}",r0,e0,m,be,p))
        for puissance in (20,40,52):
            be=np.array([[1.,1.],[0.,2.**-puissance],[0.,0.],[0.,0.]])
            reports.append(verifier(f"B_quasi_redondant_2_m{puissance}",r0,e0,m,be,p))
        for exposant in (20,50,100):
            rr=np.diag([2.**-exposant,1.,2.**exposant,2.**(2*exposant)])
            reports.append(verifier(f"R_anisotrope_2_pm{exposant}",rr,np.ones(4),m,b,p))
        rr=np.diag([2.**-50,1.,2.,3.])
        reports.append(verifier("cancellation_rang1",rr,np.ones(4),np.eye(4),np.eye(4)[:,:1],p))
    rng=np.random.default_rng(903)
    for i in range(40):
        rr=np.triu(rng.integers(-4,5,(4,4))).astype(float)
        np.fill_diagonal(rr,rng.integers(1,5,4))
        ee=2.**rng.integers(-8,9,4)
        ll=np.triu(rng.integers(-4,5,(4,4))).astype(float)
        np.fill_diagonal(ll,rng.integers(1,5,4))
        mm=ll.T@ll
        bb=rng.integers(-4,5,(4,2)).astype(float)*2.**rng.integers(-30,31,2)
        dd=np.vstack((rr/ee, 2.**-30*rng.integers(-2,3,4)))
        reports.append(verifier(f"aleatoire_{i}",rr,ee,mm,bb,16 if i%2 else 80,d=dd))
    src=Path('/tmp/vinkulum-energie-native-travail/ci/trace_complement_dirigee.py')
    report=dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),cas=reports,
        acceptes=sum(c['statut']=='valide_exactement' for c in reports),
        refuses=sum(c['statut']=='refus' for c in reports))
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    print(report['acceptes'], 'certificats vérifiés exactement;',report['refuses'],'refus')
    for c in reports:
        if c['statut']=='refus':
            print(c['nom'],c['precision'],c['type'],c['message'])


if __name__=='__main__':
    main()
