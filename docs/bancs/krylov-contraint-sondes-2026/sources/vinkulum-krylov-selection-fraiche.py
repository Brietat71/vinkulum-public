import sys,json,time
from pathlib import Path
sys.path.insert(0,'ci')
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
from vinkulum._ports.condensation_energie import CondensationEnergie
from confronte_ports_exudyn import lire,juger
from trace_complement_dirigee import certifier_complement
from krylov_contraint import KrylovContraint
from controle_complement import ControleComplement
rows=[]
for n in (32,128,512):
 p=Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee')/f'n{n}-f40.npz'
 d,m,metric,f,omega=lire(p);ni=d.shape[1]-6
 t=time.perf_counter();qr=CondensationEnergie(d,np.arange(ni),np.arange(ni,d.shape[1]),metric);tqr=time.perf_counter()-t
 mi=m[:ni,:ni].tocsr();rac=np.sqrt(mi.diagonal())
 def ap(x):
  v=x.ndim==1;y=x[:,None] if v else x
  z=rac[:,None]*qr.solve(rac[:,None]*y)
  return z[:,0] if v else z
 t=time.perf_counter();val,z=eigsh(LinearOperator((ni,ni),matvec=ap,matmat=ap,dtype=float),k=1,which='LM',tol=1e-11,ncv=20,v0=np.random.default_rng(107+n).standard_normal(ni));tm=time.perf_counter()-t
 phi=z/rac[:,None];b=mi@phi
 t=time.perf_counter();cert=certifier_complement(qr.di,qr.r,qr.echelles,mi,b);tc=time.perf_counter()-t
 k=KrylovContraint(qr,m,b,phi,cert['lambda_min'],omega[-1],blocs=4)
 c=ControleComplement(k)
 t=time.perf_counter();reps=[c.reponses(w,f) for w in omega];tr=time.perf_counter()-t
 j=juger(d,m,metric,np.array([a['champ'] for a in reps]),np.load(p.with_suffix('.ref.npy'),mmap_mode='r'))
 bs={nom:max(x for a in reps for x in a['bornes'][nom]['relatives']) for nom in ('masse','deformation')}
 row=dict(n=n,qr=tqr,mode=tm,cert=tc,krylov=k.preparation_s,controle=c.preparation_s,reponses=tr,taille=k.taille_reduite,bornes=bs,juge=j)
 rows.append(row);Path('/tmp/vinkulum-krylov-selection-fraiche.json').write_text(json.dumps(rows,allow_nan=False))
 print(n,k.taille_reduite,bs,j['accepte'],j['maxima_operateurs'],[tqr,tm,tc,k.preparation_s,c.preparation_s,tr],flush=True)
