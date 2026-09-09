import os,sys,json,time,hashlib,shutil
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
from fractions import Fraction
os.sched_setaffinity(0,{8})
out=Path('/tmp/vinkulum-inertie-controle-pilote-v1')
out.mkdir(exist_ok=False)
ci=Path('/tmp/vinkulum-energie-native-travail/ci')
(out/'sources').mkdir()
sources=('inertie_complement_dirigee.py','krylov_contraint.py','controle_complement.py','inverse_contrainte_energie.py','racine_masse_blocs.py','condensation_energie.py','confronte_ports_exudyn.py')
for f in sources: shutil.copy2(ci/f,out/'sources'/f)
shutil.copy2(__file__,out/'pilote.py')
sys.path.insert(0,str(out/'sources'))
from confronte_ports_exudyn import lire,juger
from condensation_energie import CondensationEnergie
from inertie_complement_dirigee import certifier_inertie_complement
from krylov_contraint import KrylovContraint
from controle_complement import ControleComplement
resultats=[]
for n in (32,128,512):
 d,m,metrique,forces,omega=lire(Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee')/f'n{n}-f40.npz')
 t=time.perf_counter();qr=CondensationEnergie(d,np.arange(6*n-6),np.arange(6*n-6,6*n),metrique);tqr=time.perf_counter()-t
 mi=m[qr.i][:,qr.i].tocsr();racine=np.sqrt(mi.diagonal())
 def appliquer(v):return (racine[:,None]*qr.solve(racine[:,None]*np.asarray(v).reshape(-1,1))).ravel()
 t=time.perf_counter();_,vectors=eigsh(LinearOperator(mi.shape,matvec=appliquer,dtype=float),k=1,which='LA',tol=1e-11,ncv=20,v0=np.random.default_rng(107+n).normal(size=len(qr.i)));phi=vectors/racine[:,None];b=mi@phi;tmode=time.perf_counter()-t
 t=time.perf_counter();cert=certifier_inertie_complement(d[:,qr.i],mi,b,(2*np.pi*80)**2,precision=32);tc=time.perf_counter()-t
 assert Fraction(cert['lambda_min'])>Fraction(float(omega[-1]))**2
 np.savez_compressed(out/f'n{n}-directions.npz',B=b,Phi=phi)
 (out/f'n{n}-certificat.json').write_text(json.dumps(cert,ensure_ascii=False,indent=2))
 ref=np.load(Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee')/f'n{n}-f40.ref.npy')
 for blocs in (3,4):
  t=time.perf_counter();r=KrylovContraint(qr,m,b,phi,cert['lambda_min'],omega[-1],blocs=blocs);tk=time.perf_counter()-t
  t=time.perf_counter();c=ControleComplement(r,profondeur=8);tctrl=time.perf_counter()-t
  t=time.perf_counter();reps=[c.reponses(w,forces) for w in omega];trep=time.perf_counter()-t
  champs=np.array([r['champ'] for r in reps]);j=juger(d,m,metrique,champs,ref)
  bounds={k:max((v if v is not None else float('inf')) for r in reps for v in r['bornes'][k]['relatives']) for k in ('masse','deformation')}
  resume=dict(n=n,blocs=blocs,dimension=r.taille_reduite,qr=tqr,selection=tmode,certificat=tc,krylov=tk,controle=tctrl,reponses=trep,majorants=bounds,juge=j,retours=[{k:v for k,v in r.items() if k not in ('champ','coordonnees')} for r in reps])
  def serialiser(v):return v.tolist() if hasattr(v,'tolist') else str(v)
  resultats.append(resume)
  (out/'resultats.json').write_text(json.dumps(resultats,default=serialiser,ensure_ascii=False,indent=2))
  print(n,blocs,'dim',r.taille_reduite,'cert',tc,'totalindicatif',tqr+tmode+tc+tk+tctrl+trep,'bornes',bounds,'erreurs',j['maxima_operateurs'],flush=True)
(out/'manifest.json').write_text(json.dumps({str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file()},indent=2))
