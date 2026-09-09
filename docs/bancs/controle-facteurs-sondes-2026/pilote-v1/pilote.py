import os,sys,json,time,hashlib,shutil
from pathlib import Path
os.sched_setaffinity(0,{8})
sys.dont_write_bytecode=True
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
out=Path('/tmp/vinkulum-pilote-controle-facteurs-2026-v1');out.mkdir(exist_ok=False)
ci=Path('/tmp/vinkulum-energie-native-travail/ci')
(out/'sources').mkdir()
for p in ci.glob('*.py'):
 if p.name in ('controle_facteurs.py','experience_inertie_contrainte.py','experience_krylov_contraint.py','inertie_complement_dirigee.py','krylov_contraint.py','controle_complement.py','inverse_contrainte_energie.py','racine_masse_blocs.py','trace_complement_dirigee.py','condensation_energie.py','confronte_ports_exudyn.py','reference_hcb_exudyn.py','oracle_champs_ports.py','reference_ports_precision.py','experience_reduction_native.py'):
  shutil.copy2(p,out/'sources'/p.name)
shutil.copy2(__file__,out/'pilote.py')
sys.path.insert(0,str(out/'sources'))
from confronte_ports_exudyn import lire,juger
from condensation_energie import CondensationEnergie
from inertie_complement_dirigee import certifier_inertie_complement
from krylov_contraint import KrylovContraint
from controle_complement import ControleComplement
from controle_facteurs import ControleFacteurs
from experience_inertie_contrainte import controle_oracle,jsonable,environnement
(out/'environnement.json').write_text(json.dumps(environnement('inertie_controle'),ensure_ascii=False,indent=2))
results=[]
for n in (32,128,512):
 d,m,metrique,forces,omega=lire(Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee')/f'n{n}-f40.npz')
 t=time.perf_counter();qr=CondensationEnergie(d,np.arange(6*n-6),np.arange(6*n-6,6*n),metrique);tqr=time.perf_counter()-t
 mi=m[qr.i][:,qr.i].tocsr();root=np.sqrt(mi.diagonal())
 def appliquer(v):return (root[:,None]*qr.solve(root[:,None]*np.asarray(v).reshape(-1,1))).ravel()
 t=time.perf_counter();_,vv=eigsh(LinearOperator(mi.shape,matvec=appliquer,dtype=float),k=1,which='LA',tol=1e-11,ncv=20,v0=np.random.default_rng(107+n).normal(size=len(qr.i)));phi=vv/root[:,None];b=mi@phi;ts=time.perf_counter()-t
 t=time.perf_counter();cert=certifier_inertie_complement(d[:,qr.i],mi,b,(2*np.pi*80)**2,precision=32);tc=time.perf_counter()-t
 t=time.perf_counter();r=KrylovContraint(qr,m,b,phi,cert['lambda_min'],omega[-1],blocs=4);tk=time.perf_counter()-t
 np.savez_compressed(out/f'n{n}-bases.npz',B=b,Phi=phi,base=r.base)
 (out/f'n{n}-certificat.json').write_text(json.dumps(cert,ensure_ascii=False,indent=2))
 ref=np.load(Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee')/f'n{n}-f40.ref.npy')
 ancien=None
 for variant in ('ancien','gram','qr'):
  case=out/f'n{n}-{variant}';case.mkdir()
  res=dict(n=n,variante=variant,partage_preparation=True,qr=tqr,selection=ts,certificat=tc,krylov=tk)
  try:
   t=time.perf_counter();c=ControleComplement(r) if variant=='ancien' else ControleFacteurs(r,extension=variant);prep=time.perf_counter()-t
   t=time.perf_counter();reps=[c.reponses(w,forces) for w in omega];trep=time.perf_counter()-t
   champs=np.array([v['champ'] for v in reps]);np.save(case/'champs.npy',champs)
   juge=juger(d,m,metrique,champs,ref);co=controle_oracle(d,m,champs,ref,reps)
   if variant=='ancien':ancien=champs
   bounds={k:max(v if v is not None else float('inf') for rep in reps for v in rep['bornes'][k]['relatives']) for k in ('masse','deformation')}
   res.update(statut='termine',controle_preparation=prep,reponses=trep,total_indicatif=tqr+ts+tc+tk+prep+trep,juge=juge,controle_oracle=co,majorants=bounds,champs_identiques_ancien=bool(np.array_equal(champs,ancien)),facteurs=getattr(c,'facteurs',None),retours=[{k:v for k,v in rep.items() if k not in ('champ','coordonnees')} for rep in reps])
   print(n,variant,'prep',prep,'reponses',trep,'totalindic',res['total_indicatif'],'controle',co['accepte_toutes_frequences'],'champsidentiques',res['champs_identiques_ancien'],'majorants',bounds,flush=True)
  except Exception as exc:
   res.update(statut='refus',erreur=type(exc).__name__,message=str(exc));print(n,variant,'REFUS',str(exc),flush=True)
  (case/'resultat.json').write_text(json.dumps(jsonable(res),ensure_ascii=False,indent=2,allow_nan=False))
  results.append({k:v for k,v in res.items() if k not in ('retours','controle_oracle','juge')})
  (out/'bilan.json').write_text(json.dumps(jsonable(results),ensure_ascii=False,indent=2,allow_nan=False))
(out/'manifest.json').write_text(json.dumps({str(p.relative_to(out)):dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),octets=p.stat().st_size) for p in out.rglob('*') if p.is_file()},indent=2))
