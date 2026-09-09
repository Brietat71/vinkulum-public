import sys,json,time,gzip
from pathlib import Path
sys.path.insert(0,'/tmp/vinkulum-energie-native-travail/ci')
import numpy as np
from vinkulum._ports.condensation_energie import CondensationEnergie
from confronte_ports_exudyn import lire,juger
from krylov_contraint import KrylovContraint
from controle_complement import ControleComplement
src=Path('/tmp/vinkulum-energie-native-travail/docs/bancs/complement-spectral-2026')
ancien=json.loads(gzip.decompress((src/'rapport.json.gz').read_bytes()))
out=[]
for n in (32,128,512):
 d,m,metric,f,omega=lire(src/f'entrees/n{n}-f40.npz'); ni=d.shape[1]-6
 qr=CondensationEnergie(d,np.arange(ni),np.arange(ni,d.shape[1]),metric)
 with np.load(src/f'directions/n{n}-directions.npz') as a: b,phi=a['b'][:,:1],a['phi'][:,:1]
 lam=next(c for c in ancien['cas'] if c['n']==n)['certificats'][0]['lambda_min']
 ref=np.load(f'/tmp/vinkulum-confrontation-ports-0.10.0-corrigee/n{n}-f40.ref.npy',mmap_mode='r')
 for blocs in (2,4,6):
  r=KrylovContraint(qr,m,b,phi,lam,float(omega[-1]),blocs=blocs)
  ctrl=ControleComplement(r,profondeur=8)
  t=time.perf_counter(); reps=[ctrl.reponses(w,f) for w in omega]; temps=time.perf_counter()-t
  champs=np.array([a['champ'] for a in reps]); juge=juger(d,m,metric,champs,ref)
  bornes={nom:[a['bornes'][nom]['relatives'] for a in reps] for nom in ('masse','deformation')}
  maxbornes={nom:None if any(v is None for a in b for v in a) else max(v for a in b for v in a) for nom,b in bornes.items()}
  out.append(dict(n=n,blocs=blocs,taille=r.taille_reduite,delta=ctrl.delta,bS=ctrl.borne_schur,prep_controle=ctrl.preparation_s,reponses=temps,cm=ctrl.cm,cd=ctrl.cd,borne_max=maxbornes,bornes=bornes,juge=juge,termes=ctrl.termes,marge_min=min(a['marge'] for a in reps)))
  Path('/tmp/vinkulum-pilot-controle-complement-repare.json').write_text(json.dumps(out,allow_nan=False))
  print(n,blocs,r.taille_reduite,ctrl.delta,ctrl.borne_schur,maxbornes,juge['accepte'],juge['maxima_operateurs'],ctrl.preparation_s,temps,flush=True)
