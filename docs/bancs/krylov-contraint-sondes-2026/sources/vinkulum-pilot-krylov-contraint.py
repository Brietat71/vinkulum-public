import sys, json, time
from pathlib import Path
sys.path.insert(0,'/tmp/vinkulum-energie-native-travail/ci')
import numpy as np
from vinkulum._ports.condensation_energie import CondensationEnergie
from confronte_ports_exudyn import lire, juger
from krylov_contraint import KrylovContraint
src=Path('/tmp/vinkulum-energie-native-travail/docs/bancs/complement-spectral-2026')
import gzip
ancien=json.loads(gzip.decompress((src/'rapport.json.gz').read_bytes()))
out=[]
for n in (32,128,512):
 d,m,metric,f,omega=lire(src/f'entrees/n{n}-f40.npz')
 ni=d.shape[1]-6
 qr=CondensationEnergie(d,np.arange(ni),np.arange(ni,d.shape[1]),metric)
 with np.load(src/f'directions/n{n}-directions.npz') as a: b,phi=a['b'][:,:1],a['phi'][:,:1]
 lam=next(c for c in ancien['cas'] if c['n']==n)['certificats'][0]['lambda_min']
 ref=np.load(f'/tmp/vinkulum-confrontation-ports-0.10.0-corrigee/n{n}-f40.ref.npy',mmap_mode='r')
 for blocs in (2,4,6):
  start=time.perf_counter()
  r=KrylovContraint(qr,m,b,phi,lam,float(omega[-1]),blocs=blocs)
  prep=time.perf_counter()-start
  start=time.perf_counter()
  x=np.array([r.reponses(w,f)['champ'] for w in omega])
  elapsed=time.perf_counter()-start
  j=juger(d,m,metric,x,ref)
  row=dict(n=n,blocs=blocs,taille=r.taille_reduite,prep=prep,reponses=elapsed,defaut_contrainte=r.defaut_contrainte,juge=j)
  out.append(row)
  Path('/tmp/vinkulum-pilot-krylov-contraint.json').write_text(json.dumps(out,allow_nan=False))
  print(n,blocs,r.taille_reduite,prep,elapsed,j['accepte'],j['maxima_operateurs'],flush=True)
