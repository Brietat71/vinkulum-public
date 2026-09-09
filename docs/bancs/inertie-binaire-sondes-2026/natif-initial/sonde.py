import os,sys,time,json,hashlib,shutil
from pathlib import Path
os.sched_setaffinity(0,{8});sys.dont_write_bytecode=True
import numpy as np
out=Path('/tmp/vinkulum-sonde-inertie-compilee-2026-v1');out.mkdir(exist_ok=False)
ci=Path('/tmp/vinkulum-energie-native-travail/ci');sys.path.insert(0,str(ci))
from inertie_binaire_compile import compiler,BibliothequeInertie,certifier_inertie_compilee
bib=BibliothequeInertie(compiler(out/"build"))
from confronte_ports_exudyn import lire
for name in ('inertie_binaire_compile.py','inertie_binaire_native.cpp','inertie_complement_dirigee.py'):shutil.copy2(ci/name,out/name)
shutil.copy2(__file__,out/'sonde.py')
results=[]
for n in (32,128,512):
    d,m,*_=lire(Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee')/f'n{n}-f40.npz');di=d[:,:6*n-6];mi=m[:6*n-6,:6*n-6].tocsr()
    base=np.load(Path('/tmp/vinkulum-sonde-separateurs-2026')/f'n{n}-entrees.npz');b,perm=base['B'],base['permutation']
    shutil.copy2(Path('/tmp/vinkulum-sonde-separateurs-2026')/f'n{n}-entrees.npz',out/f'n{n}-entrees.npz')
    for kind,p in [('naturel',None),('separateurs',perm)]:
        t=time.perf_counter();r=dict(n=n,ordre=kind)
        try:r.update(statut='certifie',certificat=certifier_inertie_compilee(di,mi,b,(2*np.pi*80)**2,permutation=p,bibliotheque=bib))
        except Exception as e:r.update(statut='refus',message=str(e),diagnostic=getattr(e,'bilan',{}))
        r['temps_s']=time.perf_counter()-t;results.append(r)
        print(n,kind,r['statut'],r['temps_s'],r.get('message',''),flush=True)
        (out/'resultats.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
(out/'manifest.json').write_text(json.dumps({str(p.relative_to(out)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'octets':p.stat().st_size} for p in out.rglob('*') if p.is_file()},indent=2))
