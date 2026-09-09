import os,sys,time,json,shutil,hashlib
from pathlib import Path
os.sched_setaffinity(0,{8});sys.dont_write_bytecode=True
import numpy as np
from scipy.sparse.csgraph import connected_components
out=Path('/tmp/vinkulum-sonde-separateurs-2026');out.mkdir(exist_ok=False)
ci=Path('/tmp/vinkulum-energie-native-travail/ci');sys.path.insert(0,str(ci))
from inertie_complement_dirigee import certifier_inertie_complement
from confronte_ports_exudyn import lire
shutil.copy2(__file__,out/'sonde.py');shutil.copy2(ci/'inertie_complement_dirigee.py',out/'inertie_complement_dirigee.py')
def ordre_separateurs(d,m):
    motif=d.copy().astype(bool);g=(motif.T@motif)+m.astype(bool);g=g.tocsr();g.setdiag(False);g.eliminate_zeros()
    adj=[set(map(int,g.indices[g.indptr[i]:g.indptr[i+1]])) for i in range(g.shape[0])]
    def bfs(nodes,start):
        visited={start};layers=[[start]]
        while True:
            frontier=set()
            for v in layers[-1]:frontier.update(adj[v]&nodes)
            frontier-=visited
            if not frontier:return layers,visited
            layers.append(sorted(frontier));visited|=frontier
    def rec(nodes):
        if len(nodes)<=8:return sorted(nodes)
        layers,visited=bfs(nodes,min(nodes))
        if visited!=nodes:return rec(visited)+rec(nodes-visited)
        layers,_=bfs(nodes,min(layers[-1]))
        counts=np.cumsum([len(l) for l in layers]);mid=int(np.searchsorted(counts,len(nodes)/2))
        sep=set(layers[mid]);left=set(v for layer in layers[:mid] for v in layer);right=nodes-sep-left
        return rec(left)+rec(right)+sorted(sep)
    return np.array(rec(set(range(g.shape[0]))),dtype=np.int64)
results=[]
for n in (32,128,512):
    d,m,*_=lire(Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee')/f'n{n}-f40.npz');di=d[:,:6*n-6];mi=m[:6*n-6,:6*n-6].tocsr()
    b=np.load(Path('/tmp/vinkulum-pilote-controle-facteurs-2026-v2')/f'n{n}-bases.npz')['B']
    t=time.perf_counter();perm=ordre_separateurs(di,mi);tp=time.perf_counter()-t
    np.savez_compressed(out/f'n{n}-entrees.npz',B=b,permutation=perm)
    for kind,p in [('naturel',None),('separateurs',perm)]:
      for precision in (16,20,32):
        t=time.perf_counter();r=dict(n=n,ordre=kind,precision=precision,ordre_s=tp if p is not None else 0)
        try:
            cert=certifier_inertie_complement(di,mi,b,(2*np.pi*80)**2,precision=precision,permutation=p)
            r.update(statut='certifie',certificat=cert)
        except Exception as e:r.update(statut='refus',message=str(e),diagnostic=getattr(e,'bilan',{}))
        r['temps_s']=time.perf_counter()-t;results.append(r)
        print(n,kind,precision,r['statut'],r['temps_s'],r.get('certificat',{}).get('preuve_kkt',{}).get('largeur_max'),r.get('message',''),flush=True)
        (out/'resultats.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
(out/'manifest.json').write_text(json.dumps({str(p.relative_to(out)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'octets':p.stat().st_size} for p in out.rglob('*') if p.is_file()},indent=2))
