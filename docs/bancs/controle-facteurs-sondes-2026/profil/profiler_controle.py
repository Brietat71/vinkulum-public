"""Profil exploratoire du contrôleur original ; aucun comparatif de solveurs."""
import os
THREAD_KEYS=('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','RAYON_NUM_THREADS',
             'BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')
for key in THREAD_KEYS:
    os.environ[key]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
os.sched_setaffinity(0,{8})

import cProfile
from collections import defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import inspect
import io
import json
from pathlib import Path
import platform
import pstats
import shutil
import subprocess
import sys
import time

import numpy as np
import scipy

ROOT=Path('/tmp/vinkulum-energie-native-travail')
OUT=Path(__file__).resolve().parent
SOURCE=OUT/'sources'
SOURCE.mkdir(exist_ok=False)
NAMES=('controle_complement.py','krylov_contraint.py','inverse_contrainte_energie.py',
       'racine_masse_blocs.py','condensation_energie.py','confronte_ports_exudyn.py')
for name in NAMES:
    shutil.copy2(ROOT/'ci'/name,SOURCE/name)
sys.path.insert(0,str(SOURCE))
import controle_complement as controle_module
from controle_complement import ControleComplement
from krylov_contraint import KrylovContraint
from condensation_energie import CondensationEnergie
from confronte_ports_exudyn import lire
import vinkulum
from importlib.metadata import distribution

ENTRY=Path('/tmp/vinkulum-confrontation-ports-0.10.0-corrigee/n512-f40.npz')
BASIS=ROOT/'docs/bancs/inertie-contraint-2026/essais/n512-f40-inertie_controle-passage1/bases.npz'
META=Path('/tmp/vinkulum-campagne-inertie-contrainte-2026/essais/n512-f40-inertie_controle-passage1/resultat.json')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, data):
    Path(path).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def fingerprint(responses):
    h=hashlib.sha256()
    maximum={nom:0. for nom in ('masse','deformation')}
    for response in responses:
        assert response['certification_machine'] is False
        for value in (response['champ'],response['coordonnees']):
            h.update(np.asarray(value,dtype='<f8').tobytes())
        for nom in maximum:
            b=response['bornes'][nom]
            for key in ('normes','absolues','relatives'):
                a=np.asarray(b[key],dtype='<f8')
                assert np.all(np.isfinite(a))
                h.update(a.tobytes())
            maximum[nom]=max(maximum[nom],max(b['relatives']))
    return dict(sha256_champs_coordonnees_normes_bornes=h.hexdigest(),majorants_maximaux=maximum)

def sweep(controller,omega,forces):
    return [controller.reponses(float(w),forces) for w in omega]

def measure_sweep(controller,omega,forces):
    start=time.perf_counter()
    responses=sweep(controller,omega,forces)
    duration=time.perf_counter()-start
    identity=fingerprint(responses)
    del responses
    return dict(wall_s=duration,**identity)

def main():
    d,m,metrique,forces,omega=lire(ENTRY)
    archived=json.loads(META.read_text())
    with np.load(BASIS,allow_pickle=False) as z:
        b,phi=z['B'].copy(),z['Phi'].copy()
    gamma=archived['certificat_complement']['lambda_min']
    assert hashlib.sha256(np.asarray(b,dtype='<f8').tobytes()).hexdigest()==archived['certificat_complement']['contraintes_sha256']
    assert sha(ENTRY)==archived['entree_sha256']
    shutil.copy2(BASIS,OUT/'bases_originales.npz')
    write(OUT/'certificat_offert.json',archived['certificat_complement'])
    package=Path(vinkulum.__file__).parent
    installed=package/'_ports/condensation_energie.py'
    shutil.copy2(installed,OUT/'condensation_energie_distribuee.py')
    record=distribution('vinkulum').read_text('RECORD')
    environment=dict(date_utc=datetime.now(timezone.utc).isoformat(),
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        python=sys.version,numpy=np.__version__,scipy=scipy.__version__,vinkulum=vinkulum.__version__,
        platform=platform.platform(),cpu=sorted(os.sched_getaffinity(0)),
        threads={k:os.environ.get(k) for k in THREAD_KEYS},
        sources_sha256={name:sha(SOURCE/name) for name in NAMES},
        installed_python_sha256={str(p.relative_to(package)):sha(p) for p in package.rglob('*.py')},
        extension_sha256={p.name:sha(p) for p in package.glob('*.so')},
        record_sha256=hashlib.sha256(record.encode()).hexdigest(),
        entry=dict(path=str(ENTRY),sha256=sha(ENTRY)),
        basis=dict(path=str(BASIS),sha256=sha(BASIS)),
        meta=dict(path=str(META),sha256=sha(META)))
    write(OUT/'environnement.json',environment)
    print('CALCUL CPU8 : préparation QR/Krylov/contrôle ; B/Phi et certificat offerts.',flush=True)
    times={};start=time.perf_counter()
    qr=CondensationEnergie(d,np.arange(d.shape[1]-6),np.arange(d.shape[1]-6,d.shape[1]),metrique)
    times['qr']=time.perf_counter()-start;start=time.perf_counter()
    reduction=KrylovContraint(qr,m,b,phi,gamma,omega[-1],blocs=4,max_directions=128)
    times['krylov']=time.perf_counter()-start;start=time.perf_counter()
    controller=ControleComplement(reduction,profondeur=8)
    times['controle']=time.perf_counter()-start
    print('Préparé, dimension',reduction.taille_reduite,'; quatre balayages257fréquences, six charges.',flush=True)
    warm=measure_sweep(controller,omega,forces)
    baseline=measure_sweep(controller,omega,forces)
    profile=cProfile.Profile()
    profile.enable();profiled=measure_sweep(controller,omega,forces);profile.disable()
    profile.dump_stats(str(OUT/'reponses_originales.prof'))
    stream=io.StringIO()
    statistics=pstats.Stats(profile,stream=stream).sort_stats('cumulative')
    statistics.print_stats(45)
    (OUT/'cprofile-cumulatif.txt').write_text(stream.getvalue())
    rows=[]
    for (path,line,name),(primitive,calls,self_time,cumulative,callers) in statistics.stats.items():
        rows.append(dict(path=path,line=line,function=name,calls=calls,primitive_calls=primitive,
                         self_s=self_time,cumulative_s=cumulative,
                         callers=[dict(path=k[0],line=k[1],function=k[2],stats=list(v) if isinstance(v,tuple) else v)
                                  for k,v in callers.items()]))
    rows.sort(key=lambda r:r['cumulative_s'],reverse=True)

    # Attribution exploratoire par lignes du contrôleur uniquement. Le temps
    # des appels enfants est inclus ; les lignes ne sont pas des microbenchmarks.
    code=ControleComplement.reponses.__code__
    line_times=defaultdict(int);line_counts=defaultdict(int);state={}
    def tracing(frame,event,arg):
        if frame.f_code is not code:
            return None
        now=time.perf_counter_ns()
        if event=='call':
            state[id(frame)]=(None,now)
        elif event in ('line','return'):
            previous,last=state[id(frame)]
            if previous is not None:
                line_times[previous]+=now-last
                line_counts[previous]+=1
            if event=='return':state.pop(id(frame))
            else:state[id(frame)]=(frame.f_lineno,now)
        return tracing

    norm_original=controle_module.norme
    norm_details=defaultdict(lambda:dict(calls=0,wall_ns=0))
    def norm_profiled(a):
        caller=inspect.currentframe().f_back.f_lineno
        key=(caller,tuple(a.shape))
        start=time.perf_counter_ns();result=norm_original(a);elapsed=time.perf_counter_ns()-start
        entry=norm_details[key];entry['calls']+=1;entry['wall_ns']+=elapsed
        return result
    controle_module.norme=norm_profiled
    try:
        sys.settrace(tracing)
        line_profiled=measure_sweep(controller,omega,forces)
    finally:
        sys.settrace(None)
        controle_module.norme=norm_original
    source=(SOURCE/'controle_complement.py').read_text().splitlines()
    lines=[dict(line=line,calls=line_counts[line],wall_s=ns/1e9,text=source[line-1])
           for line,ns in sorted(line_times.items(),key=lambda item:item[1],reverse=True)]
    norms=[dict(caller_line=line,shape=list(shape),calls=r['calls'],wall_s=r['wall_ns']/1e9)
           for (line,shape),r in norm_details.items()]
    for run in (baseline,profiled,line_profiled):
        assert run['sha256_champs_coordonnees_normes_bornes']==warm['sha256_champs_coordonnees_normes_bornes']
    result=dict(status='PROFIL EXPLORATOIRE NON COMPARATIF',
        notes=['B/Phi et certificat80Hz réutilisés depuis un essai antérieur.',
               'Même candidat préparé pour quatre balayages ; aucun autre solveur exécuté.',
               'Instrumentation cProfile/trace/wrapper modifie le temps ; attribution indicative.',
               'Les durées du balayage incluent uniquement les réponses ; les empreintes sont hors de wall_s.',
               'cProfile inclut toutefois le travail de fingerprint dans son périmètre global.',
               'Le profil par lignes inclut ses surcoûts ; les normes détaillées recouvrent ces mêmes lignes.'],
        dimensions=dict(physique=d.shape[1],conservee=reduction.p+reduction.r,krylov=reduction.v.shape[1],
                        frequences=len(omega),charges=forces.shape[1]),
        gamma=gamma,preparations_indicatives_s=times,
        warm=warm,baseline=baseline,cprofile_sweep=profiled,line_sweep=line_profiled,
        cprofile=rows,lines=lines,norms=norms)
    write(OUT/'profil.json',result)
    assert all(sha(ROOT/'ci'/name)==environment['sources_sha256'][name] for name in NAMES)
    print('FIN CALCUL CPU8. Identité bit à bit des champs/coordonnées/normes/majorants sur quatre balayages.',flush=True)
    print(json.dumps(dict(baseline=baseline['wall_s'],cprofile=profiled['wall_s'],trace=line_profiled['wall_s'],
                          top_lines=lines[:12],norms=norms),ensure_ascii=False),flush=True)

if __name__=='__main__':main()
