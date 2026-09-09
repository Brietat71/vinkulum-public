"""Adjoint discret : roue 0.8.2 contre produits implicites de la 0.9.0.

Processus frais, un échauffement conservé et trois répétitions, CPU fixé et
un fil demandé. Les familles mesurent un produit de pas puis des gradients
de trajectoire à un paramètre partagé ou une rigidité par élément.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import signal
import statistics
import subprocess
import sys
import time
from zipfile import ZipFile

ROOT=Path(__file__).resolve().parents[1]
THREADS={k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','RAYON_NUM_THREADS')}
CASES=[('pas',n) for n in (8,30,60,120,240)]+[(family,n) for family in ('partage','par_poutre') for n in (8,30,60,120)]
VARIANTS=('ancien','matrice','produits')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def wheel_identity(path):
    with ZipFile(path) as z:
        ext, = [s for s in z.namelist() if s.endswith('.so')]
        return dict(name=path.name,sha256=sha(path),
            extension_sha256=hashlib.sha256(z.read(ext)).hexdigest(),
            adjoint_python_sha256=hashlib.sha256(z.read('vinkulum/adjoint_temps.py')).hexdigest())


def worker(family, ne, variant, path):
    import numpy as np
    import vinkulum
    from vinkulum import Noyau, _vinkulum
    from vinkulum import adjoint_temps as at
    count=0
    def model(p):
        nonlocal count
        count+=1
        # Même modèle matériel que le pendule flexible de la vérification.
        L,b,h,rho,E=.6,.02,.004,2700.,70e9
        area=b*h; jt=b*h**3/12+h*b**3/12; dl=L/ne
        n=Noyau([0.,0.,-9.81])
        for i in range(ne+1):
            m=rho*area*dl*(.5 if i in (0,ne) else 1.)
            j=np.diag([rho*jt*dl,rho*b*h**3/12*dl,rho*h*b**3/12*dl])
            n.corps(str(i),m,j.ravel().tolist(),[i*dl,0.,0.])
        n.liaison('pivot',None,0,pa=[0.,0.,0.],bloque_r=[0,2])
        for i in range(ne):
            ei=float(p[i] if family=='par_poutre' else p[0])
            n.poutre(str(i),i,i+1,E*area,E/2.6*area,E/2.6*jt,ei)
        return n
    def rp(n,p):
        columns=[np.asarray(n.d_residu_poutre(i,'ei')) for i in range(ne)]
        return np.column_stack(columns) if family=='par_poutre' else np.sum(columns,axis=0)[:,None]
    def rtp(n,p,ud,psi):
        local=np.asarray(n.d_residu_poutres_transpose(psi[:len(ud)].tolist()))
        return local if family=='par_poutre' else np.array([sum(local)])
    p=np.full(ne if family=='par_poutre' else 1,.75)
    kwargs={} if variant=='ancien' else dict(creux=None)
    if variant=='produits':
        kwargs.update(d_residu_transpose=rtp,d_initial_transpose=lambda p,mu:np.zeros(len(p)))
    pont=at.PontNoyau(model,None if variant=='produits' else rp,**kwargs)
    initial=model(p)
    initial_rp=max(float(np.max(abs(np.asarray(initial.d_residu_poutre(i,'ei'))))) for i in range(ne))
    if initial_rp!=0.:
        raise ValueError('la preuve de dérivée initiale nulle ne tient pas')
    # p ne change ni M, ni G, ni q0/u0, et chaque R_p(q0)=0 : da0/dp=0.
    count=0; forward_seconds=0.; backward_seconds=0.; backward_calls=0
    original_forward=pont.aller
    def forward(*args):
        nonlocal forward_seconds
        start=time.perf_counter();out=original_forward(*args)
        forward_seconds+=time.perf_counter()-start
        return out
    pont.aller=forward
    method='_pas' if variant=='ancien' else '_recul'
    original_backward=getattr(pont,method)
    def backward(*args):
        nonlocal backward_seconds,backward_calls
        start=time.perf_counter();out=original_backward(*args)
        backward_seconds+=time.perf_counter()-start;backward_calls+=1
        return out
    setattr(pont,method,backward)
    if family=='pas':
        n,fr,sch=pont.aller(p,2e-4,1e-4)
        args=(n,fr[0],fr[1],np.asarray(sch[0][0]),np.asarray(sch[0][1]),np.asarray(sch[1][0]),np.asarray(sch[1][1]),p)
        mu=np.random.default_rng(918).normal(size=24*(ne+1))
        start=time.perf_counter()
        if variant=='ancien':
            a,b=pont._pas(*args);prev,grad=a.T@mu,b.T@mu
        else:
            prev,grad=pont._recul(*args,mu)
        seconds=time.perf_counter()-start
        values=np.concatenate((prev,grad)).tolist()
        final=fr[-1][1]
    else:
        def dj(n,fr):
            d=np.zeros(24*(ne+1));d[6*ne+2]=1.
            return d
        start=time.perf_counter()
        values=pont.gradient(p,.02,.001,dj).tolist()
        seconds=time.perf_counter()-start
        # Relecture indépendante hors chronomètre du gradient.
        n=model(p);n.simule(.02,.001);final=n.etat()[1]
        count-=1
    result=dict(version=vinkulum.__version__,family=family,elements=ne,variant=variant,
        core_seconds=seconds,forward_seconds=forward_seconds,backward_seconds=backward_seconds,
        backward_calls=backward_calls,build_count=count,values=values,final_positions=final,
        initial_parameter_residual_max=initial_rp,linearizations=pont.n_kcmz,
        adjoint_backward_error=getattr(pont,'residu_adjoint_max',None),
        linear_solves=getattr(pont,'n_resolutions_adjoint',None),
        extension_sha256=sha(_vinkulum.__file__),adjoint_python_sha256=sha(at.__file__),
        numpy=np.__version__,python=sys.version,affinity=sorted(os.sched_getaffinity(0)),
        threads={k:os.environ.get(k) for k in THREADS})
    Path(path).write_text(json.dumps(result,allow_nan=False))


def analyse(report):
    import numpy as np
    out=dict(definition=__doc__,metadata=report['metadata'],configurations=[],failures=[])
    runs=report['runs'];identities=[(r['family'],r['elements'],r['variant'],r['rep']) for r in runs]
    if len(set(identities))!=len(identities):raise ValueError('essai dupliqué')
    expected={(f,n,v,r) for f,n in CASES for v in VARIANTS for r in (-1,0,1,2)}
    if set(identities)!=expected or list(map(tuple,report['metadata']['cases']))!=CASES:
        raise ValueError('campagne incomplète ou configuration inconnue')
    if len(report['metadata']['affinity'])!=1:raise ValueError('CPU non fixé')
    environments=set()
    for family,n in report['metadata']['cases']:
        old=[r for r in runs if (r['family'],r['elements'],r['variant'])==(family,n,'ancien') and r['rep']>=0 and r['ok']]
        for variant in VARIANTS:
            samples=[r for r in runs if (r['family'],r['elements'],r['variant'])==(family,n,variant)]
            if sorted(r['rep'] for r in samples)!=[-1,0,1,2]:raise ValueError('répétitions ou échauffement manquants')
            for r in samples:
                if not r['ok']:
                    out['failures'].append({k:r[k] for k in ('family','elements','variant','rep','error')})
                    continue
                d=r['result']
                if (d['family'],d['elements'],d['variant'])!=(family,n,variant):raise ValueError('identité incohérente')
                version='0.8.2' if variant=='ancien' else '0.9.0'
                wheel=report['metadata']['wheels']['ancien' if variant=='ancien' else 'nouveau']
                if d['version']!=version or any(d[key]!=wheel[key] for key in ('extension_sha256','adjoint_python_sha256')):
                    raise ValueError('version ou sources chargées incorrectes')
                if d['affinity']!=report['metadata']['affinity'] or d['threads']!=THREADS:
                    raise ValueError('affinité ou parallélisme différent')
                environments.add((d['python'],d['numpy']))
                if not np.isfinite(d['values']).all() or not np.isfinite(d['final_positions']).all():raise ValueError('observable non fini')
                expected_parameters=n if family=='par_poutre' else 1
                expected_values=24*(n+1)+1 if family=='pas' else expected_parameters
                if np.shape(d['values'])!=(expected_values,) or np.shape(d['final_positions'])!=(n+1,3):raise ValueError('forme des observables incorrecte')
                if any(not np.isfinite(v) or v<=0 for v in (d['core_seconds'],d['forward_seconds'],d['backward_seconds'],r['wall_seconds'],r['rss_kib'])):
                    raise ValueError('mesure de coût invalide')
                if d['initial_parameter_residual_max']!=0.:raise ValueError('dérivée initiale non nulle')
                expected_calls=1 if family=='pas' else 20
                if d['backward_calls']!=expected_calls or d['linearizations']!=expected_calls:raise ValueError('pas adjoints manquants')
                if variant!='ancien' and not 0<=d['adjoint_backward_error']<=1e-11:raise ValueError('résidu adjoint invalide')
                if variant!='ancien' and not expected_calls<=d['linear_solves']<=4*expected_calls:raise ValueError('nombre de résolutions incohérent')
                expected_builds=1 if family=='pas' or variant=='produits' else 1+2*expected_parameters
                if d['build_count']!=expected_builds:raise ValueError('construction initiale manquante')
            good=[r for r in samples if r['ok'] and r['rep']>=0]
            c=dict(family=family,elements=n,variant=variant,successes=len(good),verified=False)
            if good:
                ds=[r['result'] for r in good]
                c.update(seconds=[d['core_seconds'] for d in ds],median_seconds=statistics.median(d['core_seconds'] for d in ds),
                    median_process_seconds=statistics.median(r['wall_seconds'] for r in good),rss_kib=[r['rss_kib'] for r in good],
                    median_backward_seconds=statistics.median(d['backward_seconds'] for d in ds),build_counts=[d['build_count'] for d in ds],
                    values=ds[0]['values'],extension_sha256=ds[0]['extension_sha256'],adjoint_python_sha256=ds[0]['adjoint_python_sha256'])
                if old:
                    ref=np.asarray(old[0]['result']['values']); scale=float(np.max(abs(ref),initial=0.))
                    error=max(float(np.max(abs(np.asarray(d['values'])-ref),initial=0.)) for d in ds)
                    pose=max(float(np.max(abs(np.asarray(d['final_positions'])-old[0]['result']['final_positions']),initial=0.)) for d in ds)
                    # Le grand covecteur de pose ne doit pas masquer une erreur
                    # sur les accélérations ou le gradient de paramètre.
                    splits=[0,6*(n+1),12*(n+1),18*(n+1),24*(n+1),len(ref)] if family=='pas' else [0,len(ref)]
                    groups=[]
                    for a,b in zip(splits[:-1],splits[1:]):
                        size=float(np.max(abs(ref[a:b]),initial=0.))
                        err=max(float(np.max(abs(np.asarray(r['result']['values'])[a:b]-ref[a:b]),initial=0.)) for r in samples if r['ok'])
                        groups.append(dict(scale=size,max_difference=err,verified=err<=2e-11+3e-9*size))
                    c.update(max_difference=error,relative_difference=error/scale if scale else None,max_position_difference=pose,
                        groups=groups,verified=len(old)==3 and len(good)==3 and all(r['ok'] for r in samples) and all(g['verified'] for g in groups) and pose<=1e-12)
                    c['speedup']=statistics.median(r['result']['core_seconds'] for r in old)/c['median_seconds']
            out['configurations'].append(c)
    if len(environments)>1:raise ValueError('environnements Python/NumPy différents')
    for n in (8,30,60,120):
        for variant in VARIANTS:
            for rep in (-1,0,1,2):
                pair=[next(r for r in runs if (r['family'],r['elements'],r['variant'],r['rep'])==(f,n,variant,rep))
                      for f in ('partage','par_poutre')]
                if all(r['ok'] for r in pair):
                    shared=pair[0]['result']['values'][0]
                    independent=sum(pair[1]['result']['values'])
                    if abs(shared-independent)>2e-11+3e-9*abs(shared):
                        raise ValueError('la somme des gradients locaux ne retrouve pas le paramètre partagé')
    out['all_verified']=not out['failures'] and all(c['verified'] for c in out['configurations'])
    return out


def main():
    if len(sys.argv)>1 and sys.argv[1]=='--worker':
        worker(sys.argv[2],int(sys.argv[3]),sys.argv[4],sys.argv[5]);return
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ancien',type=Path);p.add_argument('--nouveau',type=Path)
    p.add_argument('--ancien-roue',type=Path);p.add_argument('--nouveau-roue',type=Path)
    p.add_argument('--sortie',type=Path);p.add_argument('--archive',type=Path)
    p.add_argument('--verifier',type=Path);p.add_argument('--timeout',type=float,default=120.)
    args=p.parse_args()
    if args.verifier:
        report=json.loads(gzip.decompress(args.verifier.read_bytes()))
        if hashlib.sha256(report['source'].encode()).hexdigest()!=report['metadata']['source_sha256']:raise ValueError('source altérée')
        result=analyse(report)
        if result!=json.loads(args.verifier.with_name(args.verifier.name.replace('-essais.json.gz','.json')).read_bytes()):raise ValueError('bilan altéré')
        if not result['all_verified']:raise ValueError('comparaison non vérifiée')
        print(f"{len(report['runs'])} essais : gradients, produits, paliers et résidus vérifiés.");return
    if not all((args.ancien,args.nouveau,args.ancien_roue,args.nouveau_roue,args.sortie,args.archive)):p.error('exécutables, roues et sorties requis')
    if any(Path(str(args.archive)+suffix).exists() for suffix in ('.json','-essais.json.gz')):
        p.error('archive existante : choisir un nouveau préfixe')
    args.sortie=args.sortie.resolve();args.sortie.mkdir(exist_ok=False)
    source=Path(__file__).read_text()
    report=dict(metadata=dict(created_utc=datetime.now(timezone.utc).isoformat(),command=sys.argv,cases=CASES,
        affinity=sorted(os.sched_getaffinity(0)),threads=THREADS,
        wheels=dict(ancien=wheel_identity(args.ancien_roue),nouveau=wheel_identity(args.nouveau_roue)),
        cpu=next(s for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')),
        source_sha256=hashlib.sha256(source.encode()).hexdigest()),source=source,runs=[])
    for family,n in CASES:
        for rep in (-1,0,1,2):
            for variant in (VARIANTS if rep%2==0 else VARIANTS[::-1]):
                dest=args.sortie/f'{family}-{n}-{variant}-{rep}';dest.mkdir()
                py=args.ancien if variant=='ancien' else args.nouveau
                command=['/usr/bin/time','-f','%M','-o',str(dest/'rss.txt'),str(py),str(Path(__file__).resolve()),'--worker',family,str(n),variant,str(dest/'result.json')]
                start=time.perf_counter();timed_out=False
                proc=subprocess.Popen(command,env=dict(os.environ,**THREADS),cwd=dest,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,start_new_session=True)
                try:stdout,_=proc.communicate(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    timed_out=True;os.killpg(proc.pid,signal.SIGKILL);stdout,_=proc.communicate()
                elapsed=time.perf_counter()-start
                r=dict(family=family,elements=n,variant=variant,rep=rep,wall_seconds=elapsed,ok=proc.returncode==0 and not timed_out,
                    command=command,stdout=stdout.decode(),timed_out=timed_out)
                (dest/'stdout.txt').write_bytes(stdout)
                if r['ok']:
                    r.update(result=json.loads((dest/'result.json').read_bytes()),rss_kib=int((dest/'rss.txt').read_text()))
                else:r['error']=f'exit {proc.returncode}; timeout={timed_out}'
                report['runs'].append(r)
                print(family,n,variant,rep,r['ok'],round(elapsed,6),flush=True)
                (args.sortie/'resultats.json').write_text(json.dumps(report,allow_nan=False))
    if sha(__file__)!=report['metadata']['source_sha256']:raise ValueError('source changée pendant les mesures')
    bilan=analyse(report)
    Path(str(args.archive)+'-essais.json.gz').write_bytes(gzip.compress(json.dumps(report,ensure_ascii=False,allow_nan=False,separators=(',',':')).encode(),mtime=0))
    Path(str(args.archive)+'.json').write_text(json.dumps(bilan,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    if not bilan['all_verified']:raise ValueError('comparaison non vérifiée ; échecs conservés')


if __name__=='__main__':main()
