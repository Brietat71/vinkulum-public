"""Confrontation mesurée, sans lire les sources du solveur MBDyn.

.venv314/bin/python ci/confronte_mbdyn.py --mbdyn /chemin/mbdyn \
    --benchmarks /chemin/tests/benchmarks --sortie /tmp/confrontation

Temps de processus complet (imports, construction, calcul, sorties), trois
répétitions après échauffement ; un seul calcul à la fois, un fil demandé.
Les références sont recalculées par CHAQUE moteur à deux raffinements.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CASES = {
    'six_barres': dict(deck='6barmech/6barmech', T=3., steps=[.004,.002,.001,.0005], ref=[.00025,.000125], labels=[0,1,2,101,102], sample=.004, unit='m', target=1e-4),
    'spatial': dict(deck='srskm/srscm', T=5., steps=[.002,.001,.0005,.00025], ref=[.0000625,.00003125], labels=[100,200,300], sample=.002, unit='m', target=1e-4),
    'andrews': dict(deck='andrewssqueezer/andrewssqueezer', T=.02, steps=[8e-6,4e-6,2e-6,1e-6], ref=[2.5e-7,1.25e-7], labels=[1600], sample=.0002, unit='rad', target=1e-3),
    'princeton': dict(deck='princeton/princeton', T=.5, steps=[10,20,40], ref=[80,160], labels=[], unit='m', target=1e-5, formulation='milieu'),
    'princeton_integree': dict(deck='princeton/princeton', T=.5, steps=[4,6,8,10,12,16,20,40], ref=[80,160], labels=[], unit='m', target=1e-5, formulation='integree'),
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def princeton(nel, formulation="milieu"):
    from vinkulum import Noyau
    n=Noyau([0.,0.,0.]); length=.508; dl=length/nel
    ids=[n.corps(f'n{i}', 1., [1.,0.,0.,0.,1.,0.,0.,0.,1.], [i*dl,0.,0.]) for i in range(nel+1)]
    n.liaison('enc',None,ids[0],bloque_t=[0,1,2],bloque_r=[0,1,2])
    for i in range(nel):
        options = dict(formulation=formulation) if formulation != "milieu" else {}
        n.poutre(f'e{i}',ids[i],ids[i+1],2.84191e6,6.40131e5,3.10338,36.2794,ei3=2.42873,ga3=9.03881e5,**options)
    return n,ids


def worker(case, step, output):
    import numpy as np
    from modeles_confrontation import six_barres, spatial
    spec=CASES[case]
    static='formulation' in spec
    if static:
        n,ids=princeton(int(step),spec['formulation'])
    elif case=='andrews':
        from vinkulum.andrews import monte
        n,ids=monte()
    else:
        n,ids={'six_barres':six_barres,'spatial':spatial}[case]()
    diagnostics=[]
    start=time.perf_counter()
    if static:
        previous=0.
        for k in range(1,51):
            load=.5*(1-np.cos(np.pi*k/50))
            force=8.896*(load-previous)/np.sqrt(2)
            n.effort(ids[-1],[0.,force,force],[0.,0.,0.])
            try:
                n.statique(tol=1e-8,iters=100,strict=True)
            except Exception:
                diagnostics.append(dict(palier_charge=k,load=float(load),**n.statique_info()))
                Path(output).write_text(json.dumps(dict(static_diagnostics=diagnostics)))
                raise
            diagnostics.append(dict(palier_charge=k,load=float(load),**n.statique_info()))
            if diagnostics[-1]['statut']!='tolerance':
                raise RuntimeError('palier statique sans convergence à la tolérance')
            previous=load
        values=[n.etat()[1][ids[-1]]]; times=[.5]
    else:
        initial=n.etat()
        tr=[initial, *n.simule(spec['T'],step,tous=1)]
        times=[x[0] for x in tr]
        if case=='andrews':
            i=ids['O_F']; angles=np.unwrap([np.arctan2(x[2][i][3],x[2][i][0]) for x in tr])
            values=(angles-angles[0])[:,None]
        else:
            values=np.array([[x[1][i] for i in ids] for x in tr]).reshape(len(tr),-1)
    solve=time.perf_counter()-start
    values=np.asarray(values)
    if not static:
        if times[-1]<spec['T']-1e-10:
            raise ValueError('Vinkulum incomplete trajectory')
        grid=np.linspace(0,spec['T'],round(spec['T']/spec['sample'])+1)
        values=np.array([np.interp(grid,times,values[:,i]) for i in range(values.shape[1])]).T
        times=grid
    Path(output).write_text(json.dumps(dict(times=list(times),values=values.tolist(), solve_seconds=solve,
        phi_max=max(map(abs,n.phi()),default=0.), stats=n.stats(),static_diagnostics=diagnostics)))


def replace(s,key,value):
    pattern=r'(?m)^(\s*)'+re.escape(key)+r'\s*:[^;]*;'
    s,count=re.subn(pattern,lambda m:m[1]+key+': '+str(value)+';',s)
    if count!=1:
        raise ValueError(f'{key}: {count} declarations')
    return s


def deck(benchmarks,case,step,dest):
    spec=CASES[case]; original=benchmarks/spec['deck']; s=original.read_text()
    # Includes are mathematical model data; no solver implementation is inspected.
    for include in original.parent.iterdir():
        if include.suffix in ('.ref','.set'):
            (dest/include.name).write_bytes(include.read_bytes())
    s=replace(s,'final time',spec['T'])
    s=replace(s,'time step',.01 if 'formulation' in spec else step)
    s=re.sub(r'(?m)^\s*default output\s*:[^;]*;', '',s)
    s=re.sub(r'(?m)^\s*output precision\s*:[^;]*;', '',s)
    s=s.replace('begin: control data;','begin: control data;\n default output: none, structural nodes;\n output precision: 16;')
    if 'formulation' in spec:
        count=int(step); beams=count//2
        (dest/'princeton.set').write_text(f'''set: const integer N_BEAM_NODES = {count};
set: const integer N_BEAMS = {beams};
set: const real THETA = 45.;
set: const real P = 8.896;
set: const real EA = 2.84191e6;
set: const real GAY = 6.40131e5;
set: const real GAZ = 9.03881e5;
set: const real GJ = 3.10338;
set: const real EJY = 36.2794;
set: const real EJZ = 2.42873;
''')
        (dest/'princeton.ref').write_text('reference: 1, reference, global, null, reference, global, eye, reference, global, null, reference, global, null;\n')
        (dest/'princeton.nod').write_text('\n'.join(f'structural: {i}, static, reference, 1, {i*.508/count:.17g}, 0., 0., reference, 1, eye, reference, 1, null, reference, 1, null;' for i in range(1,count+1)))
        (dest/'princeton.elm').write_text('constitutive law: 1, 6, linear elastic generic, diag, EA, GAY, GAZ, GJ, EJY, EJZ;\n'+'\n'.join(f'beam3: {j+1}, {2*j}, reference, node, null, {2*j+1}, reference, node, null, {2*j+2}, reference, node, null, reference, 1, eye, reference, 1, same, same;' for j in range(beams)))
    (dest/'model').write_text(s)


def mb_result(path,case,step):
    import numpy as np
    spec=CASES[case]
    rows=[]
    for line in (path/'run.out').read_text().splitlines():
        fields=line.split()
        if fields and fields[0]=='Step' and fields[-1]=='1': rows.append(float(fields[2]))
    times=np.asarray(rows)
    data=np.loadtxt(path/'run.mov',ndmin=2)
    labels=[int(step)] if 'formulation' in spec else spec['labels']
    if case=='andrews':
        angles=np.unwrap(np.radians(data[data[:,0]==1600,6])) # Euler rz, degrees
        values=(angles-angles[0])[:,None]
    else:
        series=[data[data[:,0]==label,1:4] for label in labels]
        if any(len(x)!=len(times) for x in series): raise ValueError('MBDyn output/time mismatch')
        values=np.concatenate(series,axis=1)
    if len(times)==0 or len(values)!=len(times) or times[-1]<spec['T']-1e-10 or not np.all(np.diff(times)>0):
        raise ValueError('MBDyn incomplete trajectory')
    if 'formulation' in spec: grid=np.array([spec['T']])
    else: grid=np.linspace(0,spec['T'],round(spec['T']/spec['sample'])+1)
    vals=np.array([np.interp(grid,times,values[:,i]) for i in range(values.shape[1])]).T
    return dict(times=grid.tolist(),values=vals.tolist(),last_time=float(times[-1]))


def run(args,case,step,engine,rep):
    path=args.sortie/f'{case}-{step:.17g}-{engine}-{rep}'; path.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy(); env.update({key:'1' for key in ('RAYON_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
    if engine=='mbdyn':
        deck(args.benchmarks,case,step,path)
        if case=='andrews' and getattr(args,'initialisation_andrews',False):
            model=path/'model'
            model.write_text(replace(model.read_text(),'derivatives tolerance',1e-6))
        command=[str(args.mbdyn),'-f',str(path/'model'),'-o',str(path/'run')]
    else: command=[sys.executable,str(Path(__file__).resolve()),'--worker',case,str(step),str(path/'result.json')]
    command=['/usr/bin/time','-f','%M','-o',str(path/'rss.txt'),*command]
    start=time.perf_counter(); timed_out=False
    try:
        # Kill the whole process group on timeout: /usr/bin/time has a child.
        # Keep stdout even on failure and never overlap a timed-out solver
        # with the next timed measurement.
        p=subprocess.Popen(command,cwd=path,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            stdout,_=p.communicate(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            timed_out=True
            os.killpg(p.pid,signal.SIGKILL)
            stdout,_=p.communicate()
        wall=time.perf_counter()-start
        (path/'stdout.txt').write_bytes(stdout)
        if timed_out: raise TimeoutError(f'délai de {args.timeout:g} s dépassé ; groupe de processus arrêté')
        if p.returncode: raise RuntimeError(f'exit {p.returncode}: '+(path/'stdout.txt').read_text()[-1200:])
        result=mb_result(path,case,step) if engine=='mbdyn' else json.loads((path/'result.json').read_text())
        result.update(ok=True,wall_seconds=wall,rss_kib=int((path/'rss.txt').read_text().strip()), directory=str(path))
    except (Exception,) as e:
        result=dict(ok=False,error=str(e),wall_seconds=time.perf_counter()-start,directory=str(path))
        if engine=='vinkulum' and (path/'result.json').exists():
            result['diagnostic']=json.loads((path/'result.json').read_text())
    result['command']=command
    result['timed_out']=timed_out
    result['artifact_sha256']={p.name:sha(p) for p in sorted(path.iterdir()) if p.is_file()}
    (path/'measurement.json').write_text(json.dumps(result,indent=2))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mbdyn',type=Path,required=True); p.add_argument('--benchmarks',type=Path,required=True)
    p.add_argument('--sortie',type=Path,required=True); p.add_argument('--repetitions',type=int,default=3)
    p.add_argument('--cases',nargs='+',choices=CASES,default=list(CASES)); p.add_argument('--timeout',type=float,default=180)
    p.add_argument('--initialisation-andrews',action='store_true',help='Tolérance des dérivées initiales MBDyn : 1e-6 au lieu de 1e6')
    p.add_argument('--pas-supplementaires',type=json.loads,default={},help='Objet JSON cas -> liste de pas (ou intervalles Princeton)')
    p.add_argument('--pas',type=json.loads,default={},help='Remplace la grille candidate : objet JSON cas -> liste de pas')
    p.add_argument('--references',type=json.loads,default={},help='Objet JSON cas -> deux raffinements, ou cas -> {moteur: deux raffinements}')
    args=p.parse_args();
    if args.repetitions<1 or not 0<args.timeout<float('inf'):
        p.error('répétitions et délai doivent être strictement positifs')
    for case,steps in args.pas.items():
        if case not in CASES or not isinstance(steps,list) or not steps:
            p.error('grille candidate invalide')
        CASES[case]['steps']=steps
    for case,steps in args.references.items():
        grids=steps if isinstance(steps,dict) else {e:steps for e in ('vinkulum','mbdyn')}
        if case not in CASES or set(grids)!={'vinkulum','mbdyn'} or any(not isinstance(hs,list) or len(hs)!=2 for hs in grids.values()):
            p.error('références invalides')
        CASES[case]['ref_by_engine']=grids
    for case,steps in args.pas_supplementaires.items():
        if case not in CASES or any(not isinstance(h,(int,float)) or not 0<h<float("inf") for h in steps):
            p.error('cas inconnu ou pas invalide')
        CASES[case]['steps']=list(dict.fromkeys(CASES[case]['steps']+steps))
    for case in args.cases:
        spec=CASES[case]
        refs=spec.get('ref_by_engine',{e:spec['ref'] for e in ('vinkulum','mbdyn')})
        for step in spec['steps']+[h for hs in refs.values() for h in hs]:
            if isinstance(step,bool) or not isinstance(step,(int,float)) or not 0<step<float('inf'):
                p.error(f'{case}: pas invalide')
            if 'formulation' in spec and (step!=int(step) or step<2 or int(step)%2):
                p.error(f'{case}: le maillage exige un nombre pair d’intervalles, au moins 2')
        for coarse,fine in refs.values():
            if (fine<=coarse if 'formulation' in spec else fine>=coarse):
                p.error(f'{case}: les références doivent être ordonnées du grossier au fin')
    args.sortie=args.sortie.resolve(); args.sortie.mkdir(parents=True,exist_ok=False)
    args.mbdyn=args.mbdyn.resolve(); args.benchmarks=args.benchmarks.resolve()
    import vinkulum
    from vinkulum import _vinkulum
    report=dict(protocol=__doc__,cases={},metadata=dict(protocol_version=2,platform=platform.platform(),python=sys.version,
       created_utc=datetime.now(timezone.utc).isoformat(),command=[sys.executable,*sys.argv],
       repetitions=args.repetitions,timeout_seconds=args.timeout,affinity=sorted(os.sched_getaffinity(0)),
       vinkulum_version=vinkulum.__version__,vinkulum_package=str(Path(vinkulum.__file__).parent),
       initialisation_andrews=args.initialisation_andrews, pas=args.pas,pas_supplementaires=args.pas_supplementaires, references=args.references,
       cpu=Path('/proc/cpuinfo').read_text().split('model name')[1].splitlines()[0],
       harness_checkout_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
       harness_worktree_status=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip(),
       binary_sha256=sha(args.mbdyn),extension_sha256=sha(_vinkulum.__file__),
       script_sha256=sha(__file__),models_sha256=sha(ROOT/'ci/modeles_confrontation.py'),
       andrews_model_sha256=sha(Path(vinkulum.__file__).parent/'andrews.py'),
       input_includes_sha256={c:{f.name:sha(f) for f in sorted((args.benchmarks/CASES[c]['deck']).parent.iterdir()) if f.suffix in ('.set','.ref')} for c in args.cases},
       mbdyn_version=subprocess.check_output([str(args.mbdyn),'--version'],text=True,stderr=subprocess.STDOUT).strip(),
       mbdyn_checkout_commit=subprocess.check_output(['git','-C',str(args.benchmarks),'rev-parse','HEAD'],text=True).strip(),
       input_sha256={c:sha(args.benchmarks/CASES[c]['deck']) for c in args.cases}))
    for case in args.cases:
        results=[]; report['cases'][case]=dict(spec=CASES[case],runs=results)
        # Alternate engines; warm up every configuration, exclude from statistics.
        spec=CASES[case]
        refs=spec.get('ref_by_engine',{e:spec['ref'] for e in ('vinkulum','mbdyn')})
        for step in dict.fromkeys(spec['steps']+[h for hs in refs.values() for h in hs]):
            for rep in range(-1,args.repetitions if step in CASES[case]['steps'] else 1):
                for engine in (['vinkulum','mbdyn'] if rep%2 else ['mbdyn','vinkulum']):
                    if step not in spec['steps'] and step not in refs[engine]:
                        continue
                    result=run(args,case,step,engine,rep)
                    results.append(dict(step=step,engine=engine,rep=rep,**result))
                    print(case,step,engine,rep,round(result['wall_seconds'],3),result.get('error','OK'),flush=True)
                    (args.sortie/'resultats.json').write_text(json.dumps(report,indent=2))
    print(args.sortie/'resultats.json')


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--worker': worker(sys.argv[2],float(sys.argv[3]),sys.argv[4])
    else: main()
