"""Princeton : Vinkulum, MBDyn et Exudyn à précision commune.

Les exécutables concurrents sont utilisés par leurs interfaces publiques.
Chaque configuration a un échauffement conservé, puis trois processus frais.
Le juge porte sur le maximum des trois composantes de position au bout,
avec deux raffinements par variante et une marge de sensibilité aux tolérances.
Les résultats ne constituent pas une borne d'erreur mathématiquement certifiée.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = {
    'vinkulum_milieu': dict(engine='vinkulum', case='princeton', steps=[8,20,40,60,80,120], refs=[160,320]),
    'vinkulum_integree': dict(engine='vinkulum', case='princeton_integree', steps=[4,6,8,10,12,16,20,40], refs=[80,160]),
    'mbdyn': dict(engine='mbdyn', case='princeton', steps=[4,6,8,10,20,40], refs=[80,160]),
    'exudyn_newton': dict(engine='exudyn', steps=[8,20,40,50,60,80,120], refs=[160,320]),
    'exudyn_modifie': dict(engine='exudyn', steps=[8,20,40,50,60,80,120], refs=[160,320]),
    'exudyn_fast_newton': dict(engine='exudyn', steps=[8,20,40,50,60,80,120], refs=[160,320]),
    'exudyn_fast_modifie': dict(engine='exudyn', steps=[8,20,40,50,60,80,120], refs=[160,320]),
}
TARGET = 1e-5
THREADS = {k: '1' for k in ('RAYON_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def encode(obj):
    return json.dumps(obj, ensure_ascii=False, allow_nan=False, separators=(',', ':')).encode()


def exudyn_worker(nel, variant, output, tolerance):
    import math
    sys.exudynFast = variant.startswith('exudyn_fast_')
    import exudyn as exu
    import numpy as np
    from exudyn.itemInterface import (
        Beam3D, NodeRigidBodyEP, ObjectGround, MarkerBodyRigid, MarkerNodeRigid,
        GenericJoint, LoadForceVector,
    )
    sc = exu.SystemContainer()
    mbs = sc.AddSystem()
    section = exu.BeamSection()
    section.stiffnessMatrix = np.diag([2.84191e6,6.40131e5,9.03881e5,3.10338,36.2794,2.42873])
    # Positive, arbitrary inertia: these values do not enter this static problem.
    section.massPerLength = 1.
    section.inertia = np.eye(3)
    nodes = [mbs.AddNode(NodeRigidBodyEP(referenceCoordinates=[i*.508/nel,0,0,1,0,0,0])) for i in range(nel+1)]
    ground = mbs.AddObject(ObjectGround())
    mg = mbs.AddMarker(MarkerBodyRigid(bodyNumber=ground, localPosition=[0,0,0]))
    mr = mbs.AddMarker(MarkerNodeRigid(nodeNumber=nodes[0]))
    joint = mbs.AddObject(GenericJoint(markerNumbers=[mg,mr], constrainedAxes=[1]*6))
    for i in range(nel):
        mbs.AddObject(Beam3D(nodeNumbers=nodes[i:i+2], physicsLength=.508/nel, sectionData=section))
    mt = mbs.AddMarker(MarkerNodeRigid(nodeNumber=nodes[-1]))

    def force(mbs, t, vector):
        f = 8.896/math.sqrt(2)*.5*(1-math.cos(math.pi*t/.5))
        return [0,f,f]

    mbs.AddLoad(LoadForceVector(markerNumber=mt, loadVector=[0,0,0], loadVectorUserFunction=force))
    settings = exu.SimulationSettings()
    settings.linearSolverType = exu.LinearSolverType.EigenSparse
    settings.parallel.numberOfThreads = 1
    settings.solutionSettings.writeSolutionToFile = False
    settings.staticSolver.numberOfLoadSteps = 50
    settings.staticSolver.loadStepDuration = .5
    settings.staticSolver.useLoadFactor = False
    settings.staticSolver.adaptiveStep = False
    settings.staticSolver.verboseMode = 1
    settings.staticSolver.newton.absoluteTolerance = tolerance
    settings.staticSolver.newton.useModifiedNewton = variant.endswith('_modifie')
    solver = exu.MainSolverStatic()
    steps = []

    def post(mbs, t):
        if t > 0:
            steps.append(dict(t=t, converged=solver.conv.newtonConverged,
                residual=solver.conv.residual, last_residual=solver.conv.lastResidual,
                iterations=solver.it.newtonSteps,
                tip=mbs.GetNodeOutput(nodes[-1], exu.OutputVariableType.Position).tolist()))
        return True

    mbs.SetPostStepUserFunction(post)
    mbs.Assemble()
    start = time.perf_counter()
    ok = solver.SolveSystem(mbs, settings)
    elapsed = time.perf_counter()-start
    tip = mbs.GetNodeOutput(nodes[-1], exu.OutputVariableType.Position)
    root_position = mbs.GetNodeOutput(nodes[0], exu.OutputVariableType.Position)
    root_rotation = mbs.GetNodeOutput(nodes[0], exu.OutputVariableType.RotationMatrix).reshape(3,3)
    reaction = mbs.GetObjectOutput(joint, exu.OutputVariableType.ForceLocal)
    quaternions = [np.array(mbs.GetNodeOutput(n, exu.OutputVariableType.Coordinates))[3:]+[1,0,0,0] for n in nodes]
    result = dict(solver_success=ok, times=[.5], values=[tip.tolist()], solve_seconds=elapsed,
        exudyn_version=exu.config.Version(), static_diagnostics=steps,
        exudyn_module=exu.MainSolverStatic.__module__, exudyn_build=exu.config.Version(True),
        absolute_tolerance=tolerance, relative_tolerance=settings.staticSolver.newton.relativeTolerance,
        constraints_max=float(max(np.max(np.abs(root_position)),np.max(np.abs(root_rotation-np.eye(3))),
                                  max(abs(float(q@q)-1) for q in quaternions))),
        force_balance_max=float(np.max(np.abs(reaction+force(mbs,.5,None)))),
        settings=str(settings), convergence=str(solver.conv), iterations=str(solver.it),
        solver_output=str(solver.output), dimensions=dict(ode2=solver.GetODE2size(),ae=solver.GetAEsize()))
    Path(output).write_bytes(encode(result))
    if not ok or len(steps)!=50 or not all(s['converged'] for s in steps):
        raise RuntimeError('Exudyn : cinquante paliers convergés non obtenus')


def run_exudyn(args, variant, n, rep, tolerance=1e-5, role='candidate'):
    path = args.sortie/f'{variant}-{n}-{role}-{rep}'
    path.mkdir()
    command = ['/usr/bin/time','-f','%M','-o',str(path/'rss.txt'),str(args.exudyn_python),
               str(Path(__file__).resolve()),'--worker',str(n),variant,str(path/'result.json'),str(tolerance)]
    start = time.perf_counter()
    timed_out = False
    p = subprocess.Popen(command,cwd=path,env=dict(os.environ,**THREADS),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,start_new_session=True)
    try:
        stdout,_ = p.communicate(timeout=args.timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(p.pid,signal.SIGKILL)
        stdout,_ = p.communicate()
    elapsed = time.perf_counter()-start
    (path/'stdout.txt').write_bytes(stdout)
    result = json.loads((path/'result.json').read_bytes()) if (path/'result.json').exists() else {}
    result.update(ok=p.returncode==0 and not timed_out, timed_out=timed_out,
                  wall_seconds=elapsed, directory=str(path), command=command)
    if result['ok']:
        result['rss_kib'] = int((path/'rss.txt').read_text())
    else:
        result['error'] = f'exit {p.returncode}; timeout={timed_out}'
    result['artifact_sha256'] = {f.name:sha(f.read_bytes()) for f in sorted(path.iterdir()) if f.is_file()}
    return result


def validate_run(r):
    import numpy as np
    if not r['ok']:
        return
    if not np.isfinite(r['values']).all() or np.shape(r['values'])!=(1,3):
        raise ValueError('position finale invalide')
    if not r['wall_seconds']>0 or not r['rss_kib']>0:
        raise ValueError('coût invalide')
    if r['variant'].startswith('exudyn'):
        expected='exudyn.exudynCPPfast' if r['variant'].startswith('exudyn_fast_') else 'exudyn.exudynCPP'
        if 'exudyn_module' in r and r['exudyn_module']!=expected:
            raise ValueError('module Exudyn différent du réglage demandé')
        ds = r['static_diagnostics']
        if not r['solver_success'] or len(ds)!=50 or not np.allclose([d['t'] for d in ds],np.arange(1,51)*.01,rtol=0,atol=1e-12):
            raise ValueError('paliers Exudyn incomplets')
        if any(not d['converged'] or not np.isfinite(d['residual']) or d['residual']<0 for d in ds):
            raise ValueError('palier Exudyn non convergé')
        if not 0<=r['constraints_max']<=1e-8 or not 0<=r['force_balance_max']<=1e-6:
            raise ValueError('bilan mécanique Exudyn invalide')
    elif r['variant'].startswith('vinkulum'):
        ds=r['static_diagnostics']
        if len(ds)!=50 or [d['palier_charge'] for d in ds]!=list(range(1,51)) or any(
            d['statut']!='tolerance' or d['strict'] is not True or not d['tolerance_finale_atteinte']
            or d['paliers_stagnation']!=0 or d['tol']>1e-8
            or not 0<=d['residu_relatif']<=d['tol'] or not 0<=d['contraintes']<=d['tol'] for d in ds):
            raise ValueError('paliers Vinkulum stricts non vérifiés')


def analyse(report):
    from bilan_confrontation import difference
    runs = report['runs']
    variants = report['metadata']['variants']
    for r in runs:
        validate_run(r)
    identities=[(r['variant'],r['n'],r['role'],r['rep']) for r in runs]
    if len(set(identities))!=len(identities):
        raise ValueError('essai dupliqué')
    def select(v,n,role,rep=0):
        return next((r for r in runs if r['variant']==v and r['n']==n and r['role']==role and r['rep']==rep and r['ok']),None)
    refs={}; changes={}
    for v,spec in variants.items():
        pair=[select(v,n,'reference') for n in spec['refs']]
        if all(pair):
            refs[v]=pair[-1]; changes[v]=difference(*pair)
    cross=max((difference(a,b) for a in refs.values() for b in refs.values()),default=float('inf'))
    tolerance_changes={}
    exudyn_variants=[v for v,spec in variants.items() if spec['engine']=='exudyn']
    for v in exudyn_variants:
        for n in (60,160,320):
            base=select(v,n,'candidate' if n==60 else 'reference')
            probe=select(v,n,'tolerance')
            if base and probe:
                tolerance_changes[f'{v}-{n}']=difference(base,probe)
    qualified=(len(refs)==len(variants) and len(tolerance_changes)==3*len(exudyn_variants) and
               max(changes.values())<=.1*TARGET and cross<=.2*TARGET and max(tolerance_changes.values())<=.01*TARGET)
    margin=max(changes.values(),default=0)+max(tolerance_changes.values(),default=0)
    configs=[]
    for v,spec in variants.items():
        for n in spec['steps']:
            samples=[r for r in runs if r['variant']==v and r['n']==n and r['role']=='candidate' and r['rep']>=0]
            good=[r for r in samples if r['ok']]
            complete=sorted(r['rep'] for r in samples)==list(range(report['metadata']['repetitions']))
            c=dict(variant=v,n=n,attempts=len(samples),successes=len(good),repetitions_complete=complete,eligible=False)
            if good:
                envelope=max((difference(r,ref) for r in good for ref in refs.values()),default=float('inf'))
                c.update(median_seconds=statistics.median(r['wall_seconds'] for r in good),
                    wall_seconds=[r['wall_seconds'] for r in good],rss_kib=[r['rss_kib'] for r in good],
                    final_values=good[0]['values'][0],repeat_difference=max(difference(good[0],r) for r in good),
                    error_envelope=envelope,estimated_error_with_margin=envelope+margin,
                    eligible=qualified and complete and len(good)==len(samples) and envelope+margin<=TARGET)
                if all('solve_seconds' in r for r in good):
                    c['median_solve_seconds']=statistics.median(r['solve_seconds'] for r in good)
            configs.append(c)
    best={v:min((c for c in configs if c['variant']==v and c['eligible']),key=lambda c:c['median_seconds'],default=None) for v in variants}
    return dict(definition=__doc__,target=TARGET,reference_changes=changes,reference_cross_gap=cross,
                tolerance_changes=tolerance_changes,references_qualified=qualified,margin=margin,
                configurations=configs,best=best,total_runs=len(runs),warmups=sum(r['rep']<0 for r in runs),
                failures=[dict(variant=r['variant'],n=r['n'],rep=r['rep'],role=r['role'],error=r['error']) for r in runs if not r['ok']])


def archive(report, prefix):
    sources={name:(ROOT/name).read_text() for name in ('ci/confronte_exudyn.py','ci/confronte_mbdyn.py','ci/modeles_confrontation.py','ci/bilan_confrontation.py')}
    if {n:sha(s.encode()) for n,s in sources.items()}!=report['metadata']['scripts_sha256']:
        raise ValueError('sources modifiées pendant les mesures')
    contents={}; entries=[]
    for r in report['runs']:
        embedded={}; omitted={}
        for name,digest in r['artifact_sha256'].items():
            data=(Path(r['directory'])/name).read_bytes()
            if sha(data)!=digest:
                raise ValueError('artefact de calcul modifié')
            if name=='run.mov':
                omitted[name]=digest
            else:
                contents[digest]=data.decode(); embedded[name]=digest
        entries.append(dict(directory=r['directory'],embedded=embedded,omitted=omitted))
    artifacts={'.json':encode(analyse(report)), '-essais.json.gz':gzip.compress(encode(report),mtime=0),
        '-journaux.json.gz':gzip.compress(encode(dict(sources=sources,contents=contents,entries=entries)),mtime=0)}
    manifest=dict(files={Path(str(prefix)+s).name:dict(sha256=sha(b),bytes=len(b)) for s,b in artifacts.items()},
                  scripts_sha256=report['metadata']['scripts_sha256'])
    artifacts['-manifest.json']=encode(manifest)
    for suffix,data in artifacts.items():
        path=Path(str(prefix)+suffix)
        with path.open('xb') as f:
            f.write(data)
    verify(prefix)


def verify(prefix):
    manifest=json.loads(Path(str(prefix)+'-manifest.json').read_bytes())
    for name,info in manifest['files'].items():
        data=(prefix.parent/name).read_bytes()
        if sha(data)!=info['sha256'] or len(data)!=info['bytes']:
            raise ValueError('empreinte d’archive invalide')
    report=json.loads(gzip.decompress(Path(str(prefix)+'-essais.json.gz').read_bytes()))
    if analyse(report)!=json.loads(Path(str(prefix)+'.json').read_bytes()):
        raise ValueError('bilan différent du recalcul')
    logs=json.loads(gzip.decompress(Path(str(prefix)+'-journaux.json.gz').read_bytes()))
    if {n:sha(s.encode()) for n,s in logs['sources'].items()}!=manifest['scripts_sha256'] or manifest['scripts_sha256']!=report['metadata']['scripts_sha256']:
        raise ValueError('sources non liées aux mesures')
    if any(sha(s.encode())!=h for h,s in logs['contents'].items()) or len(logs['entries'])!=len(report['runs']):
        raise ValueError('journaux invalides')
    for r,e in zip(report['runs'],logs['entries']):
        if e['directory']!=r['directory'] or dict(e['embedded'],**e['omitted'])!=r['artifact_sha256'] or any(h not in logs['contents'] for h in e['embedded'].values()):
            raise ValueError('journal non lié à son essai')
    print(f"{len(report['runs'])} essais : bilan, paliers, journaux et empreintes vérifiés.")


def main():
    if len(sys.argv)>1 and sys.argv[1]=='--worker':
        exudyn_worker(int(sys.argv[2]),sys.argv[3],sys.argv[4],float(sys.argv[5])); return
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--verifier',type=Path)
    p.add_argument('--sortie',type=Path)
    p.add_argument('--archive',type=Path)
    p.add_argument('--exudyn-python',type=Path)
    p.add_argument('--exudyn-wheel',type=Path)
    p.add_argument('--mbdyn',type=Path)
    p.add_argument('--benchmarks',type=Path)
    p.add_argument('--repetitions',type=int,default=3)
    p.add_argument('--timeout',type=float,default=60)
    args=p.parse_args()
    if args.verifier:
        verify(args.verifier); return
    if not all((args.sortie,args.archive,args.exudyn_python,args.exudyn_wheel,args.mbdyn,args.benchmarks)) or args.repetitions<3 or not 0<args.timeout<float('inf'):
        p.error('sortie, archive, exécutables, roue et entrées requis ; au moins trois répétitions')
    import numpy as np
    import scipy
    import vinkulum
    from confronte_mbdyn import run
    args.sortie=args.sortie.resolve(); args.sortie.mkdir(exist_ok=False)
    identity_code="import exudyn, numpy, scipy, sys, json, pathlib, hashlib; print(json.dumps(dict(version=exudyn.config.Version(),python=sys.version,numpy=numpy.__version__,scipy=scipy.__version__,extensions={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pathlib.Path(exudyn.__file__).parent.glob('*.so')})))"
    identity=json.loads(subprocess.check_output([str(args.exudyn_python),'-c',identity_code],env=dict(os.environ,**THREADS),text=True))
    release=json.loads((ROOT/f'docs/bancs/version-{vinkulum.__version__}.json').read_text())
    extensions=list(Path(vinkulum.__file__).parent.glob('*.so'))
    if len(extensions)!=1 or sha(extensions[0].read_bytes())!=release['roue']['extension_sha256']:
        raise ValueError('la roue Vinkulum ne correspond pas à la livraison')
    metadata=dict(created_utc=datetime.now(timezone.utc).isoformat(),command=sys.argv,variants=VARIANTS,
        repetitions=args.repetitions,timeout=args.timeout,affinity=sorted(os.sched_getaffinity(0)),
        platform=platform.platform(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,threads=THREADS,
        vinkulum_version=vinkulum.__version__,vinkulum_extension_sha256=sha(extensions[0].read_bytes()),vinkulum_wheel=release['roue'],
        exudyn=identity,exudyn_wheel=dict(name=args.exudyn_wheel.name,sha256=sha(args.exudyn_wheel.read_bytes())),
        mbdyn_binary_sha256=sha(args.mbdyn.read_bytes()),
        cpu=next(s for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')),
        harness_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scripts_sha256={n:sha((ROOT/n).read_bytes()) for n in ('ci/confronte_exudyn.py','ci/confronte_mbdyn.py','ci/modeles_confrontation.py','ci/bilan_confrontation.py')})
    report=dict(metadata=metadata,runs=[])
    # Each variant's grids are explicit; process order reverses every repetition.
    configs=[(v,n,'candidate') for v,s in VARIANTS.items() for n in s['steps']]
    for rep in range(-1,args.repetitions):
        ordered=configs[rep%len(configs):]+configs[:rep%len(configs)]
        if rep%2: ordered=list(reversed(ordered))
        for v,n,role in ordered:
            spec=VARIANTS[v]
            r=run_exudyn(args,v,n,rep) if spec['engine']=='exudyn' else run(args,spec['case'],n,spec['engine'],rep)
            r.update(variant=v,n=n,role=role,rep=rep); report['runs'].append(r)
            print(v,n,rep,r['ok'],round(r['wall_seconds'],6),flush=True)
    for v,spec in VARIANTS.items():
        for n in spec['refs']:
            for rep in (-1,0):
                r=run_exudyn(args,v,n,rep,role='reference') if spec['engine']=='exudyn' else run(args,spec['case'],n,spec['engine'],rep)
                r.update(variant=v,n=n,role='reference',rep=rep); report['runs'].append(r)
                print(v,n,'reference',rep,r['ok'],flush=True)
    for v,spec in VARIANTS.items():
        if spec['engine']!='exudyn':
            continue
        for n in (60,160,320):
            r=run_exudyn(args,v,n,0,tolerance=3e-6,role='tolerance')
            r.update(variant=v,n=n,role='tolerance',rep=0); report['runs'].append(r)
            print(v,n,'tolerance',r['ok'],flush=True)
    (args.sortie/'resultats.json').write_bytes(encode(report))
    archive(report,args.archive)


if __name__=='__main__':
    main()
