"""Mesure le prototype Rust de poutre mixte contre le continuum de Cosserat.

Un processus à la fois, un échauffement et trois essais par configuration.
L'erreur porte sur le bout ET sur 961 points de la ligne moyenne. Les
échecs restent dans la campagne ; aucune modification des tolérances.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import statistics
import subprocess
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--sortie',type=Path,required=True)
    p.add_argument('--timeout',type=float,default=60.)
    args=p.parse_args()
    args.sortie=args.sortie.resolve(); args.sortie.mkdir(parents=True,exist_ok=False)
    args.executable=args.executable.resolve()
    payload=args.reference.read_bytes()
    reference=json.loads(gzip.decompress(payload) if args.reference.suffix=='.gz' else payload)
    fine=reference['references'][-1]
    positions=np.asarray(fine['positions'])
    if positions.shape!=(961,3) or reference['ecart_raffinement_m']>1e-9:
        raise ValueError('référence continue invalide')
    sources=[Path('examples/poutre_mixte.rs'),*sorted(Path('examples/poutre_mixte').glob('*.rs')),
             Path('src/ad.rs'),Path('ci/mesure_poutre_mixte.py'),Path('ci/reference_cosserat.py')]
    report=dict(definition=__doc__,metadata=dict(executable_sha256=sha(args.executable),
        reference_sha256=sha(args.reference),sources_sha256={str(path):sha(ROOT/path) for path in sources},
        cpu=Path('/proc/cpuinfo').read_text().split('model name')[1].splitlines()[0],
        affinity=sorted(os.sched_getaffinity(0)),python=platform.python_version(),
        repetitions=3,timeout_s=args.timeout),essais=[],configurations=[])
    env=dict(os.environ,RAYON_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for order,counts in [(1,[2,4,8,16]),(2,[1,2,3,4,8])]:
        samples=[]
        for count in counts:
            samples=[]
            for rep in range(-1,3):
                directory=args.sortie/f'ordre-{order}-elements-{count}-essai-{rep}'
                directory.mkdir()
                command=['/usr/bin/time','-f','%M','-o',str(directory/'rss.txt'),str(args.executable),
                         str(order),str(count),'princeton','50','profil']
                start=time.perf_counter()
                process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=env,start_new_session=True)
                timed_out=False
                try:
                    stdout,stderr=process.communicate(timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    timed_out=True; os.killpg(process.pid,signal.SIGKILL)
                    stdout,stderr=process.communicate()
                wall=time.perf_counter()-start
                (directory/'stdout.json').write_bytes(stdout); (directory/'stderr.txt').write_bytes(stderr)
                result=dict(ordre=order,elements=count,repetition=rep,temps_total_s=wall,ok=False,
                    delai_depasse=timed_out,code=process.returncode,command=command,stderr=stderr.decode(),
                    stdout_sha256=sha(directory/'stdout.json'),stderr_sha256=sha(directory/'stderr.txt'))
                if process.returncode==0 and not timed_out:
                    numerical=json.loads(stdout)
                    line=np.asarray(numerical['positions'])
                    if line.shape!=positions.shape or not np.isfinite(line).all():
                        raise ValueError('ligne moyenne invalide')
                    if numerical['residu_max_paliers']>=1e-8 or numerical['residu_interne']>=1e-9:
                        raise ValueError('équilibre déclaré sans résidu suffisant')
                    if not np.allclose(line[-1],numerical['position'],rtol=0.,atol=1e-14):
                        raise ValueError('position au bout incohérente')
                    result.update(ok=True,resultat=numerical,
                        rss_kib=int((directory/'rss.txt').read_text().splitlines()[-1]),
                        erreur_bout_m=float(np.max(np.abs(line[-1]-positions[-1]))),
                        erreur_ligne_m=float(np.max(np.abs(line-positions))),
                        erreur_reaction_moment_nm=float(np.max(np.abs(np.asarray(numerical['reaction_racine'][3:])+fine['moment_racine']))))
                report['essais'].append(result)
                if rep>=0: samples.append(result)
                print(order,count,rep,round(wall,3),result.get('erreur_bout_m',result['stderr'][-220:]),flush=True)
                (args.sortie/'resultats.json').write_text(json.dumps(report,indent=2)+'\n')
            good=[s for s in samples if s['ok']]
            config=dict(ordre=order,elements=count,reussites=len(good),essais=len(samples))
            if good:
                config.update(temps_median_s=statistics.median(s['temps_total_s'] for s in good),
                    temps_calcul_median_s=statistics.median(s['resultat']['temps_s'] for s in good),
                    erreur_bout_m=max(s['erreur_bout_m'] for s in good),
                    erreur_ligne_m=max(s['erreur_ligne_m'] for s in good),
                    rss_max_kib=max(s['rss_kib'] for s in good),
                    ecart_repetitions_m=max(float(np.max(np.abs(np.asarray(s['resultat']['positions'])-good[0]['resultat']['positions']))) for s in good))
            report['configurations'].append(config)
            (args.sortie/'resultats.json').write_text(json.dumps(report,indent=2)+'\n')
    print(args.sortie/'resultats.json')


if __name__=='__main__': main()
