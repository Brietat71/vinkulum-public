"""Comparaison alternée : Newton redondant mobile, témoins et échecs conservés."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from bilan_contraintes import sha

ROOT = Path(__file__).resolve().parents[1]
CAS = dict(double=[8, 16, 32, 64], mixte=[16, 32, 64], parallelogrammes=[1, 4, 8],
           articulee=[8, 128], spatiale=[128], chaine=[128], rotors=[128], boucle=[512],
           princeton_milieu=[10, 20], princeton_integree=[10, 20])
GRANDS = dict(double=[128, 256, 512, 1024], mixte=[128, 256, 512])
WORKERS = ['ci/diagnostic_newton_redondant.py', 'ci/diagnostic_redondances.py',
           'ci/diagnostic_contraintes.py', 'ci/diagnostic_chaine_articulee.py',
           'ci/confronte_mbdyn.py']


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ancien-pythonpath', required=True, type=Path)
    p.add_argument('--sortie', required=True, type=Path)
    args = p.parse_args()
    cpu = min(os.sched_getaffinity(0))
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    envs, binaries = {}, {}
    for version in ('avant', 'apres'):
        e = env.copy()
        if version == 'avant':
            e['PYTHONPATH'] = str(args.ancien_pythonpath.resolve())
        else:
            e.pop('PYTHONPATH', None)
        path = subprocess.check_output([sys.executable, '-c',
            'from vinkulum import _vinkulum; print(_vinkulum.__file__)'], env=e, text=True).strip()
        envs[version] = e
        binaries[version] = dict(path=path, sha256=sha(path))
    assert binaries['avant']['sha256'] != binaries['apres']['sha256']
    out = dict(cpu=cpu, python=platform.python_version(), plateforme=platform.platform(),
               repetitions=3, fils_demandes=1, delai_processus_s=90, binaires=binaries,
               programme_sha256=sha(__file__), workers_sha256={f: sha(ROOT/f) for f in WORKERS},
               cas=CAS, grandes_tailles_actuelles=GRANDS, mesures=[])
    jobs = [(f, n, False) for f, sizes in CAS.items() for n in sizes]
    jobs += [(f, n, True) for f, sizes in GRANDS.items() for n in sizes]
    for rep in range(3):
        for family, size, large in jobs:
            versions = ['apres'] if large else (['avant', 'apres'] if rep % 2 == 0 else ['apres', 'avant'])
            for version in versions:
                cmd = ['taskset', '-c', str(cpu), sys.executable, str(ROOT/WORKERS[0]), family, str(size)]
                try:
                    c = subprocess.run(cmd, env=envs[version], capture_output=True, text=True, timeout=90)
                    d = json.loads(c.stdout) if c.returncode == 0 else dict(erreur=c.stderr)
                    d.update(code_sortie=c.returncode, stderr=c.stderr)
                except subprocess.TimeoutExpired as e:
                    d = dict(erreur='délai de 90 s dépassé', code_sortie=None,
                             stdout=e.stdout.decode() if isinstance(e.stdout, bytes) else e.stdout,
                             stderr=e.stderr.decode() if isinstance(e.stderr, bytes) else e.stderr)
                d.update(famille=family, taille=size, version=version, repetition=rep,
                         grande_taille_sans_temoin=large)
                out['mesures'].append(d)
                args.sortie.write_text(json.dumps(out, ensure_ascii=False, separators=(',', ':'), allow_nan=False)+'\n')
                print(family, size, version, rep, d.get('secondes'), d['erreur'], flush=True)


if __name__ == '__main__':
    main()
