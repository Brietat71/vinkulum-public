"""Mesures alternées de la projection redondante, avec témoins hors du nouveau chemin."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from bilan_contraintes import sha

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ancien-pythonpath', required=True, type=Path)
    p.add_argument('--sortie', required=True, type=Path)
    a = p.parse_args()
    cpu = min(os.sched_getaffinity(0))
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    envs, binaries = {}, {}
    for version in ('avant', 'apres'):
        e = env.copy()
        if version == 'avant':
            e['PYTHONPATH'] = str(a.ancien_pythonpath.resolve())
        else:
            e.pop('PYTHONPATH', None)
        path = subprocess.check_output([sys.executable, '-c',
               'from vinkulum import _vinkulum; print(_vinkulum.__file__)'], env=e, text=True).strip()
        envs[version] = e
        binaries[version] = dict(path=path, sha256=sha(path))
    assert binaries['avant']['sha256'] != binaries['apres']['sha256']
    workers = ['ci/diagnostic_redondances.py', 'ci/diagnostic_contraintes.py',
               'ci/diagnostic_chaine_articulee.py']
    cases = {'boucle': [16, 32, 64, 128, 256, 512], 'spatiale': [16, 32, 64, 128],
             'double': [8, 16, 32, 64], 'chaine': [128], 'rotors': [128], 'articulee': [128]}
    out = dict(cpu=cpu, python=platform.python_version(), plateforme=platform.platform(),
               repetitions=3, fils_demandes=1, delai_processus_s=90, binaires=binaries,
               programme_sha256=sha(__file__), workers_sha256={f: sha(ROOT/f) for f in workers},
               cas=cases, grandes_tailles_actuelles=[1024, 2048], mesures=[])
    for rep in range(3):
        for family, sizes in dict(cases, grande=[1024, 2048]).items():
            for nb in sizes:
                versions = ['apres'] if family == 'grande' else (['avant', 'apres'] if rep%2 == 0 else ['apres', 'avant'])
                for version in versions:
                    if family == 'articulee':
                        command = [str(ROOT/workers[2]), '--worker', str(nb)]
                    else:
                        command = [str(ROOT/workers[0]), 'boucle' if family == 'grande' else family, str(nb)]
                    cmd = ['taskset', '-c', str(cpu), sys.executable, *command]
                    try:
                        c = subprocess.run(cmd, env=envs[version], capture_output=True, text=True, timeout=90)
                        r = json.loads(c.stdout) if c.returncode == 0 else dict(erreur=c.stderr)
                        r.update(code_sortie=c.returncode, stderr=c.stderr)
                    except subprocess.TimeoutExpired:
                        r = dict(erreur='délai de 90 s dépassé')
                    r.update(famille=family, corps=nb, version=version, repetition=rep)
                    out['mesures'].append(r)
                    a.sortie.write_text(json.dumps(out, ensure_ascii=False, separators=(',', ':'), allow_nan=False)+'\n')
                    print(family, nb, version, rep, r.get('secondes'), r['erreur'], flush=True)


if __name__ == '__main__':
    main()
