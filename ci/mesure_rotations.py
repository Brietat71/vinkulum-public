"""Compare les trajectoires énergie–moment dans des processus isolés alternés."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from bilan_contraintes import sha

ROOT = Path(__file__).resolve().parents[1]
WORKER = ROOT/'ci/diagnostic_rotation.py'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='append', required=True, metavar='NOM=PYTHONPATH')
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    cpu = min(os.sched_getaffinity(0))
    envs, binaries = {}, {}
    for spec in args.version:
        name, path = spec.split('=', 1)
        if not name or name in envs:
            p.error('chaque version doit porter un nom unique non vide')
        envs[name] = dict(os.environ, PYTHONPATH=str(Path(path).resolve()),
                         RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        binary = subprocess.check_output([sys.executable, '-c',
            'from vinkulum import _vinkulum; print(_vinkulum.__file__)'], env=envs[name], text=True).strip()
        binaries[name] = dict(path=binary, sha256=sha(binary))
    out = dict(cpu=cpu, python=platform.python_version(), plateforme=platform.platform(),
               repetitions=3, fils_demandes=1, delai_processus_s=60,
               binaires=binaries, programme_sha256=sha(__file__), worker_sha256=sha(WORKER), mesures=[])
    args.sortie.parent.mkdir(parents=True, exist_ok=True)
    for rep in range(3):
        for name in list(envs)[::1 if rep % 2 == 0 else -1]:
            try:
                c = subprocess.run(['taskset', '-c', str(cpu), sys.executable, str(WORKER)],
                                   cwd=ROOT, env=envs[name], capture_output=True, text=True, timeout=60)
                d = json.loads(c.stdout) if c.returncode == 0 else dict(erreur=c.stderr, stdout=c.stdout)
                d.update(code_sortie=c.returncode, stderr=c.stderr)
            except subprocess.TimeoutExpired as e:
                d = dict(erreur='délai de 60 s dépassé', code_sortie=None,
                         stdout=e.stdout.decode() if isinstance(e.stdout, bytes) else e.stdout,
                         stderr=e.stderr.decode() if isinstance(e.stderr, bytes) else e.stderr)
            d.update(version=name, repetition=rep)
            out['mesures'].append(d)
            args.sortie.write_text(json.dumps(out, ensure_ascii=False, separators=(',', ':'), allow_nan=False)+'\n')
            print(name, rep, d['code_sortie'], flush=True)


if __name__ == '__main__':
    main()
