"""Mesures isolées et alternées de binaires, avec états et échecs conservés."""
import argparse
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from bilan_contraintes import sha

ROOT = Path(__file__).resolve().parents[1]
WORKERS = ['ci/diagnostic_cascade.py', 'ci/diagnostic_newton_redondant.py',
           'ci/diagnostic_redondances.py', 'ci/diagnostic_contraintes.py',
           'ci/diagnostic_chaine_articulee.py', 'ci/confronte_mbdyn.py']
FAMILLES = {'cascade', 'double', 'mixte', 'parallelogrammes', 'boucle', 'spatiale',
            'chaine', 'rotors', 'articulee', 'princeton_milieu', 'princeton_integree'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='append', required=True, metavar='NOM=PYTHONPATH')
    p.add_argument('--cas', action='append', required=True, metavar='FAMILLE:TAILLE')
    p.add_argument('--repetitions', type=int, default=3)
    p.add_argument('--delai', type=float, default=90.)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    if args.repetitions < 1 or args.delai <= 0:
        p.error('répétitions et délai doivent être positifs')
    cases = []
    for case in args.cas:
        family, size = case.split(':')
        if family not in FAMILLES or int(size) <= 0:
            p.error(f'cas invalide : {case}')
        cases.append((family, int(size)))
    cpu = min(os.sched_getaffinity(0))
    base = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    envs, binaries = {}, {}
    for spec in args.version:
        name, path = spec.split('=', 1)
        if not name or name in envs:
            p.error('chaque version doit porter un nom unique non vide')
        envs[name] = dict(base, PYTHONPATH=str(Path(path).resolve()))
        binary = subprocess.check_output([sys.executable, '-c',
            'from vinkulum import _vinkulum; print(_vinkulum.__file__)'],
            env=envs[name], text=True).strip()
        binaries[name] = dict(path=binary, sha256=sha(binary))
    out = dict(cpu=cpu, python=platform.python_version(), plateforme=platform.platform(),
               repetitions=args.repetitions, fils_demandes=1, delai_processus_s=args.delai,
               binaires=binaries, programme_sha256=sha(__file__),
               workers_sha256={f: sha(ROOT/f) for f in WORKERS}, cas=cases, mesures=[])
    args.sortie.parent.mkdir(parents=True, exist_ok=True)
    for rep in range(args.repetitions):
        versions = list(envs)
        if rep % 2:
            versions.reverse()
        for family, size in cases:
            worker = [WORKERS[0], str(size)] if family == 'cascade' else [WORKERS[1], family, str(size)]
            for version in versions:
                cmd = ['taskset', '-c', str(cpu), sys.executable, *worker]
                try:
                    c = subprocess.run(cmd, cwd=ROOT, env=envs[version], capture_output=True,
                                       text=True, timeout=args.delai)
                    d = json.loads(c.stdout) if c.returncode == 0 else dict(erreur=c.stderr, stdout=c.stdout)
                    d.update(code_sortie=c.returncode, stderr=c.stderr)
                except subprocess.TimeoutExpired as e:
                    d = dict(erreur=f'délai de {args.delai:g} s dépassé', code_sortie=None,
                             stdout=e.stdout.decode() if isinstance(e.stdout, bytes) else e.stdout,
                             stderr=e.stderr.decode() if isinstance(e.stderr, bytes) else e.stderr)
                d.update(famille=family, taille=size, version=version, repetition=rep)
                out['mesures'].append(d)
                args.sortie.write_text(json.dumps(out, ensure_ascii=False, separators=(',', ':'), allow_nan=False)+'\n')
                print(family, size, version, rep, d.get('secondes'), d['erreur'], flush=True)


if __name__ == '__main__':
    main()
