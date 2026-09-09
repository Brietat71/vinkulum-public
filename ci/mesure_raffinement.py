"""Statique : comparaison alternée, références mécaniques et repères tournés."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
WORKERS = ['ci/diagnostic_cascade.py', 'ci/diagnostic_newton_redondant.py',
           'ci/diagnostic_redondances.py', 'ci/diagnostic_contraintes.py',
           'ci/diagnostic_chaine_articulee.py', 'ci/confronte_mbdyn.py']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version', action='append', required=True, metavar='NOM=PYTHONPATH')
    p.add_argument('--sortie', type=Path, required=True)
    p.add_argument('--repetitions', type=int, default=3)
    p.add_argument('--cpu', type=int, default=8)
    args = p.parse_args()
    if args.cpu not in os.sched_getaffinity(0) or args.repetitions < 1:
        p.error('CPU disponible et nombre de répétitions positif requis')
    environments, versions = {}, {}
    probe = ('import json, sys, numpy, scipy; from vinkulum import _vinkulum, __version__; '
             'print(json.dumps(dict(extension=_vinkulum.__file__, version=__version__, '
             'python=sys.version, numpy=numpy.__version__, scipy=scipy.__version__)))')
    for spec in args.version:
        name, path = spec.split('=', 1)
        if not name or name in environments:
            p.error('chaque version doit avoir un nom unique non vide')
        env = dict(os.environ, PYTHONPATH=str(Path(path).resolve()), RAYON_NUM_THREADS='1',
                   OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
        for key in ('VINKULUM_TRACE', 'VINKULUM_EXPORT_STATIQUE'):
            env.pop(key, None)
        info = json.loads(subprocess.check_output([sys.executable, '-c', probe],
                                                  env=env, text=True))
        info['extension_sha256'] = sha(info['extension'])
        environments[name], versions[name] = env, info
    if len({v['extension_sha256'] for v in versions.values()}) != len(versions):
        p.error('les extensions doivent être distinctes')
    cases = [(family, n) for family in ('cascade', 'cascade_tournee') for n in (8, 16, 32)]
    cases += [(family, n) for family in ('double', 'mixte', 'parallelogrammes') for n in (8, 32)]
    cases += [(family, n) for family in ('chaine', 'rotors', 'articulee', 'spatiale', 'boucle')
              for n in (32, 128)]
    cases += [(family, n) for family in ('princeton_milieu', 'princeton_integree') for n in (10, 40)]
    out = dict(cpu=args.cpu, plateforme=platform.platform(), fils_demandes=1,
               repetitions=args.repetitions, delai_processus_s=90,
               versions=versions, programme_sha256=sha(__file__),
               workers_sha256={f: sha(ROOT/f) for f in WORKERS},
               sources_candidat_sha256={str(f.relative_to(ROOT)): sha(f)
                                        for f in sorted((ROOT/'src').rglob('*.rs'))}, mesures=[])
    args.sortie.parent.mkdir(parents=True, exist_ok=True)
    for rep in range(args.repetitions):
        for family, size in cases:
            if family.startswith('cascade'):
                worker = WORKERS[0]
                arguments = [str(size)] + (['--repere-tourne'] if family.endswith('tournee') else [])
            else:
                worker = WORKERS[1]
                arguments = [family, str(size)]
            for name in list(environments)[::1 if rep % 2 == 0 else -1]:
                cmd = ['taskset', '-c', str(args.cpu), sys.executable, str(ROOT/worker), *arguments]
                try:
                    c = subprocess.run(cmd, cwd=ROOT, env=environments[name], text=True,
                                       capture_output=True, timeout=out['delai_processus_s'])
                    try:
                        d = json.loads(c.stdout)
                    except json.JSONDecodeError:
                        d = dict(erreur='sortie JSON absente ou invalide', stdout=c.stdout)
                    d.update(code_sortie=c.returncode, stderr=c.stderr)
                    if c.returncode and not d.get('erreur'):
                        d['erreur'] = f'processus terminé avec le code {c.returncode}'
                except subprocess.TimeoutExpired as e:
                    d = dict(erreur='délai de 90 s dépassé', code_sortie=None,
                             stdout=e.stdout.decode() if isinstance(e.stdout, bytes) else e.stdout,
                             stderr=e.stderr.decode() if isinstance(e.stderr, bytes) else e.stderr)
                d.update(famille=family, taille=size, version=name, repetition=rep)
                out['mesures'].append(d)
                args.sortie.write_text(json.dumps(out, ensure_ascii=False, separators=(',', ':'),
                                                 allow_nan=False)+'\n')
                print(name, family, size, rep, d.get('secondes'), d.get('erreur'), flush=True)


if __name__ == '__main__':
    main()
