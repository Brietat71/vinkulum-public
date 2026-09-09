"""Compare deux extensions sur quatre topologies, dans des processus isolés."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ancien-pythonpath', required=True, type=Path)
    p.add_argument('--sortie', required=True, type=Path)
    p.add_argument('--repetitions', type=int, default=3)
    args = p.parse_args()
    cpu = min(os.sched_getaffinity(0))
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    probe = 'import vinkulum._vinkulum as m; print(m.__file__)'
    environments, binaries = {}, {}
    for version in ('avant', 'apres'):
        e = env.copy()
        if version == 'avant':
            e['PYTHONPATH'] = str(args.ancien_pythonpath.resolve())
        else:
            e.pop('PYTHONPATH', None)
        path = Path(subprocess.check_output([sys.executable, '-c', probe], env=e, text=True).strip())
        binaries[version] = dict(path=str(path), sha256=sha(path))
        environments[version] = e
    if binaries['avant']['sha256'] == binaries['apres']['sha256']:
        raise RuntimeError('Les deux extensions sont identiques')
    workers = ['ci/diagnostic_contraintes.py', 'ci/diagnostic_chaine_articulee.py']
    out = dict(cpu=cpu, python=platform.python_version(), plateforme=platform.platform(),
               fils_demandes=1, delai_par_processus_s=90, repetitions=args.repetitions,
               binaires=binaries, programme_sha256=sha(__file__),
               workers_sha256={f: sha(ROOT/f) for f in workers}, mesures=[])
    families = dict(chaine=[32, 64, 128, 256, 512], rotors=[32, 64, 128, 256, 512],
                    boucle=[8, 16, 32, 64, 128, 256], articulee=[8, 16, 32, 64, 128, 256])
    for rep in range(args.repetitions):
        for family, sizes in families.items():
            for nb in sizes:
                if family == 'articulee':
                    worker = ROOT/workers[1]
                    arguments = ['--worker', str(nb)]
                else:
                    worker = ROOT/workers[0]
                    arguments = ['--worker', family, str(nb)]
                versions = ('avant', 'apres') if rep % 2 == 0 else ('apres', 'avant')
                for version in versions:
                    cmd = ['taskset', '-c', str(cpu), sys.executable, str(worker), *arguments]
                    try:
                        c = subprocess.run(cmd, env=environments[version], capture_output=True,
                                           text=True, timeout=out['delai_par_processus_s'])
                        r = json.loads(c.stdout) if c.returncode == 0 else dict(erreur=c.stderr)
                        r.update(code_sortie=c.returncode, stderr=c.stderr)
                    except subprocess.TimeoutExpired:
                        r = dict(erreur='délai de 90 s dépassé')
                    r.update(famille=family, corps=nb, repetition=rep, version=version)
                    out['mesures'].append(r)
                    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2,
                                                     allow_nan=False)+'\n')
                    print(version, family, nb, rep, r.get('secondes'), r['erreur'], flush=True)


if __name__ == '__main__':
    main()
