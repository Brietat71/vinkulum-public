"""Contrôle croisé des temps statiques : versions alternées, même CPU Linux."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ancien-pythonpath', required=True, type=Path)
    p.add_argument('--sortie', required=True, type=Path)
    p.add_argument('--tailles', nargs='+', type=int, default=[240, 480])
    args = p.parse_args()
    cpu = min(os.sched_getaffinity(0))
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    probe = 'import vinkulum._vinkulum as m; print(m.__file__)'
    environnements = {}
    binaries = {}
    for version in ('avant', 'apres'):
        e = env.copy()
        if version == 'avant':
            e['PYTHONPATH'] = str(args.ancien_pythonpath.resolve())
        else:
            e.pop('PYTHONPATH', None)
        path = Path(subprocess.check_output([sys.executable, '-c', probe], env=e, text=True).strip())
        binaries[version] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        environnements[version] = e
    worker = Path(__file__).with_name('diagnostic_assemblage_local.py').resolve()
    out = dict(cpu=cpu, binaires=binaries, mesures=[], bilan=[],
               programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               worker_sha256=hashlib.sha256(worker.read_bytes()).hexdigest())
    def save():
        args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n')
    for repetition in range(3):
        for ne in args.tailles:
            versions = ('avant', 'apres') if repetition % 2 == 0 else ('apres', 'avant')
            for version in versions:
                cmd = ['taskset', '-c', str(cpu), sys.executable, str(worker), '--worker-statique', str(ne)]
                result = subprocess.run(cmd, env=environnements[version], text=True,
                                        capture_output=True, timeout=120)
                if result.returncode:
                    out['erreur'] = dict(commande=cmd, version=version, stderr=result.stderr,
                                         code_sortie=result.returncode)
                    save()
                    raise RuntimeError(out['erreur'])
                r = json.loads(result.stdout)
                r.update(version=version, repetition=repetition)
                out['mesures'].append(r)
                save()
                if r['erreur'] is not None or not r['statique_rapports'] or not all(
                        s['strict'] and s['statut'] == 'tolerance' for s in r['statique_rapports']):
                    raise RuntimeError(r)
                print(version, ne, r['seconds'], r['rss_max_kib'], flush=True)
    for ne in args.tailles:
        groupes = [[r for r in out['mesures'] if r['intervalles'] == ne and r['version'] == v]
                   for v in ('avant', 'apres')]
        times = [statistics.median(r['seconds'] for r in g) for g in groupes]
        out['bilan'].append(dict(intervalles=ne, avant_s=times[0], apres_s=times[1],
                                 gain=times[0]/times[1],
                                 rss_mio=[statistics.median(r['rss_max_kib'] for r in g)/1024 for g in groupes]))
    save()


if __name__ == '__main__':
    main()
