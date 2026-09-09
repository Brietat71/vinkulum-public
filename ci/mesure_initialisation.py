"""Compare les initialisations dans des processus neufs, alternés et isolés."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--avant', required=True)
    p.add_argument('--apres', default=sys.executable)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    affinite = sorted(os.sched_getaffinity(0))
    cpu = affinite[0]
    os.sched_setaffinity(0, {cpu})
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', RAYON_NUM_THREADS='1')
    configurations = [(cas, nb, tourne) for cas in ('analyse', 'initialisation')
                      for nb in (1, 8, 16, 32, 48) for tourne in (False, True)]
    configurations += [('libre', nb, False) for nb in (64, 512, 2048)]
    records = []
    bilan = []
    out = dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), machine=platform.platform(),
               cpu=cpu, affinite_disponible=affinite, fils=1, repetitions=3, alternance='AB / BA / AB',
               sources_sha256={str(f.relative_to(ROOT)):sha(f) for f in sorted((ROOT/'src').rglob('*.rs'))},
               versions_sources_sha256={f:sha(ROOT/f) for f in ('Cargo.toml', 'Cargo.lock', 'pyproject.toml')},
               programme_sha256=sha(__file__), diagnostic_sha256=sha(ROOT/'ci/diagnostic_initialisation.py'),
               resultats=records, bilan=bilan)
    for cas, nb, tourne in configurations:
        essais = []
        for repetition in range(3):
            for label in (('avant', 'apres') if repetition % 2 == 0 else ('apres', 'avant')):
                cmd = [getattr(args, label), str(ROOT/'ci/diagnostic_initialisation.py'), cas,
                       '--corps' if cas == 'libre' else '--cellules', str(nb)]
                if tourne:
                    cmd.append('--tourne')
                proc = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)
                if proc.returncode:
                    raise RuntimeError(f'{cmd}: {proc.returncode}\n{proc.stdout}\n{proc.stderr}')
                record = json.loads(proc.stdout)
                record.update(version_comparaison=label, repetition=repetition+1)
                essais.append(record); records.append(record)
                args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
        ligne = dict(cas=cas, taille=nb, tourne=tourne)
        for label in ('avant', 'apres'):
            subset = [r for r in essais if r['version_comparaison'] == label]
            errors = [r['erreur'] for r in subset if r['erreur'] is not None]
            controles = [r['controles'] for r in subset]
            valide = not errors and all(all((v if isinstance(v, bool) else v == nb if k == 'mobilites' else v <= 5e-10)
                for k, v in c.items()) for c in controles)
            ligne[label] = dict(valide=valide, erreurs=errors, secondes_mediane=statistics.median(r['secondes'] for r in subset),
                                rss_kib_max=max(r['rss_kib'] for r in subset), controles=controles)
        ligne['rapport_temps_avant_sur_apres'] = (ligne['avant']['secondes_mediane']/ligne['apres']['secondes_mediane']
            if ligne['avant']['valide'] and ligne['apres']['valide'] else None)
        bilan.append(ligne)
        print(json.dumps(ligne, ensure_ascii=False), flush=True)
    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    if not all(l['apres']['valide'] for l in bilan):
        raise SystemExit('Des cas candidats échouent : résultats conservés.')


if __name__ == '__main__':
    main()
