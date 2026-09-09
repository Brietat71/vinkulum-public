"""Compare deux binaires sur des équilibres avec contrôles indépendants.

Les versions alternent, sur un CPU et un fil par bibliothèque. Chaque
processus chauffe un modèle puis chronomètre la résolution d'un modèle neuf.
Construction, référence indépendante et import Python sont hors chronomètre.
Ce protocole ne compare pas les coûts d'usage complets de MBDyn ou Simpack.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def corpus():
    cases = []
    for form in ('milieu', 'integree'):
        base = dict(famille='poutre', formulation=form)
        cases += [base | dict(n=n, rampe=r) for n in (8, 20, 60) for r in (False, True)]
        cases += [base | dict(n=20, angle=a) for a in (0., 30., 60., 90.)]
        cases += [base | dict(n=10, tourne=True, unites=u) for u in
                  ([1., 1.], [1e-3, 1.], [1e3, 1.], [1., 1e-3], [1., 1e3])]
        cases += [base | dict(n=8, moment=True, tourne=q) for q in (False, True)]
    cases += [dict(famille='cable', n=n) for n in (40, 52, 54, 60, 80, 120)]
    cases += [dict(famille='cascade', n=n, tourne=q) for n in (1, 4, 16, 32) for q in (False, True)]
    cases += [dict(famille='contact', n=n) for n in (1, 32)]
    cases += [dict(famille='sans_equilibre', n=1)]
    return cases


def empreinte(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def environnement(package=None):
    env = os.environ.copy()
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                'RAYON_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[key] = '1'
    env.pop('VINKULUM_TRACE', None)
    env.pop('PYTHONPATH', None)
    if package:
        env['PYTHONPATH'] = str(package)
    return env


def identite(python, package):
    script = ('import json,platform,hashlib,pathlib,numpy,scipy,vinkulum; '
              'from vinkulum import _vinkulum as v; '
              'print(json.dumps(dict(python=platform.python_version(),numpy=numpy.__version__, '
              'scipy=scipy.__version__,version=vinkulum.__version__,extension=v.__file__, '
              'sha256=hashlib.sha256(pathlib.Path(v.__file__).read_bytes()).hexdigest())))')
    return json.loads(subprocess.check_output([str(python), '-c', script],
                      env=environnement(package), cwd='/tmp', text=True))


def controle(result):
    """Seuils fixés pour ce corpus ; aucune équivalence déduite du temps seul."""
    if 'worker_error' in result:
        return [result['worker_error']]
    issues = []
    if result['spec']['famille'] == 'sans_equilibre':
        if not result['erreur'] or result['restaure'] is not True:
            issues.append('Un problème sans équilibre doit échouer et restaurer son état.')
        return issues
    if result['erreur']:
        return [result['erreur']]
    for r in result['rapports']:
        if r['statut'] != 'tolerance' or r['paliers_stagnation']:
            issues.append('Tolérance stricte non atteinte à tous les paliers.')
        if r['residu_libre'] > r['tol'] * r['echelle_force'] or r['contraintes'] > r['tol']:
            issues.append('Résidu ou contrainte au-dessus de la tolérance publiée.')
    for key, value in result['controles'].items():
        # Forces/moments et géométrie sont contrôlés séparément. L'arc
        # analytique n'est demandé qu'à la formulation qui le reproduit.
        limit = 2e-6 if key.endswith(('_N', '_Nm')) or key == 'reactions_hors_plan' else 1e-7
        if not value <= limit:
            issues.append(f'{key}={value:.9g} > {limit}')
    return issues


def bilan(records):
    summary = []
    for case in records:
        row = dict(spec=case['spec'], versions={})
        for version, samples in case['versions'].items():
            failures = [dict(repetition=i, raisons=controle(s)) for i, s in enumerate(samples) if controle(s)]
            durations = [s['temps_s'] for s in samples if 'temps_s' in s]
            row['versions'][version] = dict(controles_reussis=not failures,
                echecs=failures, mediane_s=statistics.median(durations) if durations else None,
                evaluations=[sum(r['evaluations'] for r in s.get('rapports', [])) for s in samples],
                tentatives=[sum(r['tentatives'] for r in s.get('rapports', [])) for s in samples])
        a, b = (row['versions'].get(k) for k in ('reference', 'candidat'))
        row['ecarts_etats'] = []
        if a and b and a['controles_reussis'] and b['controles_reussis']:
            for old, new in zip(case['versions']['reference'], case['versions']['candidat'], strict=True):
                if old.get('erreur') or new.get('erreur'):
                    continue
                def ecart(x, y):
                    if isinstance(x, list):
                        return max((ecart(u, v) for u, v in zip(x, y, strict=True)), default=0.)
                    return abs(x-y)
                row['ecarts_etats'].append(dict(position=ecart(old['position'], new['position']),
                                                rotation=ecart(old['rotation'], new['rotation'])))
        row['etats_equivalents'] = all(e['position'] <= 1e-8 and e['rotation'] <= 1e-7
                                      for e in row['ecarts_etats']) if row['ecarts_etats'] else None
        if a and b and a['controles_reussis'] and b['controles_reussis'] and a['mediane_s'] and b['mediane_s']:
            row['temps_candidat_sur_reference'] = b['mediane_s']/a['mediane_s']
        summary.append(row)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--candidat', type=Path, required=True)
    p.add_argument('--paquet-candidat', type=Path)
    p.add_argument('--sortie', type=Path, required=True)
    p.add_argument('--repetitions', type=int, default=3)
    p.add_argument('--delai', type=float, default=90.)
    p.add_argument('--cpu', type=int, default=min(os.sched_getaffinity(0)))
    args = p.parse_args()
    if args.repetitions < 1 or args.delai <= 0:
        p.error('repetitions et delai doivent être positifs')
    os.sched_setaffinity(0, {args.cpu})
    # Résoudre le lien de bin/python ferait perdre l'environnement virtuel.
    versions = dict(reference=(args.reference.absolute(), None),
                    candidat=(args.candidat.absolute(), args.paquet_candidat.resolve() if args.paquet_candidat else None))
    sources = ['ci/diagnostic_globalisation.py', 'ci/diagnostic_cascade.py', 'ci/mesure_globalisation.py',
               'src/lib.rs', 'src/tangent.rs', 'src/ad.rs']
    output = dict(metadata=dict(date_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        platform=platform.platform(), cpu=args.cpu, repetitions=args.repetitions, delai_s=args.delai,
        sources={s:empreinte(ROOT/s) for s in sources},
        binaires={k:identite(*v) for k, v in versions.items()}), cas=[])
    args.sortie.parent.mkdir(parents=True, exist_ok=True)
    for index, spec in enumerate(corpus()):
        case = dict(spec=spec, versions={k:[] for k in versions});output['cas'].append(case)
        for repetition in range(args.repetitions):
            order = list(versions)
            if (index+repetition) % 2:
                order.reverse()
            for version in order:
                python, package = versions[version]
                cmd = [str(python), str(ROOT/'ci/diagnostic_globalisation.py'), json.dumps(spec), '--echauffement']
                try:
                    process = subprocess.run(cmd, env=environnement(package), cwd='/tmp',
                        capture_output=True, text=True, timeout=args.delai)
                    result = json.loads(process.stdout) if process.returncode == 0 else dict(worker_error=f'code {process.returncode}')
                    result['stderr'] = process.stderr
                except subprocess.TimeoutExpired as e:
                    result = dict(worker_error=f'delai {args.delai}s', stdout=str(e.stdout), stderr=str(e.stderr))
                case['versions'][version].append(result)
        row = bilan([case])[0]
        print(index+1, spec, {k:dict(temps=v['mediane_s'], ok=v['controles_reussis'], evaluations=v['evaluations'])
                             for k, v in row['versions'].items()}, flush=True)
        output['bilan'] = bilan(output['cas'])
        args.sortie.write_text(json.dumps(output, ensure_ascii=False, allow_nan=False, separators=(',', ':'))+'\n')
    failed = [r for r in output['bilan'] if not r['versions']['candidat']['controles_reussis'] or r['etats_equivalents'] is False]
    print('Cas candidats hors critères :', len(failed), '/', len(output['cas']))
    raise SystemExit(bool(failed))


if __name__ == '__main__':
    main()
