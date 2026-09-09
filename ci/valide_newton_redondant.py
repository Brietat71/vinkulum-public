"""Contrôles physiques, échecs du témoin, comparaison et empreintes du Newton réduit."""
import argparse
import datetime
import json
from pathlib import Path
import shutil
import statistics

import numpy as np
from vinkulum import _vinkulum
from bilan_contraintes import etat, sha, strict
from mesure_newton_redondant import CAS, GRANDS, WORKERS

ROOT = Path(__file__).resolve().parents[1]
BANCS = ROOT/'docs/bancs'
BASELINE = '32d333a80b62c87a2b05f99e02e93cc77a5e87a1820cd4cf2d5d1029590875e0'


def controles(sample):
    if 'controles' in sample:
        return sample['controles']
    return {k: v for k, v in sample.items() if (
        k.startswith('erreur_') or k == 'orthogonalite_autocontraintes') and v is not None}


def bilan(data):
    assert data['programme_sha256'] == sha(ROOT/'ci/mesure_newton_redondant.py')
    assert data['cas'] == CAS and data['grandes_tailles_actuelles'] == GRANDS
    assert data['repetitions'] == 3 and data['fils_demandes'] == 1
    assert data['delai_processus_s'] == 90
    assert data['workers_sha256'] == {f: sha(ROOT/f) for f in WORKERS}
    assert data['binaires']['avant']['sha256'] == BASELINE
    assert data['binaires']['apres']['sha256'] == sha(_vinkulum.__file__)
    expected = {(f, n, v, r) for f, sizes in CAS.items() for n in sizes
                for v in ('avant', 'apres') for r in range(3)}
    expected |= {(f, n, 'apres', r) for f, sizes in GRANDS.items() for n in sizes for r in range(3)}
    samples = {(s['famille'], s['taille'], s['version'], s['repetition']): s for s in data['mesures']}
    assert len(samples) == len(data['mesures']) == 141 and samples.keys() == expected
    failures = set()
    for key, s in samples.items():
        f, n, v, _ = key
        assert s['code_sortie'] == 0 and s['secondes'] > 0 and s['rss_kib'] > 0, key
        etat(s)
        if s['erreur'] is not None:
            assert f.startswith('princeton_') and v == 'avant', (key, s['erreur'])
            assert s['rapport']['statut'] == 'echec' and s['rapport']['strict']
            assert not s['rapport']['tolerance_finale_atteinte']
            assert s['rapport']['residu_relatif'] > s['rapport']['tol']
            assert s['rapport']['etat_restaure'] and s['etat'] == s['etat_initial']
            failures.add(key)
            continue
        strict(s['rapport'])
        assert s['reactions'] is not None
        assert all(np.isfinite(values).all() for _, values in s['reactions'])
        checks = controles(s)
        assert checks and all(np.isfinite(v) and v >= 0 for v in checks.values()), key
        for name, value in checks.items():
            threshold = 1e-12 if name in ('orthogonalite_autocontraintes', 'erreur_pose_max') else 1e-8
            if name == 'erreur_reactions_absolue':
                # Grandes réactions de moment de la chaîne soudée : publier
                # l'écart absolu et contrôler son erreur relative à la charge.
                threshold = 1e-9*max(9.81*n, 9.81*.2*n*(n-1)/2)
            assert value <= threshold, (key, name, value, threshold)
    expected_failures = {(f, n, 'avant', r) for f in ('princeton_milieu', 'princeton_integree')
                         for n in (10, 20) for r in range(3)}
    assert failures == expected_failures
    comparison = []
    for f, n in sorted({(f, n) for f, n, _, _ in expected}):
        row = dict(famille=f, taille=n, versions={})
        for v in ('avant', 'apres'):
            group = [s for (family, size, version, _), s in samples.items() if (family, size, version) == (f, n, v)]
            if not group:
                continue
            times = [s['secondes'] for s in group]
            checks = [controles(s) for s in group]
            row['versions'][v] = dict(mediane_s=statistics.median(times), etendue_s=[min(times), max(times)],
                rss_mio=statistics.median(s['rss_kib'] for s in group)/1024,
                statuts=[s['rapport']['statut'] for s in group], evaluations=[s['rapport']['evaluations'] for s in group],
                tentatives=[s['rapport']['tentatives'] for s in group],
                controles_max={k: max(c[k] for c in checks) for k in checks[0]})
        if len(row['versions']) == 2 and not f.startswith('princeton_'):
            row['gain'] = row['versions']['avant']['mediane_s']/row['versions']['apres']['mediane_s']
            errors, reaction_error = [0., 0.], 0.
            for i in range(3):
                for j in range(3):
                    a, b = samples[f, n, 'avant', i], samples[f, n, 'apres', j]
                    for k, (x, y) in enumerate(zip(etat(a), etat(b))):
                        errors[k] = max(errors[k], float(np.max(np.abs(x-y))))
                    for (na, ra), (nb, rb) in zip(a['reactions'], b['reactions'], strict=True):
                        assert na == nb
                        reaction_error = max(reaction_error, float(np.max(np.abs(np.asarray(ra)-rb))))
            assert max(errors) <= 1e-10, (f, n, errors)
            row.update(ecart_position_compensee_m=errors[0], ecart_rotation=errors[1],
                       ecart_reactions_SI=reaction_error)
        comparison.append(row)
    return comparison


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mesures-seules', action='store_true')
    p.add_argument('--journal-ci', type=Path)
    p.add_argument('--journal-tests', type=Path)
    p.add_argument('--journal-negatif', type=Path)
    args = p.parse_args()
    source = BANCS/'newton-redondant-mesures.json'
    comparison = bilan(json.loads(source.read_text()))
    (BANCS/'newton-redondant-comparaison.json').write_text(json.dumps(dict(
        entree_sha256=sha(source), comparaison=comparison), ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print('141 essais : 129 équilibres stricts, 12 échecs du témoin Princeton conservés.')
    if args.mesures_seules:
        return
    assert args.journal_ci and args.journal_tests and args.journal_negatif
    journal = args.journal_ci.read_text()
    for phrase in ('42 passed; 0 failed', '41/41 cas', '46/46 cas', '9/9 cas', '185 entrées', 'CI locale OK'):
        assert phrase in journal, phrase
    tests = args.journal_tests.read_text()
    assert 'Ran 64 tests' in tests and '\nOK\n' in tests
    negative = args.journal_negatif.read_text()
    assert 'FAILED (errors=1)' in negative and 'pas convergé en 100 itérations' in negative
    for src, suffix in [(args.journal_ci, 'ci.log'), (args.journal_tests, 'tests-python.log'),
                        (args.journal_negatif, 'princeton-negatif.log')]:
        dest = BANCS/f'newton-redondant-{suffix}'
        if src.resolve() != dest.resolve():
            shutil.copyfile(src, dest)
    current = sha(_vinkulum.__file__)
    princeton = json.loads((BANCS/'newton-redondant-princeton-strict.json').read_text())
    assert princeton['extension_sha256'] == current
    assert {(c['intervalles'], c['formulation']) for c in princeton['princeton']} == {
        (n, f) for n in (10, 20, 40, 60) for f in ('milieu', 'integree')}
    reports = [r for c in princeton['princeton'] for r in c['appels']]
    assert len(reports) == 400 and all(c['tolerance_tous_paliers'] for c in princeton['princeton'])
    for r in reports:
        assert r['erreur'] is None
        strict(r['rapport'])
    assert len(princeton['barres']) == 3
    assert all(r['erreur'] is None and r['rapport']['statut'] == 'tolerance' and r['rapport']['strict']
               for r in princeton['barres'])
    cable = json.loads((BANCS/'newton-redondant-cable.json').read_text())
    assert cable['extension_sha256'] == current
    assert all(r['statut'] == 'echec' and not r['tol_originale_tenue'] for r in cable['cas']['integree']['essais'].values())
    assert all(r['statut'] == 'ok' and r['tol_originale_tenue'] for r in cable['cas']['milieu']['essais'].values())
    files = list((ROOT/'src').rglob('*.rs')) + [ROOT/f for f in (
        'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'README.md', 'python/vinkulum/test_noyau.py',
        'ci/mesure_newton_redondant.py', 'ci/valide_newton_redondant.py', 'ci/bilan_contraintes.py',
        'ci/diagnostic_positions_compensees.py', 'ci/diagnostic_precision_statique.py',
        'ci/diagnostic_statut_statique.py', 'ci/diagnostic_cable.py', 'ci/diagnostic_poutres.py',
        'ci/local.sh', 'docs/NEWTON_REDONDANT.md', 'docs/OBJECTIF_MBDYN.md', *WORKERS)]
    artifacts = [f for f in BANCS.glob('newton-redondant-*') if f.name != 'newton-redondant-validation.json']
    manifest = dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        temoin_commit='d93593cc751d3aeaf523c0c11507d6a4caf4143f',
        binaires_sha256=dict(avant=BASELINE, courant=current),
        sources_sha256={str(f.relative_to(ROOT)): sha(f) for f in sorted(set(files))},
        artefacts_sha256={str(f.relative_to(ROOT)): sha(f) for f in sorted(artifacts)},
        ci=dict(commande='ci/local.sh --bancs', code_sortie=0, rust=42, verification=41,
                bancs=46, contacts=9, api=185,
                python_complementaire=dict(commande='python -c "from vinkulum.verification import _regressions; _regressions()"',
                    tests=64, code_sortie=0, raison='Test Princeton à encastrement double ajouté après la campagne générale.')),
        controles=dict(essais=141, equilibres_stricts=129, echecs_temoin_princeton=12, paliers_princeton_stricts=400),
        limites=['Réduction exacte limitée aux lignes identiques ou opposées, avec écarts compatibles.',
                 'Autres dépendances : replis SVD conservés ; parallélogrammes fermés mesurés.',
                 'Surcoût de détection sans lignes équivalentes publié.',
                 'Tolérance relative au déséquilibre initial libre, pas garantie absolue universelle.',
                 'Câble condensé extrême encore en échec ; aucune nouvelle exécution MBDyn ou Simpack.'])
    (BANCS/'newton-redondant-validation.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print('42 tests Rust, 64 Python, campagnes et 400 paliers Princeton validés ; manifeste écrit.')


if __name__ == '__main__':
    main()
