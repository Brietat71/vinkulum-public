"""Lie mesures, régressions, prototype conservé et CI aux sources validées."""
import datetime
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import tempfile

from vinkulum import _vinkulum, __version__
from bilan_contraintes import etat, sha, strict
from bilan_raffinement import bilan

ROOT = Path(__file__).resolve().parents[1]
BANCS = ROOT/'docs/bancs'
BASE = '5d0559ca3766792e361d9c28a3131c59db636491'


def lit(name):
    return json.loads((BANCS/f'raffinement-{name}.json').read_text())


def sauve(name, data):
    (BANCS/f'raffinement-{name}.json').write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def main():
    current = sha(_vinkulum.__file__)
    assert __version__ == '0.7.1'
    baseline = json.loads((BANCS/'version-0.7.0.json').read_text())['extension_ci_sha256']
    for name in ('prototype-mesures', 'mesures'):
        data = lit(name)
        assert data['versions']['avant']['extension_sha256'] == baseline
        if name == 'mesures':
            assert data['versions']['apres']['extension_sha256'] == current
            report = bilan(data)
        else:
            with tempfile.TemporaryDirectory(prefix='vinkulum-raffinement-sources-') as temp:
                dest = Path(temp)
                for file in data['sources_candidat_sha256']:
                    if file in ('src/equilibrage.rs', 'src/raffinement.rs'):
                        continue
                    target = dest/file
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(subprocess.check_output(['git', 'show', f'{BASE}:{file}'], cwd=ROOT))
                subprocess.run(['git', 'apply', str(BANCS/'raffinement-prototype.patch')], cwd=dest, check=True)
                report = bilan(data, source_root=dest)
        assert report['essais'] == report['equilibres_stricts'] == 156
        assert not report['echecs_conserves']
        report.update(mesures_sha256=sha(BANCS/f'raffinement-{name}.json'),
                      validateur_sha256=sha(ROOT/'ci/bilan_raffinement.py'))
        sauve('prototype-bilan' if name.startswith('prototype') else 'bilan', report)

    princeton = lit('princeton')
    assert princeton['extension_sha256'] == current
    for path, digest in princeton['worktree_lanceur_sha256'].items():
        assert sha(ROOT/path) == digest, path
    cases = princeton['princeton']
    assert len(cases) == 8 and {(c['intervalles'], c['formulation']) for c in cases} == {
        (n, f) for n in (10, 20, 40, 60) for f in ('milieu', 'integree')}
    for case in cases:
        assert case['tolerance_tous_paliers'] and case['charge_finale_acceptee'] == 1.
        assert len(case['appels']) == 50
        etat(dict(etat=case['etat_final']))
        for index, call in enumerate(case['appels'], 1):
            assert call['palier_externe'] == index and call['erreur'] is None
            strict(call['rapport'])
    assert {r['element'] for r in princeton['barres']} == {'milieu', 'integree', 'super'}
    assert len(princeton['barres']) == 3
    for r in princeton['barres']:
        assert r['erreur'] is None and r['rapport']['statut'] == 'tolerance'
        assert r['rapport']['strict'] and r['rapport']['tol'] == 1e-14
        assert r['rapport']['residu_relatif'] <= 1e-14 and r['rapport']['contraintes'] <= 1e-14
        etat(r)
        state = r['etat']
        elongation = sum(Fraction(state[k][1][0])-Fraction(state[k][0][0]) for k in (1, 6))-1
        assert elongation == Fraction(1, 2**60)

    diagnostic = lit('diagnostic')
    assert diagnostic['programme_sha256'] == sha(ROOT/'ci/diagnostic_raffinement.py')
    for path, digest in diagnostic['archives_sha256'].items():
        assert sha(ROOT/path) == digest, path
    assert diagnostic['dimension'] == 1216 and diagnostic['physiques'] == 576
    assert diagnostic['brut']['rang_au_seuil'] == diagnostic['blocs']['rang_au_seuil'] == 1120
    assert diagnostic['jacobi']['rang_au_seuil'] > 1120
    assert diagnostic['brut']['asymetrie_max'] > .05
    assert diagnostic['lsmr_jacobi']['code_arret'] == diagnostic['lsmr_brut']['code_arret'] == 7
    trials = diagnostic['raffinement_blocs']
    assert not trials[0]['tolerance_atteinte'] and trials[0]['residus'][-1] > trials[0]['residus'][0]
    assert all(r['tolerance_atteinte'] for r in trials[1:])

    ci = (BANCS/'raffinement-ci.log').read_text()
    for phrase in ('57 passed; 0 failed', '41/41 cas', '46/46 cas', '9/9 cas', '185 entrées', 'CI locale OK'):
        assert phrase in ci, phrase
    tests = (BANCS/'raffinement-tests-python.log').read_text()
    assert 'Ran 67 tests' in tests and '\nOK\n' in tests
    files = list((ROOT/'src').rglob('*.rs')) + list((ROOT/'ci').glob('*raffinement.py')) + [ROOT/f for f in (
        'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'README.md', 'python/vinkulum/test_noyau.py',
        'ci/diagnostic_cascade.py', 'ci/local.sh', 'docs/RAFFINEMENT_STATIQUE.md', 'docs/OBJECTIF_MBDYN.md')]
    artifacts = [p for p in BANCS.glob('raffinement-*') if p.name != 'raffinement-validation.json']
    sauve('validation', dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        base_commit=BASE, version=__version__, extension_sha256=current,
        sources_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(set(files))},
        artefacts_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(artifacts)},
        ci=dict(commande='ci/local.sh --bancs', code_sortie=0, rust=57, python=67,
                verification=41, bancs=46, contacts=9, api=185),
        equilibres_comparaison=156, equilibres_prototype=156, paliers_princeton=400, barres_precision=3,
        portee='Comparaison interne avec la roue 0.7.0. Aucun nouveau calcul MBDyn ou Simpack.'))
    print('Validés : 312 équilibres alternés, 400 paliers stricts, 3 barres précises, CI et empreintes.')


if __name__ == '__main__':
    main()
