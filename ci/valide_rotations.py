"""Vérifie les campagnes de rotations et lie leurs preuves au binaire installé."""
import datetime
import gzip
import json
from pathlib import Path
import subprocess
import tempfile

from vinkulum import _vinkulum
from bilan_contraintes import sha, strict
from bilan_factorisations import bilan as statique
from bilan_rotations import bilan as dynamique
from valide_echelle_svd import valide as echelles

ROOT = Path(__file__).resolve().parents[1]
BANCS = ROOT/'docs/bancs'
AVANT = '74c52b7bd4d8aaa9743d84c2845296634bf11e81fdde26c3f56686c85a60e882'
POLAIRE = 'bbf620524ef647ccae3d9450d7a1ecdf65676d1d68efcbce42d46e13ee5c8acb'
INERTIE = '71d323ee30e9c282d8dbb508c9a52664164f4be01806b1c1230dfe7937621492'
CAS = [('cascade', 1), ('cascade', 4), ('parallelogrammes', 1), ('parallelogrammes', 8),
       ('double', 64), ('mixte', 64), ('articulee', 128), ('spatiale', 128), ('chaine', 128),
       ('rotors', 128), ('boucle', 512), ('princeton_milieu', 10), ('princeton_integree', 10),
       ('cascade', 8), ('cascade', 16), ('cascade', 32), ('parallelogrammes', 16), ('parallelogrammes', 32)]


def lit(name):
    path = BANCS/f'rotations-{name}'
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes())


def sauve(name, data):
    (BANCS/f'rotations-{name}.json').write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def main():
    current = sha(_vinkulum.__file__)
    reports = {}
    specs = {
        'statique-premiere': (statique, dict(avant=AVANT, apres=POLAIRE)),
        'statique': (statique, dict(avant=AVANT, apres=current)),
        'dynamique-premiere': (dynamique, dict(avant=AVANT, apres=POLAIRE)),
        'dynamique-carres': (dynamique, dict(avant=AVANT, polaire=POLAIRE, apres=current)),
        'dynamique-inertie': (dynamique, dict(avant=AVANT, polaire=POLAIRE, carres=current, apres=INERTIE)),
    }
    for name, (validate, binaries) in specs.items():
        data = lit(name+'.json.gz')
        assert {k: v['sha256'] for k, v in data['binaires'].items()} == binaries
        if validate is statique:
            assert {tuple(c) for c in data['cas']} == set(CAS) and len(data['cas']) == len(CAS)
            assert data['delai_processus_s'] == 90.
        report = validate(data)
        if validate is statique:
            assert report['essais'] == report['equilibres_stricts'] == 108
            assert not report['echecs_restaures'] and not report['delais_depasses']
        reports[name] = report
        report.update(entree_sha256=sha(BANCS/f'rotations-{name}.json.gz'),
                      programme_sha256=sha(ROOT/f'ci/{"bilan_factorisations" if validate is statique else "bilan_rotations"}.py'))
        sauve(name+'-bilan', report)

    numerical = echelles(lit('echelles-rust.json'))
    numerical.update(entree_sha256=sha(BANCS/'rotations-echelles-rust.json'),
                     programme_sha256=sha(ROOT/'ci/valide_echelle_svd.py'))
    sauve('echelles-lapack', numerical)
    executable = lit('echelles-executable.json')
    assert executable['sources_sha256'] == {
        str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT/'src').rglob('*.rs'))}
    assert executable['exemple_sha256'] == sha(ROOT/'examples/diagnostic_echelle_svd.rs')
    assert executable['executable_sha256'] == sha(ROOT/'target/release/examples/diagnostic_echelle_svd')
    for stage, python_count in [('premiere', 65), ('finale', 66)]:
        log = (BANCS/f'rotations-ci-{stage}.log').read_text()
        for phrase in ('52 passed; 0 failed', '41/41 cas', '46/46 cas', '9/9 cas', '185 entrées', 'CI locale OK'):
            assert phrase in log, (stage, phrase)
        tests = (BANCS/f'rotations-tests-python-{stage}.log').read_text()
        assert f'Ran {python_count} tests' in tests and '\nOK\n' in tests

    princeton = lit('princeton.json')
    assert princeton['extension_sha256'] == current
    assert {(c['intervalles'], c['formulation']) for c in princeton['princeton']} == {
        (n, f) for n in (10, 20, 40, 60) for f in ('milieu', 'integree')}
    calls = [r for c in princeton['princeton'] for r in c['appels']]
    assert len(calls) == 400 and all(c['tolerance_tous_paliers'] for c in princeton['princeton'])
    for call in calls:
        assert call['erreur'] is None
        strict(call['rapport'])
    assert len(princeton['barres']) == 3
    assert all(r['erreur'] is None and r['rapport']['statut'] == 'tolerance' and r['rapport']['strict']
               for r in princeton['barres'])
    cable = lit('cable.json')
    assert cable['extension_sha256'] == current
    assert all(r['statut'] == 'echec' and not r['tol_originale_tenue'] for r in cable['cas']['integree']['essais'].values())
    assert all(r['statut'] == 'ok' and r['tol_originale_tenue'] for r in cable['cas']['milieu']['essais'].values())

    prototypes = lit('prototypes.json')
    for name, binary in [('premiere', POLAIRE), ('inertie', INERTIE)]:
        proto = prototypes[name]
        assert proto['extension_sha256'] == binary
        with tempfile.TemporaryDirectory(prefix='vinkulum-rotation-source-') as temp:
            dest = Path(temp)
            for file in proto['sources_rust_sha256']:
                target = dest/file
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes((ROOT/file).read_bytes())
            subprocess.run(['git', 'apply', '--unidiff-zero', str(ROOT/proto['patch'])], cwd=dest, check=True)
            assert {f: sha(dest/f) for f in proto['sources_rust_sha256']} == proto['sources_rust_sha256']

    files = list((ROOT/'src').rglob('*.rs')) + list((ROOT/'ci').glob('*.py')) + [ROOT/f for f in (
        'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'README.md', 'python/vinkulum/test_noyau.py',
        'examples/diagnostic_echelle_svd.rs', 'ci/local.sh', 'docs/ROTATIONS_POLAIRES.md', 'docs/OBJECTIF_MBDYN.md')]
    artifacts = [p for p in BANCS.glob('rotations-*') if p.name != 'rotations-validation.json']
    manifest = dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        temoin_commit='f445fc363ee59512531629297f093724220a8ab8',
        binaires_sha256=dict(avant=AVANT, premiere=POLAIRE, courant=current, prototype_inertie_retire=INERTIE),
        sources_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(set(files))},
        artefacts_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(artifacts)},
        ci=dict(commande='ci/local.sh --bancs', code_sortie=0, rust=52, python=66, verification=41, bancs=46, contacts=9, api=185),
        essais_statiques_stricts=216, simulations_dynamiques=432, paliers_princeton_stricts=400,
        mesures='Processus neufs alternés, même CPU, fils demandés à un ; trois répétitions. Les temps de CI et Princeton ne servent pas à comparer les vitesses.',
        limites=['Projection polaire conditionnée au voisinage de SO(3), sans certificat par intervalles.',
                 'SVD : une perte lors de la normalisation ou une valeur singulière non représentable provoque un refus.',
                 'Écarts entre repères à long terme encore présents ; référence temporelle indépendante limitée à une seconde.',
                 'Des surcoûts demeurent. La préparation des inverses ralentit le témoin testé et a été retirée.',
                 'Câble condensé extrême encore en échec ; aucune nouvelle exécution MBDyn ou Simpack.'])
    sauve('validation', manifest)
    print('Rotations : 216 équilibres stricts, 432 simulations, cinq échelles SVD et empreintes validés.')


if __name__ == '__main__':
    main()
