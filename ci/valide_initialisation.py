"""Vérifie la provenance des mesures et de la livraison 0.7.2."""
import argparse
import base64
import csv
import datetime
from email.parser import BytesParser
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tomllib
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BANCS = ROOT/'docs/bancs'
BASE = 'afcc8c573ef3cff53258f9499ba4084ed8d74888'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--roue', type=Path, required=True)
    p.add_argument('--controle-roue', type=Path, required=True)
    args = p.parse_args()
    from vinkulum import _vinkulum, __version__
    assert __version__ == '0.7.2'
    version = '0.7.2'
    cfg = tomllib.loads((ROOT/'Cargo.toml').read_text())
    py = tomllib.loads((ROOT/'pyproject.toml').read_text())
    lock = tomllib.loads((ROOT/'Cargo.lock').read_text())
    assert cfg['package']['version'] == py['project']['version'] == version
    assert next(p['version'] for p in lock['package'] if p['name'] == 'vinkulum') == version
    ancien = tomllib.loads(subprocess.check_output(['git', 'show', f'{BASE}:Cargo.lock'], cwd=ROOT, text=True))
    for pack in lock['package']:
        if pack['name'] == 'vinkulum':
            pack['version'] = '0.7.1'
    assert lock == ancien, 'dépendances modifiées'
    current = sha(_vinkulum.__file__)
    witness = json.loads((BANCS/'version-0.7.1.json').read_text())['extension_ci_sha256']
    data = json.loads((BANCS/'initialisation-mesures.json').read_text())
    assert len(data['resultats']) == 138 and len(data['bilan']) == 23
    assert all(l['apres']['valide'] for l in data['bilan'])
    for path, digest in data['sources_sha256'].items():
        assert sha(ROOT/path) == digest, path
    for path, digest in data['versions_sources_sha256'].items():
        assert sha(ROOT/path) == digest, path
    assert data['programme_sha256'] == sha(ROOT/'ci/mesure_initialisation.py')
    assert data['diagnostic_sha256'] == sha(ROOT/'ci/diagnostic_initialisation.py')
    for r in data['resultats']:
        avant = r['version_comparaison'] == 'avant'
        assert r['version'] == ('0.7.1' if avant else version)
        assert r['binaire_sha256'] == (witness if avant else current)
        if not avant:
            assert r['erreur'] is None and r['secondes'] > 0 and r['rss_kib'] > 0
            c = r['controles']
            if r['cas'] != 'libre':
                assert c['etat_inchange']
            if r['cas'] == 'analyse':
                assert c['mobilites'] == r['cellules']
            for k, value in c.items():
                if k.startswith('erreur_'):
                    assert 0 <= value <= 5e-10, (r['cas'], k, value)
    refs = json.loads((BANCS/'initialisation-0.7.2-references.json').read_text())
    assert {r['cas'] for r in refs} == {'courbure', 'commande', 'decollement'}
    for r in refs:
        assert r['version'] == version and r['binaire_sha256'] == current
        if r['cas'] == 'courbure':
            assert all(v['erreur_relative'] < 2e-14 for v in r['resultats'])
        elif r['cas'] == 'commande':
            assert r['erreur'] < 1e-13
        else:
            assert all(v['erreur'] < 1e-12 for v in r['resultats'])
    ci = (BANCS/'initialisation-ci.log').read_text()
    for token in ('66 passed; 0 failed', '41/41 cas', '46/46 cas', '9/9 cas', '185 entrées', 'CI locale OK'):
        assert token in ci, token
    tests = (BANCS/'initialisation-tests-python.log').read_text()
    assert 'Ran 72 tests' in tests and '\nOK\n' in tests
    fresh = json.loads(args.controle_roue.read_text())
    assert fresh['version'] == version and fresh['extension_sha256'] == current
    assert fresh['import_hors_depot'] and fresh['exemple_readme_execute']
    assert fresh['readme_sha256'] == sha(ROOT/'README.md')
    fresh_log = (BANCS/'initialisation-roue-verification.log').read_text()
    assert '41/41 cas' in fresh_log
    fresh_tests = (BANCS/'initialisation-roue-tests-python.log').read_text()
    assert 'Ran 72 tests' in fresh_tests and '\nOK\n' in fresh_tests
    with zipfile.ZipFile(args.roue) as z:
        metadata = BytesParser().parsebytes(z.read(f'vinkulum-{version}.dist-info/METADATA'))
        assert metadata['Version'] == version and metadata['Requires-Python'] == '>=3.14'
        assert metadata.get_payload(decode=True).decode('utf8') == (ROOT/'README.md').read_text()+'\n'
        record = z.read(f'vinkulum-{version}.dist-info/RECORD').decode()
        count = 0
        for path, digest, size in csv.reader(io.StringIO(record)):
            if not digest:
                continue
            payload = z.read(path)
            expected = 'sha256='+base64.urlsafe_b64encode(hashlib.sha256(payload).digest()).decode().rstrip('=')
            assert digest == expected and int(size) == len(payload), path
            count += 1
        extension = next(f for f in z.namelist() if f.endswith('.so') and '_vinkulum' in f)
        assert hashlib.sha256(z.read(extension)).hexdigest() == current
    sources = list((ROOT/'src').rglob('*.rs')) + [ROOT/f for f in (
        'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'README.md', 'docs/API.md',
        'docs/VERSION_0.7.2.md', 'docs/INITIALISATION.md', 'docs/OBJECTIF_MBDYN.md',
        'ci/local.sh', 'ci/diagnostic_initialisation.py', 'ci/mesure_initialisation.py',
        'ci/valide_initialisation.py', 'ci/controle_roue_initialisation.py', 'python/vinkulum/test_noyau.py')]
    artifacts = list(BANCS.glob('initialisation-*'))
    manifest = dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), version=version,
        tag_annote=f'v{version}', base_commit=BASE, dependances_inchangees=True,
        extension_ci_sha256=current, roue=dict(nom=args.roue.name, sha256=sha(args.roue),
            octets=args.roue.stat().st_size, entrees_record_verifiees=count),
        readme_roue_conforme=True, controle_roue=fresh,
        verification_roue=dict(verification=41, python=72),
        ci=dict(commande='ci/local.sh --bancs', code_sortie=0, rust=66, python=72,
                verification=41, bancs=46, contacts=9, api=185),
        comparaison=dict(essais=138, configurations=23, candidats_valides=69,
                         echecs_temoins=sum(r['erreur'] is not None for r in data['resultats'])),
        sources_sha256={str(f.relative_to(ROOT)):sha(f) for f in sorted(sources)},
        artefacts_sha256={str(f.relative_to(ROOT)):sha(f) for f in sorted(artifacts)},
        portee='Initialisation et analyse. Aucune nouvelle exécution MBDyn ou Simpack.')
    (BANCS/f'version-{version}.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print('Sources, mesures, références, CI et roue 0.7.2 vérifiées.')


if __name__ == '__main__':
    main()
