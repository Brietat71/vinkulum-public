"""Rassemble les preuves numériques, mécaniques et les empreintes de cette étape."""
import datetime
import gzip
import hashlib
import json
from pathlib import Path

from vinkulum import _vinkulum
from bilan_contraintes import sha, strict
from bilan_factorisations import bilan
from valide_svd import valide

ROOT = Path(__file__).resolve().parents[1]
BANCS = ROOT/'docs/bancs'
AVANT = 'f7efce93e24837af8f7b3ae721533c4fd7be704522ccc68f6257bd7c157cd96a'
SVD = 'ac8b4c9090108668b73ed29d09a30aeb16fab7ff3a6bcd2052435190814870cf'
COD = '8077da2c497d560ed3d64924c48dd8d0079b0c99049feb8ec42b48d395945735'
BLOCS = '66b30802b8e82bff3287e7ce25107aa264cd9af4db47de76c9b3cff51c9d2fc3'
COMMUNS = [('cascade', 1), ('cascade', 4), ('parallelogrammes', 1), ('parallelogrammes', 8),
           ('double', 64), ('mixte', 64), ('articulee', 128), ('spatiale', 128), ('chaine', 128),
           ('rotors', 128), ('boucle', 512), ('princeton_milieu', 10), ('princeton_integree', 10)]
GRANDS = [('cascade', 8), ('cascade', 16), ('cascade', 32), ('parallelogrammes', 16), ('parallelogrammes', 32)]
ABLATION = [('cascade', n) for n in (4, 8, 16, 32)] + [('parallelogrammes', n) for n in (8, 16, 32)] + [
            ('double', 64), ('spatiale', 128), ('princeton_milieu', 10)]


def lit(suffix):
    return json.loads((BANCS/f'factorisations-{suffix}.json').read_text())


def main():
    courant = sha(_vinkulum.__file__)
    specs = {
        'ablation': (dict(svd=SVD, cod=COD), ABLATION),
        'blocs': (dict(avant=AVANT, apres=BLOCS), COMMUNS+[(f, n) for f, n in GRANDS if (f, n) != ('cascade', 32)]),
        'mesures': (dict(avant=AVANT, apres=courant), COMMUNS),
        'grands': (dict(apres=courant), GRANDS),
    }
    reports = {}
    for name, (binaries, cases) in specs.items():
        data = lit(name)
        assert {k: v['sha256'] for k, v in data['binaires'].items()} == binaries
        assert {tuple(c) for c in data['cas']} == set(cases)
        assert data['delai_processus_s'] == 90
        report = bilan(data)
        report.update(entree_sha256=sha(BANCS/f'factorisations-{name}.json'),
                      programme_sha256=sha(ROOT/'ci/bilan_factorisations.py'))
        (BANCS/f'factorisations-{name}-bilan.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n')
        reports[name] = {k: v for k, v in report.items() if k not in ('comparaison', 'programme_sha256', 'entree_sha256')}
    assert reports['ablation']['equilibres_stricts'] == 48
    assert len(reports['ablation']['echecs_restaures']) == 12 and not reports['ablation']['delais_depasses']
    for name in ('mesures', 'grands'):
        assert not reports[name]['echecs_restaures'] and not reports[name]['delais_depasses']
    assert len(reports['blocs']['echecs_restaures']) == 9 and len(reports['blocs']['delais_depasses']) == 3

    startup = lit('demarrage')
    for name, binary in [('blocs', BLOCS), ('finale', courant)]:
        item = startup[name]
        trace = gzip.decompress((BANCS/f'factorisations-{name}-openat.log.gz').read_bytes())
        assert hashlib.sha256(trace).hexdigest() == item['trace_sha256']
        assert item['extension_sha256'] == binary
        lines = trace.decode().splitlines()
        assert item['ouvertures_cache'] == sum('/sys/devices/system/cpu/' in l and '/cache' in l for l in lines)
        strict(item['rapport'])
    assert startup['blocs']['ouvertures_cache'] > 0 and startup['finale']['ouvertures_cache'] == 0

    numerical = valide(lit('svd-rust'))
    archived = lit('svd-lapack')
    assert numerical == {k: v for k, v in archived.items() if k not in ('entree_sha256', 'programme_sha256')}
    assert archived['programme_sha256'] == sha(ROOT/'ci/valide_svd.py')
    assert archived['entree_sha256'] == sha(BANCS/'factorisations-svd-rust.json')

    ci = (BANCS/'factorisations-ci.log').read_text()
    for phrase in ('48 passed; 0 failed', '41/41 cas', '46/46 cas', '9/9 cas', '185 entrées', 'CI locale OK'):
        assert phrase in ci, phrase
    tests = (BANCS/'factorisations-tests-python.log').read_text()
    assert 'Ran 65 tests' in tests and '\nOK\n' in tests
    negative = (BANCS/'factorisations-negatif.log').read_text()
    assert 'FAILED (errors=1)' in negative and 'pas convergé en 30 itérations' in negative

    princeton = lit('princeton')
    assert princeton['extension_sha256'] == courant
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
    cable = lit('cable')
    assert cable['extension_sha256'] == courant
    assert all(r['statut'] == 'echec' and not r['tol_originale_tenue'] for r in cable['cas']['integree']['essais'].values())
    assert all(r['statut'] == 'ok' and r['tol_originale_tenue'] for r in cable['cas']['milieu']['essais'].values())

    files = list((ROOT/'src').rglob('*.rs')) + list((ROOT/'ci').glob('*.py')) + [ROOT/f for f in (
        'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'README.md', 'python/vinkulum/test_noyau.py',
        'examples/diagnostic_svd.rs', 'ci/local.sh', 'docs/FACTORISATIONS_ORTHOGONALES.md', 'docs/OBJECTIF_MBDYN.md')]
    artifacts = [p for p in BANCS.glob('factorisations-*') if p.name != 'factorisations-validation.json']
    manifest = dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        temoin_commit='6d38aa435c38360df23d0e48315300d3e6cb48ba',
        binaires_sha256=dict(avant=AVANT, svd_seule=SVD, prototype_cod=COD, prototype_blocs=BLOCS, courant=courant),
        sources_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(set(files))},
        artefacts_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(artifacts)},
        ci=dict(commande='ci/local.sh --bancs', code_sortie=0, rust=48, python=65, verification=41, bancs=46, contacts=9, api=185),
        campagnes=reports, paliers_princeton_stricts=400, controles_svd=numerical['controles'],
        limites=['Le repli de rang général reste dense.',
                 'La borne de rang concerne le triangle calculé, sans certificat par intervalles des données.',
                 'Les échecs et délais dépassés des témoins ne sont pas des gains de vitesse.',
                 'Tolérance relative au déséquilibre initial libre ; pas de garantie absolue universelle.',
                 'Câble condensé extrême encore en échec ; aucune nouvelle exécution MBDyn ou Simpack.'])
    (BANCS/'factorisations-validation.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print('Factorisations : preuves algébriques, mécaniques, campagnes et empreintes validées.')


if __name__ == '__main__':
    main()
