"""Lie les contrôles des contraintes aux binaires, sources et journaux mesurés."""
import argparse
import datetime
import json
import math
from pathlib import Path
import re
import shutil

import numpy as np
from vinkulum import _vinkulum
from bilan_contraintes import bilan, sha, strict

ROOT = Path(__file__).resolve().parents[1]
BANCS = ROOT/'docs/bancs'


def lit(nom):
    return json.loads((BANCS/f'contraintes-{nom}.json').read_text())


def directions(before, after):
    assert before['programme_sha256'] == after['programme_sha256'] == sha(
        ROOT/'ci/diagnostic_contraintes_proches.py')
    assert [x['delta'] for x in before['cas']] == [x['delta'] for x in after['cas']] == [
        1., 1e-3, 1e-6, 1e-9, 1e-12, 0.]
    out = []
    for a, b in zip(before['cas'], after['cas']):
        delta = b['delta']
        assert a['etat_initial'] == a['etat_final'] == b['etat_initial'] == b['etat_final']
        if delta == 0:
            assert all(x['erreur'] is not None and x['rapport']['statut'] == 'echec'
                       and x['rapport']['etat_restaure'] for x in (a, b))
            continue
        assert b['erreur'] is None
        strict(b['rapport'])
        angle = math.atan(delta)
        refs = [2.-1./math.tan(angle), 1./math.sin(angle)]
        reactions = dict(b['reactions'])
        got = [reactions['x'][0], reactions['proche'][0]]
        np.testing.assert_allclose(got, refs, rtol=2e-15, atol=1e-12)
        residual = max(abs(got[0]+got[1]*math.cos(angle)-2.), abs(got[1]*math.sin(angle)-1.))
        assert residual <= 1e-8
        if delta <= 1e-6:
            assert a['erreur'] is not None and a['rapport']['statut'] == 'echec'
        else:
            assert a['erreur'] is None
            strict(a['rapport'])
        out.append(dict(delta=delta, avant=a['rapport']['statut'], apres=b['rapport']['statut'],
                        residu_avant=a['rapport']['residu_libre'],
                        residu_apres=b['rapport']['residu_libre'],
                        residu_calcule_independamment_N=residual))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--journal-ci', required=True, type=Path)
    args = p.parse_args()
    journal = args.journal_ci.read_text()
    for phrase in ('35 passed; 0 failed', '41/41 cas', '46/46 cas', '9/9 cas',
                   '185 entrées', 'CI locale OK'):
        assert phrase in journal, phrase
    folders = re.findall(r'journaux : ([^\n]+)', journal)
    tests = [f for f in Path(folders[0]).glob('*.log') if 'Ran 60 tests' in f.read_text()]
    assert len(tests) == 1 and '\nOK\n' in tests[0].read_text()
    shutil.copyfile(args.journal_ci, BANCS/'contraintes-ci.log')
    shutil.copyfile(tests[0], BANCS/'contraintes-tests-python.log')
    negative = (BANCS/'contraintes-negatif-avant.log').read_text()
    assert 'Ran 1 test' in negative and 'FAILED (failures=1)' in negative

    current_sha = sha(_vinkulum.__file__)
    data = lit('alternees')
    assert data['binaires']['apres']['sha256'] == current_sha
    # Empreinte du témoin conservé au commit d'assemblage local.
    baseline_sha = '66f36f1af74288d1f5f198cc02025663c8f90d977311e355661f0926aa8815ca'
    assert data['binaires']['avant']['sha256'] == baseline_sha
    comparison = bilan(data)
    recorded = lit('comparaison')
    assert recorded['comparaison'] == comparison
    assert recorded['entree_sha256'] == sha(BANCS/'contraintes-alternees.json')
    assert recorded['programme_sha256'] == sha(ROOT/'ci/bilan_contraintes.py')
    before, after = lit('proches-avant'), lit('proches-apres')
    assert before['extension_sha256'] == baseline_sha and after['extension_sha256'] == current_sha
    close = directions(before, after)

    princeton = lit('princeton-strict')
    assert princeton['extension_sha256'] == current_sha
    assert {(s['intervalles'], s['formulation']) for s in princeton['princeton']} == {
        (n, f) for n in (10, 20, 40, 60) for f in ('milieu', 'integree')}
    reports = [r for s in princeton['princeton'] for r in s['appels']]
    assert len(reports) == 400
    assert all(s['tolerance_tous_paliers'] for s in princeton['princeton'])
    for sample in reports:
        assert sample['erreur'] is None
        strict(sample['rapport'])
    assert len(princeton['barres']) == 3 and all(s['erreur'] is None and
        s['rapport']['statut'] == 'tolerance' for s in princeton['barres'])
    cable = lit('cable')
    assert cable['extension_sha256'] == current_sha
    assert all(s['statut'] == 'echec' and not s['tol_originale_tenue']
               for s in cable['cas']['integree']['essais'].values())
    assert all(s['statut'] == 'ok' and s['tol_originale_tenue']
               for s in cable['cas']['milieu']['essais'].values())
    timings = lit('princeton-alterne')
    assert {v: x['sha256'] for v, x in timings['binaires'].items()} == {
        'avant': baseline_sha, 'apres': current_sha}
    assert timings['programme_sha256'] == sha(ROOT/'ci/mesure_statique_alternee.py')
    assert timings['worker_sha256'] == sha(ROOT/'ci/diagnostic_assemblage_local.py')
    assert len(timings['mesures']) == 12
    for n in (240, 480):
        group = [r for r in timings['mesures'] if r['intervalles'] == n]
        assert {(r['version'], r['repetition']) for r in group} == {
            (v, rep) for v in ('avant', 'apres') for rep in range(3)}
        for r in group:
            assert r['erreur'] is None
            for s in r['statique_rapports']:
                strict(s)
        tips = [np.asarray(r['position']) for r in group]
        assert max(float(np.max(np.abs(a-b))) for a in tips for b in tips) <= 1e-10

    files = list((ROOT/'src').rglob('*.rs'))+list((ROOT/'examples').rglob('*.rs'))
    files += [ROOT/f for f in ('Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'README.md',
        'docs/CONTRAINTES_CREUSES.md', 'docs/OBJECTIF_MBDYN.md', 'python/vinkulum/test_noyau.py',
        'ci/diagnostic_contraintes.py', 'ci/diagnostic_contraintes_proches.py',
        'ci/diagnostic_chaine_articulee.py', 'ci/mesure_contraintes_alternee.py',
        'ci/bilan_contraintes.py', 'ci/valide_contraintes.py', 'ci/diagnostic_positions_compensees.py',
        'ci/diagnostic_precision_statique.py', 'ci/diagnostic_statut_statique.py',
        'ci/diagnostic_cable.py', 'ci/diagnostic_poutres.py', 'ci/confronte_mbdyn.py',
        'ci/mesure_statique_alternee.py', 'ci/diagnostic_assemblage_local.py', 'ci/local.sh')]
    artifacts = [f for f in BANCS.glob('contraintes-*') if f.name != 'contraintes-validation.json']
    manifest = dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        normalisation_archives='JSON alterné compacté sans changement des valeurs ; '
                                'espaces terminaux retirés du journal négatif.',
        temoin_commit='4b79a761c47d3de2d74018808258b38d9c5958db',
        binaires_sha256=dict(avant=baseline_sha, courant=current_sha),
        sources_sha256={str(f.relative_to(ROOT)): sha(f) for f in sorted(files)},
        artefacts_sha256={str(f.relative_to(ROOT)): sha(f) for f in sorted(artifacts)},
        ci=dict(commande='ci/local.sh --bancs', code_sortie=0, rust=35, python=60,
                verification=41, bancs_rapides=46, contacts=9, api=185),
        controles=dict(equilibres_alternes=132, topologies=4, proches=close,
                       directions_confondues_refusees=True, paliers_princeton_stricts=400,
                       barres_allongement_exact=3, controle_negatif_ancien_binaire=1,
                       cable_condense_original_et_lineaire='echec'),
        limites=['Comparaisons entre versions Vinkulum ; aucun nouveau classement externe.',
                 'Repli dense par composante pour les projections redondantes.',
                 'Matrices publiques et autres projecteurs restent denses.',
                 'Critères statiques et échelles physiques inchangés ; pas de garantie universelle de rang.',
                 'Câble condensé extrême en échec ; couverture Simpack non acquise.'])
    (BANCS/'contraintes-validation.json').write_text(json.dumps(manifest, ensure_ascii=False,
                                                             indent=2, allow_nan=False)+'\n')
    print('Validation : 35 Rust, 60 Python, 132 équilibres comparés, 400 paliers stricts.')


if __name__ == '__main__':
    main()
