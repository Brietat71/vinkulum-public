"""Vérifie les archives d'assemblage et lie les résultats aux sources et journaux."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil

import numpy as np
from vinkulum import _vinkulum
from bilan_assemblage_local import bilan
from mesure_dynamique_locale import modele

ROOT = Path(__file__).resolve().parents[1]
BANCS = ROOT / 'docs/bancs'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lit(nom):
    return json.loads((BANCS / f'assemblage-local-{nom}.json').read_text())


def ecrit(nom, valeur):
    (BANCS / f'assemblage-local-{nom}.json').write_text(
        json.dumps(valeur, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def dynamique(avant, apres):
    for key in ('cpu', 'pas_s', 'duree_s', 'programme_sha256'):
        assert avant[key] == apres[key], key
    assert avant['cas'].keys() == apres['cas'].keys()
    assert apres['programme_sha256'] == sha(ROOT / 'ci/mesure_dynamique_locale.py')
    out = {}
    labels = ('temps_s', 'position_haute_m', 'orientation', 'vitesse_m_s',
              'vitesse_angulaire_rad_s', 'inflow_m_s', 'position_basse_m')
    for nom, a in avant['cas'].items():
        b = apres['cas'][nom]
        assert len(a['echantillons']) == len(b['echantillons']) == 5
        errors, exact = {}, True
        for index, label in enumerate(labels):
            values = []
            for s in a['echantillons']:
                for t in b['echantillons']:
                    x, y = (np.asarray(e['etat'][index], dtype=np.float64) for e in (s, t))
                    assert x.shape == y.shape and np.isfinite(x).all() and np.isfinite(y).all()
                    values.append(float(np.max(np.abs(x-y), initial=0)))
                    exact = exact and x.tobytes() == y.tobytes()
            errors[label] = max(values)
        assert max(errors.values()) < 2e-14, (nom, errors)
        out[nom] = dict(avant_s=a['mediane_s'], apres_s=b['mediane_s'],
                        variation_temps_pourcent=100*(b['mediane_s']/a['mediane_s']-1),
                        stats_avant=[s['stats'] for s in a['echantillons']],
                        stats_apres=[s['stats'] for s in b['echantillons']],
                        ecarts_max_par_composante=errors,
                        etat_exporte_identique_bit_a_bit=exact)
    assert out['uniforme']['etat_exporte_identique_bit_a_bit']
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--journal-ci', type=Path, required=True)
    args = parser.parse_args()
    journal = args.journal_ci.read_text()
    for message in ('31 passed; 0 failed', '41/41 cas', '46/46 cas',
                    '9/9 cas', '185 entrées', 'CI locale OK'):
        assert message in journal, message
    dossiers = re.findall(r'journaux : ([^\n]+)', journal)
    assert dossiers
    tests = [f for f in Path(dossiers[0]).glob('*.log') if 'Ran 59 tests' in f.read_text()]
    assert len(tests) == 1 and '\nOK\n' in tests[0].read_text()
    shutil.copyfile(args.journal_ci, BANCS / 'assemblage-local-ci.log')
    shutil.copyfile(tests[0], BANCS / 'assemblage-local-tests-python.log')
    assert 'FAILED (failures=2)' in (BANCS / 'assemblage-local-regressions-avant.log').read_text()

    before, after, alternate = map(lit, ('avant', 'apres', 'alterne'))
    current_sha = sha(_vinkulum.__file__)
    assert after['extension_sha256'] == current_sha
    assert alternate['binaires']['apres']['sha256'] == current_sha
    assert before['extension_sha256'] == alternate['binaires']['avant']['sha256']
    assert after['programme_sha256'] == sha(ROOT / 'ci/diagnostic_assemblage_local.py')
    assert after['modeles_et_controles_sha256'] == sha(ROOT / 'python/vinkulum/test_noyau.py')
    assert alternate['programme_sha256'] == sha(ROOT / 'ci/mesure_statique_alternee.py')
    assert alternate['worker_sha256'] == after['programme_sha256']
    calculated = bilan(before, after)
    recorded = lit('comparaison')
    for key, value in calculated.items():
        assert recorded[key] == value, key
    for path, digest in recorded['entrees_sha256'].items():
        assert sha(ROOT / path) == digest, path
    groups = []
    for nom, source in [('avant', before), ('apres', after)]:
        groups.append(dict(source, statique=[dict(c, code_sortie=0)
                      for c in alternate['mesures'] if c['version'] == nom]))
    alt_comparison = bilan(*groups)['statique']
    assert all(c['ecart_positions_m'] == 0 for c in alt_comparison)
    assert all(c['iterations_last_attempt'] == 36
               and c['statique_rapports'][0]['evaluations'] == 48
               and c['statique_rapports'][0]['tentatives'] == 2 for c in alternate['mesures'])

    for nom in ('princeton-strict', 'cable', 'dynamique-apres', 'tangentes-princeton'):
        assert lit(nom)['extension_sha256'] == current_sha, nom
    assert lit('statique-princeton')['metadata']['extension_sha256'] == current_sha
    da, db = map(lit, ('dynamique-avant', 'dynamique-apres'))
    assert da['extension_sha256'] == before['extension_sha256']
    dyn = dynamique(da, db)
    p = lit('princeton-strict')
    assert {(c['intervalles'], c['formulation']) for c in p['princeton']} == {
        (n, f) for n in (10, 20, 40, 60) for f in ('milieu', 'integree')}
    reports = [x for c in p['princeton'] for x in c['appels']]
    assert len(reports) == 400 and all(c['tolerance_tous_paliers'] for c in p['princeton'])
    assert all(x['erreur'] is None and x['rapport']['statut'] == 'tolerance'
               and x['rapport']['strict'] and x['rapport']['tol'] == 1e-8
               and x['rapport']['paliers_stagnation'] == 0
               and x['rapport']['tolerance_finale_atteinte'] for x in reports)
    assert len(p['barres']) == 3 and all(x['erreur'] is None
             and x['rapport']['statut'] == 'tolerance' for x in p['barres'])
    cable = lit('cable')['cas']
    assert all(x['statut'] == 'echec' and not x['tol_originale_tenue']
               for x in cable['integree']['essais'].values())
    assert all(x['statut'] == 'ok' and x['tol_originale_tenue']
               for x in cable['milieu']['essais'].values())

    # La référence douze directions provient du binaire intermédiaire archivé.
    # Recalculer le courant ; ne pas lui attribuer les anciennes mesures.
    uniform = lit('uniforme')
    assert uniform['dual12']['extension_sha256'] == lit('dynamique-dual12')['extension_sha256']
    for path, digest in uniform['modeles_sha256'].items():
        assert sha(ROOT / path) == digest, path
    kc = modele('uniforme').k_c_m_z()[:2]
    assert np.asarray(uniform['dual12']['KC'], dtype=np.float64).tobytes() == np.asarray(kc, dtype=np.float64).tobytes()
    uniform['courant'] = dict(extension_sha256=current_sha, KC=kc)
    uniform['identique_bit_a_bit'] = True
    ecrit('uniforme', uniform)

    sources = list((ROOT / 'src').rglob('*.rs')) + [ROOT / f for f in (
        'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'python/vinkulum/test_noyau.py',
        'ci/diagnostic_assemblage_local.py', 'ci/bilan_assemblage_local.py',
        'ci/mesure_dynamique_locale.py', 'ci/mesure_statique_alternee.py',
        'ci/valide_assemblage_local.py', 'ci/mesure_statique.py', 'ci/mesure_tangentes.py',
        'ci/confronte_mbdyn.py', 'ci/diagnostic_positions_compensees.py',
        'ci/diagnostic_cable.py', 'ci/local.sh', 'README.md',
        'docs/ASSEMBLAGE_LOCAL.md', 'docs/OBJECTIF_MBDYN.md')]
    artifacts = [f for f in BANCS.glob('assemblage-local-*')
                 if f.name != 'assemblage-local-validation.json']
    manifest = dict(
        date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        temoin_commit='49bdaaa2f249dc008a8e93d6747b3752d8d1644e',
        binaires_sha256=dict(avant=before['extension_sha256'],
            intermediaire_dual12=uniform['dual12']['extension_sha256'], courant=current_sha),
        sources_sha256={str(f.relative_to(ROOT)): sha(f) for f in sorted(sources)},
        artefacts_sha256={str(f.relative_to(ROOT)): sha(f) for f in sorted(artifacts)},
        ci=dict(commande='PY=.venv314/bin/python PYO3_PYTHON=.venv314/bin/python '
                'UV_CACHE_DIR=/tmp/vinkulum-uv-cache ci/local.sh --bancs',
                code_sortie=0, rust=31, python=59, verification=41,
                bancs_rapides=46, contacts=9, api=185),
        controles=dict(statique_alternee=alt_comparison,
            portee_positions_statiques='Trois coordonnées du bout ; équilibre et contraintes '
                                       'vérifiés sur tous les degrés libres.',
            dynamique=dyn, paliers_princeton_tolerance=400, paliers_stagnation=0,
            barres_allongement_exact=3,
            negative_ancienne_extension=dict(tests=2, echecs_attendus=2),
            uniforme_dual9_dual12_identique_bit_a_bit=True,
            cable_integre=dict(originale='echec', lineaire='echec',
                fleche_reference_lineaire_m=cable['integree']['fleche_reference_m'])),
        limites=['Comparaison avant/après Vinkulum ; aucun nouveau classement MBDyn ou Simpack.',
            'G, GG^T, matrices publiques et replis singuliers restent denses.',
            'Blocs de poutre entre rotations encore approchés ; états internes aérodynamiques '
            'figés dans les dérivées.',
            'Coût supplémentaire des jacobiens spatiaux publié dans les mesures dynamiques.',
            'Câble condensé extrême toujours en échec ; validité physique non rétablie.'])
    ecrit('validation', manifest)
    print('Validation : 31 Rust, 59 Python, 400 paliers stricts, empreintes liées aux archives.')


if __name__ == '__main__':
    main()
