"""Contrôle des équilibres, comparaison et manifeste des réactions redondantes."""
import argparse
import datetime
import json
from pathlib import Path
import re
import shutil
import statistics

import numpy as np
from vinkulum import _vinkulum
from bilan_contraintes import etat, sha, strict

ROOT = Path(__file__).resolve().parents[1]
BANCS = ROOT/'docs/bancs'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--journal-ci', required=True, type=Path)
    a = p.parse_args()
    data = json.loads((BANCS/'redondances-mesures.json').read_text())
    assert data['programme_sha256'] == sha(ROOT/'ci/mesure_redondances.py')
    for f, h in data['workers_sha256'].items():
        assert sha(ROOT/f) == h
    baseline = 'ac19229e9cce57fc0001af0f12cc32ea2f0b0d258b52abbba269267ae252e01d'
    current = sha(_vinkulum.__file__)
    assert data['binaires']['avant']['sha256'] == baseline
    assert data['binaires']['apres']['sha256'] == current
    expected = {(f,n,v,r) for f, ns in dict(boucle=[16,32,64,128,256,512],
        spatiale=[16,32,64,128], double=[8,16,32,64], chaine=[128], rotors=[128],
        articulee=[128]).items() for n in ns for v in ('avant','apres') for r in range(3)}
    expected |= {('grande',n,'apres',r) for n in (1024,2048) for r in range(3)}
    samples = {(s['famille'],s['corps'],s['version'],s['repetition']): s for s in data['mesures']}
    assert len(samples) == len(data['mesures']) == 108 and samples.keys() == expected
    for (family, nb, version, _), s in samples.items():
        assert s['erreur'] is None and s['code_sortie'] == 0, (family,nb,version,s)
        strict(s['rapport'])
        etat(s)
        assert s['secondes'] > 0 and s['rss_kib'] > 0
        if family in ('boucle', 'grande', 'chaine', 'rotors'):
            assert s['erreur_pose_max'] <= 1e-12
            assert s['erreur_equilibre_independant'] <= 1e-8
            if family in ('boucle','grande'):
                assert np.isfinite(s['orthogonalite_autocontraintes'])
                if family == 'boucle':
                    assert s['orthogonalite_autocontraintes'] <= 1e-12
        elif family == 'articulee':
            for key in ('erreur_moment_independant_Nm','erreur_reactions_N','erreur_fermeture_m'):
                assert s[key] <= 1e-8
        else:
            controls = s['controles']
            assert controls and all(np.isfinite(v) for v in controls.values())
            for key, value in controls.items():
                if key == 'orthogonalite_autocontraintes':
                    assert value <= 1e-12, (family,nb,version,key,value)
                elif key == 'erreur_reactions_reference_SI':
                    assert value <= 1e-7, (family,nb,version,key,value)
                else:
                    assert value <= 1e-8, (family,nb,version,key,value)
    comparison = []
    for family, nb in sorted({(f,n) for f,n,_,_ in expected}):
        row = dict(famille=family, corps=nb, versions={})
        for version in ('avant','apres'):
            group = [s for (f,n,v,r),s in samples.items() if (f,n,v)==(family,nb,version)]
            if not group:
                continue
            times = [s['secondes'] for s in group]
            controls = [dict(s.get('controles',{}), **{k:v for k,v in s.items() if (
                k.startswith('erreur_') or k == 'orthogonalite_autocontraintes') and v is not None}) for s in group]
            row['versions'][version] = dict(mediane_s=statistics.median(times),
                etendue_s=[min(times),max(times)], rss_mio=statistics.median(s['rss_kib'] for s in group)/1024,
                controles_max={k:max(c[k] for c in controls) for k in controls[0]},
                evaluations=[s['rapport']['evaluations'] for s in group])
        if family != 'grande':
            row['gain'] = row['versions']['avant']['mediane_s']/row['versions']['apres']['mediane_s']
            errors = [0.,0.]
            for i in range(3):
                for j in range(3):
                    x,y = [etat(samples[family,nb,v,r]) for v,r in (('avant',i),('apres',j))]
                    for k in range(2):
                        errors[k] = max(errors[k],float(np.max(np.abs(x[k]-y[k]))))
            assert max(errors) <= 1e-10, (family,nb,errors)
            row['ecart_position_compensee_m'],row['ecart_rotation'] = errors
        comparison.append(row)
    journal = a.journal_ci.read_text()
    for phrase in ('38 passed; 0 failed','41/41 cas','46/46 cas','9/9 cas','185 entrées','CI locale OK'):
        assert phrase in journal, phrase
    folder = re.findall(r'journaux : ([^\n]+)',journal)[0]
    tests = [f for f in Path(folder).glob('*.log') if 'Ran 61 tests' in f.read_text()]
    assert len(tests)==1 and '\nOK\n' in tests[0].read_text()
    shutil.copyfile(a.journal_ci,BANCS/'redondances-ci.log')
    shutil.copyfile(tests[0],BANCS/'redondances-tests-python.log')
    princeton = json.loads((BANCS/'redondances-princeton-strict.json').read_text())
    assert princeton['extension_sha256']==current
    assert len(princeton['princeton'])==8 and {(c['intervalles'],c['formulation']) for c in princeton['princeton']}=={
        (n,f) for n in (10,20,40,60) for f in ('milieu','integree')}
    reports = [r for c in princeton['princeton'] for r in c['appels']]
    assert len(reports)==400 and all(c['tolerance_tous_paliers'] for c in princeton['princeton'])
    for r in reports:
        assert r['erreur'] is None
        strict(r['rapport'])
    assert len(princeton['barres'])==3 and all(r['erreur'] is None and
        r['rapport']['statut']=='tolerance' and r['rapport']['strict'] for r in princeton['barres'])
    cable = json.loads((BANCS/'redondances-cable.json').read_text())
    assert cable['extension_sha256']==current
    assert all(r['statut']=='echec' and not r['tol_originale_tenue']
               for r in cable['cas']['integree']['essais'].values())
    assert all(r['statut']=='ok' and r['tol_originale_tenue']
               for r in cable['cas']['milieu']['essais'].values())
    intermediate = json.loads((BANCS/'redondances-intermediaire.json').read_text())
    assert intermediate['programme_sha256']==data['programme_sha256']
    assert intermediate['workers_sha256']==data['workers_sha256']
    assert intermediate['binaires']['apres']['sha256']=='5c89332696f5455f7e61c1462277f06129eb5a8ef301ffd71aed481b7906a8fa'
    result = dict(entree_sha256=sha(BANCS/'redondances-mesures.json'), comparaison=comparison)
    (BANCS/'redondances-comparaison.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    files = list((ROOT/'src').rglob('*.rs'))+[ROOT/f for f in ('Cargo.toml','Cargo.lock',
        'pyproject.toml','README.md','python/vinkulum/test_noyau.py','ci/diagnostic_redondances.py',
        'ci/diagnostic_contraintes.py','ci/diagnostic_chaine_articulee.py','ci/mesure_redondances.py',
        'ci/valide_redondances.py','ci/bilan_contraintes.py','ci/diagnostic_positions_compensees.py',
        'ci/diagnostic_precision_statique.py','ci/diagnostic_statut_statique.py','ci/confronte_mbdyn.py',
        'ci/diagnostic_cable.py','ci/diagnostic_poutres.py','python/vinkulum/campagne.py',
        'ci/local.sh','docs/REDONDANCES_CREUSES.md','docs/OBJECTIF_MBDYN.md')]
    artifacts = [f for f in BANCS.glob('redondances-*') if f.name!='redondances-validation.json']
    manifest = dict(date_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        temoin_commit='6a09df2079885da8540f4e8e7aaa8aa3969c7c5b',
        binaires_sha256=dict(avant=baseline,courant=current),
        intermediaire_sans_borne_de_rang_sha256=intermediate['binaires']['apres']['sha256'],
        traces_double8=dict(commande='VINKULUM_TRACE=1 RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 '
                            'OMP_NUM_THREADS=1 .venv314/bin/python ci/diagnostic_redondances.py double 8',
                            portee='Diagnostic séparé, hors des mesures alternées avec affinité CPU.',
                            binaires_sha256=dict(avant=baseline,intermediaire=intermediate['binaires']['apres']['sha256'],apres=current)),
        sources_sha256={str(f.relative_to(ROOT)):sha(f) for f in sorted(files)},
        artefacts_sha256={str(f.relative_to(ROOT)):sha(f) for f in sorted(artifacts)},
        ci=dict(commande='ci/local.sh --bancs',code_sortie=0,rust=38,python=61,verification=41,
                bancs=46,contacts=9,api=185),
        controles=dict(equilibres_stricts=108,paliers_princeton_stricts=400),
        limites=['Gain de projection de blocs larges de rang plein, entre versions Vinkulum.',
                 'Rang réellement déficient : SVD dense conservée, surcoût du QR tenté publié.',
                 'Newton augmenté singulier reste dense ; autres projecteurs inchangés.',
                 'Norme minimale moins précise à certaines grandes tailles : erreurs publiées.',
                 'Aucune nouvelle exécution MBDyn ou Simpack.'])
    (BANCS/'redondances-validation.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('108 équilibres, 400 paliers stricts ; 38 Rust et 61 Python ; manifeste écrit.')


if __name__ == '__main__':
    main()
