"""Contrôle les équilibres indépendants avant de comparer temps et mémoire."""
import argparse
import json
from pathlib import Path
import statistics

import numpy as np
from bilan_contraintes import etat, sha, strict

ROOT = Path(__file__).resolve().parents[1]
CONTROLES = {
    'cascade': {'erreur_position_reference_m', 'erreur_rotation_reference', 'erreur_fermeture_m',
                'reactions_hors_plan', 'erreur_forces_N', 'erreur_moments_Nm'},
    'double': {'erreur_reactions_N', 'ecart_doublons', 'erreur_moments_Nm', 'erreur_fermeture_m'},
    'mixte': {'erreur_reactions_N', 'ecart_copies', 'erreur_moments_Nm', 'erreur_fermeture_m'},
    'parallelogrammes': {'erreur_pose_reference_m', 'erreur_rotation_reference', 'erreur_forces_N',
                        'erreur_moments_Nm', 'erreur_reactions_N', 'reactions_hors_plan', 'erreur_fermeture_m'},
    'chaine': {'erreur_reactions_absolue', 'erreur_pose_max', 'erreur_equilibre_independant'},
    'rotors': {'erreur_pose_max', 'erreur_equilibre_independant'},
    'articulee': {'erreur_moment_independant_Nm', 'erreur_reactions_N', 'erreur_fermeture_m'},
    'spatiale': {'erreur_forces_N', 'erreur_moments_Nm', 'erreur_reactions_reference_SI', 'orthogonalite_autocontraintes'},
    'boucle': {'erreur_pose_max', 'erreur_equilibre_independant', 'orthogonalite_autocontraintes'},
    'princeton': {'erreur_reactions_N', 'erreur_moments_Nm', 'ecart_copies'},
}


def controles(s):
    if 'controles' in s:
        return s['controles']
    return {k: v for k, v in s.items() if v is not None and
            (k.startswith('erreur_') or k == 'orthogonalite_autocontraintes')}


def bilan(data, source_root=ROOT):
    assert data['repetitions'] == 3 and data['fils_demandes'] == 1
    assert data['programme_sha256'] == sha(ROOT/'ci/mesure_raffinement.py')
    for f, digest in data['workers_sha256'].items():
        assert sha(ROOT/f) == digest, f
    for f, digest in data['sources_candidat_sha256'].items():
        assert sha(source_root/f) == digest, f
    assert data['versions']['avant']['version'] == '0.7.0'
    assert data['versions']['apres']['version'] == '0.7.1'
    for name in ('numpy', 'scipy', 'python'):
        assert data['versions']['avant'][name] == data['versions']['apres'][name], name
    cases = {(f, n) for f in ('cascade', 'cascade_tournee') for n in (8, 16, 32)}
    cases |= {(f, n) for f in ('double', 'mixte', 'parallelogrammes') for n in (8, 32)}
    cases |= {(f, n) for f in ('chaine', 'rotors', 'articulee', 'spatiale', 'boucle') for n in (32, 128)}
    cases |= {(f, n) for f in ('princeton_milieu', 'princeton_integree') for n in (10, 40)}
    expected = {(f, n, v, r) for f, n in cases for v in ('avant', 'apres') for r in range(3)}
    samples = {(s['famille'], s['taille'], s['version'], s['repetition']): s for s in data['mesures']}
    assert samples.keys() == expected and len(samples) == len(data['mesures']) == 156
    failures = []
    for key, s in samples.items():
        if s.get('erreur') is not None:
            failures.append(dict(cas=key, erreur=s['erreur']))
            assert key[2] == 'avant', (key, s['erreur'])
            if s.get('code_sortie') == 0:
                assert s['rapport']['statut'] == 'echec' and s['rapport']['etat_restaure']
                assert s['etat'] == s['etat_initial']
            continue
        assert s['code_sortie'] == 0 and np.isfinite(s['secondes']) and s['secondes'] > 0 and s['rss_kib'] > 0, key
        strict(s['rapport'])
        etat(s)
        assert s['reactions'] is not None
        assert all(np.isfinite(values).all() for _, values in s['reactions'])
        checks = controles(s)
        assert checks.keys() == CONTROLES[key[0].split('_')[0]], key
        assert all(np.isfinite(x) and x >= 0 for x in checks.values()), key
        for name, value in checks.items():
            tol = 1e-12 if name in ('erreur_pose_max', 'orthogonalite_autocontraintes') else 1e-8
            if name == 'erreur_reactions_absolue':
                n = s['corps']
                tol = 1e-9 * max(9.81*n, 9.81*.2*n*(n-1)/2)
            assert value <= tol, (key, name, value, tol)
    rows = []
    for family, size in sorted(cases):
        row = dict(famille=family, taille=size, versions={})
        groups = {}
        for version in ('avant', 'apres'):
            group = [samples[family, size, version, r] for r in range(3)]
            groups[version] = group
            if any(s.get('erreur') is not None for s in group):
                row['versions'][version] = dict(erreurs=[s.get('erreur') for s in group])
                continue
            times = [s['secondes'] for s in group]
            checks = [controles(s) for s in group]
            row['versions'][version] = dict(mediane_s=statistics.median(times),
                etendue_s=[min(times), max(times)], rss_mio=statistics.median(s['rss_kib'] for s in group)/1024,
                evaluations=[s['rapport']['evaluations'] for s in group],
                tentatives=[s['rapport']['tentatives'] for s in group],
                controles_max={k: max(c[k] for c in checks) for k in checks[0]})
        if all('mediane_s' in row['versions'][v] for v in ('avant', 'apres')):
            row['gain'] = row['versions']['avant']['mediane_s']/row['versions']['apres']['mediane_s']
            errors, reactions = [0., 0.], 0.
            for a in groups['avant']:
                for b in groups['apres']:
                    for i, (x, y) in enumerate(zip(etat(a), etat(b), strict=True)):
                        errors[i] = max(errors[i], float(np.max(np.abs(x-y))))
                    for (na, ra), (nb, rb) in zip(a['reactions'], b['reactions'], strict=True):
                        assert na == nb
                        reactions = max(reactions, float(np.max(np.abs(np.asarray(ra)-rb))))
            row.update(ecart_position_compensee_m=errors[0], ecart_rotation=errors[1],
                       ecart_reactions_SI=reactions)
        rows.append(row)
    return dict(essais=len(samples), equilibres_stricts=len(samples)-len(failures),
                echecs_conserves=failures, comparaisons=rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mesures', type=Path)
    p.add_argument('sortie', type=Path)
    args = p.parse_args()
    result = bilan(json.loads(args.mesures.read_text()))
    result.update(mesures_sha256=sha(args.mesures), validateur_sha256=sha(__file__))
    args.sortie.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(result['equilibres_stricts'], 'équilibres stricts sur', result['essais'])
    for row in result['comparaisons']:
        before, after = row['versions']['avant'], row['versions']['apres']
        if 'gain' in row:
            print(row['famille'], row['taille'], f"{before['mediane_s']:.6f} → {after['mediane_s']:.6f} s",
                  f"×{row['gain']:.3f}", 'Mio', before['rss_mio'], after['rss_mio'])


if __name__ == '__main__':
    main()
