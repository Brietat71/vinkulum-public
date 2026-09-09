"""Valide les essais orthogonaux, les échecs restaurés et les contrôles physiques."""
import argparse
import json
from pathlib import Path
import statistics

import numpy as np

from bilan_contraintes import etat, sha, strict
from mesure_svd import ROOT, WORKERS


def controles(sample):
    if 'controles' in sample:
        return sample['controles']
    return {k: v for k, v in sample.items() if (
        k.startswith('erreur_') or k == 'orthogonalite_autocontraintes') and v is not None}


def bilan(data):
    assert data['programme_sha256'] == sha(ROOT/'ci/mesure_svd.py')
    assert data['workers_sha256'] == {f: sha(ROOT/f) for f in WORKERS}
    assert data['repetitions'] == 3 and data['fils_demandes'] == 1
    expected = {(f, n, v, r) for f, n in data['cas'] for v in data['binaires'] for r in range(3)}
    samples = {(s['famille'], s['taille'], s['version'], s['repetition']): s for s in data['mesures']}
    assert samples.keys() == expected and len(samples) == len(data['mesures'])
    failures, timeouts = [], []
    for key, s in samples.items():
        if s['code_sortie'] is None:
            assert s['version'] in ('avant', 'svd')
            assert s['erreur'] == f"délai de {data['delai_processus_s']:g} s dépassé"
            timeouts.append(key)
            continue
        assert s['code_sortie'] == 0, (key, s.get('erreur'))
        assert np.isfinite(s['secondes']) and s['secondes'] > 0 and s['rss_kib'] > 0
        etat(s)
        if s['erreur'] is not None:
            assert s['version'] in ('avant', 'svd'), (key, s['erreur'])
            r = s['rapport']
            assert r['strict'] and r['statut'] == 'echec' and not r['tolerance_finale_atteinte']
            assert r['etat_restaure'] and s['etat'] == s['etat_initial']
            failures.append(key)
            continue
        strict(s['rapport'])
        assert s['reactions'] is not None
        assert all(np.isfinite(values).all() for _, values in s['reactions'])
        checks = controles(s)
        assert checks and all(np.isfinite(v) and v >= 0 for v in checks.values())
        for name, value in checks.items():
            limit = 1e-12 if name in ('orthogonalite_autocontraintes', 'erreur_pose_max') else 1e-8
            if name == 'erreur_reactions_absolue':
                n = s['taille']
                limit = 1e-9*max(9.81*n, 9.81*.2*n*(n-1)/2)
            assert value <= limit, (key, name, value, limit)
    rows = []
    for f, n in data['cas']:
        row = dict(famille=f, taille=n, versions={})
        for v in data['binaires']:
            group = [samples[f, n, v, r] for r in range(3)]
            complete = all(s['code_sortie'] == 0 for s in group)
            times = [s['secondes'] for s in group if s['code_sortie'] == 0]
            checks = [controles(s) for s in group]
            row['versions'][v] = dict(mediane_s=statistics.median(times) if complete else None,
                etendue_s=[min(times), max(times)] if complete else None,
                rss_mio=statistics.median(s['rss_kib'] for s in group)/1024 if complete else None,
                statuts=[s.get('rapport', {}).get('statut', 'delai_depasse') for s in group],
                evaluations=[s.get('rapport', {}).get('evaluations') for s in group],
                tentatives=[s.get('rapport', {}).get('tentatives') for s in group],
                controles_max={k: max(c[k] for c in checks if k in c) for k in sorted(set().union(*checks))})
        versions = list(data['binaires'])
        if len(versions) == 2 and all(samples[f, n, v, r]['erreur'] is None for v in versions for r in range(3)):
            va, vb = versions
            row['gain_medianes'] = row['versions'][va]['mediane_s']/row['versions'][vb]['mediane_s']
            row['gains_par_repetition'] = [samples[f, n, va, r]['secondes']/samples[f, n, vb, r]['secondes'] for r in range(3)]
            errors, reaction_error = [0., 0.], 0.
            for i in range(3):
                for j in range(3):
                    a, b = samples[f, n, va, i], samples[f, n, vb, j]
                    for k, (x, y) in enumerate(zip(etat(a), etat(b))):
                        errors[k] = max(errors[k], float(np.max(np.abs(x-y))))
                    for (na, ra), (nb, rb) in zip(a['reactions'], b['reactions'], strict=True):
                        assert na == nb
                        reaction_error = max(reaction_error, float(np.max(np.abs(np.asarray(ra)-rb))))
            assert max(errors) <= 1e-8, (f, n, errors)
            row.update(ecart_position_compensee_m=errors[0], ecart_rotation=errors[1], ecart_reactions_SI=reaction_error)
        rows.append(row)
    return dict(essais=len(samples), equilibres_stricts=len(samples)-len(failures)-len(timeouts),
                echecs_restaures=failures, delais_depasses=timeouts, comparaison=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('entree', type=Path)
    p.add_argument('sortie', type=Path)
    args = p.parse_args()
    out = bilan(json.loads(args.entree.read_text()))
    out.update(entree_sha256=sha(args.entree), programme_sha256=sha(__file__))
    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(f"{out['essais']} essais : {out['equilibres_stricts']} équilibres stricts, "
          f"{len(out['echecs_restaures'])} échecs restaurés, {len(out['delais_depasses'])} délais dépassés.")
