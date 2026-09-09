"""Recalcule invariants, ordre temporel, écarts entre repères et temps mesurés."""
import argparse
import json
from pathlib import Path
import statistics

import numpy as np
from scipy.spatial.transform import Rotation

from bilan_contraintes import sha

ROOT = Path(__file__).resolve().parents[1]
CADRES = ('initial', 'spatial', 'materiel', 'combine')
PAS = ((1., .004), (1., .002), (1., .001), (20., .001))


def controle(run):
    assert run['code_sortie'] == 0 and run['rss_kib'] > 0, run.get('erreur')
    samples = {(s['cadre'], s['duree'], s['pas']): s for s in run['essais']}
    assert len(samples) == len(run['essais']) == 16
    assert samples.keys() == {(c, t, h) for c in CADRES for t, h in PAS}
    j, r0, w0 = (np.asarray(run[k]) for k in ('inertie', 'rotation_initiale', 'omega_initial'))
    wb0 = r0.T@w0
    e0, l0 = .5*wb0@j@wb0, r0@j@wb0
    ref = run['reference']
    assert ref['rtol'] == 2e-13 and ref['atol'] == 2e-15 and ref['evaluations'] > 0
    rref, wref = np.asarray(ref['rotation']).reshape(3, 3), np.asarray(ref['omega'])
    maximum = dict(energie_relative=0., moment_relatif=0., orthogonalite=0., determinant=0.)
    for (frame, duration, h), s in samples.items():
        assert np.isfinite(s['secondes']) and s['secondes'] > 0.
        assert abs(s['etat_final'][0]-duration) < 1e-10
        states = s['echantillons']
        times = np.asarray([x[0] for x in states])
        assert len(states) >= 200 and np.all(np.diff(times) > 0.)
        assert abs(times[-1]-duration) < 1e-10
        r = np.asarray([x[1] for x in states]).reshape(-1, 3, 3)
        w = np.asarray([x[2] for x in states])
        assert np.isfinite(r).all() and np.isfinite(w).all()
        wb = np.einsum('tji,tj->ti', r, w)
        energy = .5*np.einsum('ti,ij,tj->t', wb, j, wb)
        momentum = np.einsum('tij,jk,tk->ti', r, j, wb)
        checks = dict(energie_relative=float(np.max(np.abs(energy-e0)/e0)),
                      moment_relatif=float(np.max(np.linalg.norm(momentum-l0, axis=1))/np.linalg.norm(l0)),
                      orthogonalite=float(np.max(np.linalg.norm(r.transpose(0, 2, 1)@r-np.eye(3), axis=(1, 2)))),
                      determinant=float(np.max(np.abs(np.linalg.det(r)-1.))))
        for key, value in checks.items():
            assert abs(value-s['controles'][key]) < 3e-15
            assert value < (3e-15 if key in ('orthogonalite', 'determinant') else 1e-11), (frame, key, value)
            maximum[key] = max(maximum[key], value)
        left, right = np.asarray(run['cadres'][frame])
        final = s['etat_final']
        assert np.max(np.abs(final[1])) == np.max(np.abs(final[3])) == 0.
        assert np.linalg.norm(left.T@np.asarray(final[2][0]).reshape(3, 3)@right.T-r[-1]) < 1e-14
        assert np.linalg.norm(left.T@np.asarray(final[4][0])-w[-1]) < 1e-13
        if duration == 1.:
            errors = dict(erreur_rotation_reference_rad=float(Rotation.from_matrix(rref.T@r[-1]).magnitude()),
                          erreur_omega_reference=float(np.linalg.norm(w[-1]-wref)))
            for key, value in errors.items():
                assert abs(value-s[key]) < 1e-13
                assert 0. < value < 7000*h*h
    orders, frames = {}, []
    for c in CADRES:
        orders[c] = {}
        for key in ('erreur_rotation_reference_rad', 'erreur_omega_reference'):
            errors = [samples[c, 1., h][key] for h in (.004, .002, .001)]
            values = np.log2(np.asarray(errors[:-1])/errors[1:])
            assert np.all((1.98 < values) & (values < 2.02)), (c, key, values)
            orders[c][key] = dict(erreurs=errors, ordres=values.tolist())
    for duration, h in PAS:
        base = samples['initial', duration, h]['echantillons']
        diffs = np.zeros(2)
        for c in CADRES[1:]:
            states = samples[c, duration, h]['echantillons']
            assert len(states) == len(base)
            for a, b in zip(states, base, strict=True):
                assert a[0] == b[0]
                diffs = np.maximum(diffs, [np.linalg.norm(np.asarray(a[k])-b[k]) for k in (1, 2)])
        assert diffs[0] < 1e-6 and diffs[1] < 1e-5, (duration, h, diffs)
        frames.append(dict(duree=duration, pas=h, ecart_rotation_frobenius=float(diffs[0]),
                           ecart_omega=float(diffs[1])))
    return dict(controles_max=maximum, ordres=orders, ecarts_reperes=frames)


def bilan(data):
    assert data['repetitions'] == 3 and data['fils_demandes'] == 1 and data['delai_processus_s'] == 60
    assert data['programme_sha256'] == sha(ROOT/'ci/mesure_rotations.py')
    assert data['worker_sha256'] == sha(ROOT/'ci/diagnostic_rotation.py')
    runs = {(r['version'], r['repetition']): r for r in data['mesures']}
    assert len(runs) == len(data['mesures']) == len(data['binaires'])*3
    assert runs.keys() == {(v, i) for v in data['binaires'] for i in range(3)}
    checks = {v: [controle(runs[v, i]) for i in range(3)] for v in data['binaires']}
    rows = []
    for c in CADRES:
        for t, h in PAS:
            row = dict(cadre=c, duree=t, pas=h, versions={})
            for v in data['binaires']:
                group = [next(s for s in runs[v, i]['essais'] if (s['cadre'], s['duree'], s['pas']) == (c, t, h))
                         for i in range(3)]
                seconds = [s['secondes'] for s in group]
                row['versions'][v] = dict(mediane_s=statistics.median(seconds),
                    etendue_s=[min(seconds), max(seconds)], secondes=seconds,
                    rss_mio=statistics.median(runs[v, i]['rss_kib'] for i in range(3))/1024)
            first = next(iter(data['binaires']))
            row['gains_medianes'] = {v: row['versions'][first]['mediane_s']/r['mediane_s']
                                     for v, r in row['versions'].items() if v != first}
            rows.append(row)
    return dict(processus=len(runs), simulations=16*len(runs), controles=checks, comparaison=rows,
                portee='Corps libre ; invariants recalculés aux instantanés sauvegardés, référence temporelle à 1 s.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('entree', type=Path)
    p.add_argument('sortie', type=Path)
    args = p.parse_args()
    out = bilan(json.loads(args.entree.read_text()))
    out.update(entree_sha256=sha(args.entree), programme_sha256=sha(__file__))
    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(f"{out['simulations']} simulations : invariants, ordre deux et changements de repères contrôlés.")
