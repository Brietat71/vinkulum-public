"""Vérifie la physique et compare les coûts de la campagne alternée."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def strict(r):
    assert r['strict'] and r['statut'] == 'tolerance' and r['tol'] == 1e-8, r
    assert r['tolerance_finale_atteinte'] and r['paliers_stagnation'] == 0, r
    assert r['residu_relatif'] <= r['tol'] and r['contraintes'] <= r['tol'], r


def etat(r):
    s = r['etat']
    assert len(s) == 7
    position = np.asarray(s[1], dtype=np.longdouble)+np.asarray(s[6], dtype=np.longdouble)
    rotation = np.asarray(s[2], dtype=np.float64)
    for a in (position, rotation, np.asarray(s[3]), np.asarray(s[4])):
        assert np.isfinite(a).all()
    assert np.max(np.abs(s[3])) == np.max(np.abs(s[4])) == 0
    return position, rotation


def bilan(data):
    assert data['repetitions'] == 3 and data['fils_demandes'] == 1
    assert data['programme_sha256'] == sha(ROOT/'ci/mesure_contraintes_alternee.py')
    for f, digest in data['workers_sha256'].items():
        assert digest == sha(ROOT/f), f
    cases = {(f, n) for f in ('chaine', 'rotors') for n in (32, 64, 128, 256, 512)} | {
             (f, n) for f in ('boucle', 'articulee') for n in (8, 16, 32, 64, 128, 256)}
    expected = {(f, n, v, r) for f, n in cases for v in ('avant', 'apres') for r in range(3)}
    samples = {(r['famille'], r['corps'], r['version'], r['repetition']): r for r in data['mesures']}
    assert len(samples) == len(data['mesures']) and samples.keys() == expected
    for (family, nb, version, _), sample in samples.items():
        assert sample['erreur'] is None and sample['code_sortie'] == 0, (family, nb, version, sample)
        strict(sample['rapport'])
        assert np.isfinite(sample['secondes']) and sample['secondes'] > 0 and sample['rss_kib'] > 0
        etat(sample)
        assert np.isfinite(np.asarray([v for _, v in sample['reactions']])).all()
        if family == 'articulee':
            for key in ('erreur_moment_independant_Nm', 'erreur_reactions_N', 'erreur_fermeture_m'):
                assert sample[key] <= 1e-8, (family, nb, version, key, sample[key])
        else:
            assert sample['erreur_pose_max'] <= 1e-12
            assert sample['erreur_equilibre_independant'] <= 1e-8, (family, nb, version, sample)
            if family == 'boucle':
                assert np.isfinite(sample['orthogonalite_autocontraintes'])
                # Le témoin perd des chiffres sur cette répartition indéterminée.
                # Son erreur reste publiée ; l'équilibre à 1e-8 est exigé des deux.
                if version == 'apres':
                    assert sample['orthogonalite_autocontraintes'] <= 1e-12
            if family == 'chaine':
                # Référence de la plus grande composante : moment à l'encastrement, en SI.
                scale = max(9.81*nb, 9.81*.2*nb*(nb-1)/2)
                assert sample['erreur_reactions_absolue']/scale <= 1e-9
    out = []
    for family, nb in sorted(cases):
        groups = [[samples[family, nb, v, rep] for rep in range(3)] for v in ('avant', 'apres')]
        times = [statistics.median(s['secondes'] for s in g) for g in groups]
        errors = [0., 0.]
        exact = True
        reaction_error = 0.
        for a in groups[0]:
            for b in groups[1]:
                ea, eb = etat(a), etat(b)
                exact &= a['etat'] == b['etat']
                for i in range(2):
                    errors[i] = max(errors[i], float(np.max(np.abs(ea[i]-eb[i]))))
                ra, rb = [np.asarray([v for _, v in x['reactions']]) for x in (a, b)]
                reaction_error = max(reaction_error, float(np.max(np.abs(ra-rb))))
        assert max(errors) <= 1e-10, (family, nb, errors)
        controls = ('erreur_moment_independant_Nm', 'erreur_reactions_N', 'erreur_fermeture_m') if (
            family == 'articulee') else ('erreur_equilibre_independant', 'erreur_pose_max',
                                        'erreur_reactions_absolue', 'orthogonalite_autocontraintes')
        out.append(dict(famille=family, corps=nb, avant_s=times[0], apres_s=times[1], gain=times[0]/times[1],
                        etendues_s=[[min(s['secondes'] for s in g), max(s['secondes'] for s in g)] for g in groups],
                        rss_median_mio=[statistics.median(s['rss_kib'] for s in g)/1024 for g in groups],
                        ecart_position_compensee_m=errors[0], ecart_rotation=errors[1],
                        etat_exporte_identique=exact, ecart_reactions_max_SI=reaction_error,
                        evaluations=[[s['rapport']['evaluations'] for s in g] for g in groups],
                        controles_max=[{k: max(s[k] for s in g) for k in controls if g[0][k] is not None}
                                       for g in groups]))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('entree', type=Path)
    p.add_argument('sortie', type=Path)
    args = p.parse_args()
    out = dict(entree_sha256=sha(args.entree), programme_sha256=sha(__file__),
               comparaison=bilan(json.loads(args.entree.read_text())))
    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print('132 équilibres stricts vérifiés ; poses et contrôles physiques comparés.')


if __name__ == '__main__':
    main()
