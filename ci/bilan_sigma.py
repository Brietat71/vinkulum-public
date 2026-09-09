"""Contrôle l'archive sigma et rend les compromis mesurés, sans interpolation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import median


def bilan(path):
    raw = Path(path).read_bytes()
    data = json.loads(raw)
    rows = data['mesures']
    cases = data['references']
    assert data['repetitions'] == 3 and len(cases) == 7 and len(rows) == 252
    expected_sigmas = {0., .665, 1.}
    assert set(data['sigma']) == expected_sigmas
    summary = []
    for name in cases:
        selected = [r for r in rows if r['cas'] == name]
        steps = sorted({r['h'] for r in selected}, reverse=True)
        assert len(steps) == 4
        samples = set()
        for h in steps:
            for sigma in sorted(expected_sigmas):
                group = [r for r in selected if r['h'] == h and r['sigma'] == sigma]
                assert len(group) == 3 and {r['repetition'] for r in group} == {0, 1, 2}
                assert all('echec' not in r for r in group), (name, h, sigma, group)
                assert all(r['adapt'][0] == 0 for r in group), (name, h, sigma, 'rejeu')
                assert all(r['phi_max_final'] < 1e-9 for r in group), (name, h, sigma, 'contrainte')
                samples.update(r['nombre_sorties'] for r in group)
                error = {}
                for quantity in group[0]['erreurs']:
                    values = [r['erreurs'][quantity] for r in group]
                    assert all(math.isfinite(v) and v >= 0. for v in values)
                    assert max(values)-min(values) < 1e-14, (name, h, sigma, 'non déterministe')
                    error[quantity] = max(values)
                times = [r['temps_s'] for r in group]
                assert all(math.isfinite(t) and t > 0. for t in times)
                summary.append(dict(cas=name, h=h, sigma=sigma, erreurs=error,
                                    mediane_s=median(times), min_s=min(times), max_s=max(times),
                                    newton=group[0]['stats'][0]))
        assert len(samples) == 1, (name, 'sorties différentes')
        if 'ecart_h' in cases[name]:
            # La variation de référence doit rester sous 2 % de la plus
            # petite erreur mesurée, dans chacune des quatre grandeurs.
            for key in cases[name]['ecart_h']:
                best = min(r['erreurs'][key] for r in selected)
                assert cases[name]['ecart_h'][key] < .02*best, (name, key, 'référence')
                assert cases[name]['ecart_sigma'][key] < .02*best, (name, key, 'référence sigma')
    # Seuils de lecture de la grille APRÈS exploration, pas des seuils
    # préenregistrés de preuve ni une optimisation continue du pas.
    thresholds = [('affine', {'orientation_rad': 1e-5}),
                  ('libre', {'orientation_rad': 1e-3}),
                  ('plan', {'position_m': 1e-4}),
                  ('toupie', {'orientation_rad': 1e-3}),
                  ('toupie', {'orientation_rad': 1e-3, 'norme_reaction_N': .005}),
                  ('axes_commandes', {'rotation_rad_s': 1e-4}),
                  ('branche_flexible', {'position_m': 1e-6}),
                  ('boucle_spatiale', {'position_m': 1e-4})]
    choices = []
    for name, limits in thresholds:
        results = {}
        for sigma in sorted(expected_sigmas):
            valid = [r for r in summary if r['cas'] == name and r['sigma'] == sigma
                     and all(r['erreurs'][k] <= limit for k, limit in limits.items())]
            results[str(sigma)] = min(valid, key=lambda r: r['mediane_s']) if valid else None
        choices.append(dict(cas=name, seuils=limits, choix=results))
    return dict(archive_sha256=hashlib.sha256(raw).hexdigest(), essais=252,
                echec=0, configurations=summary, seuils_exploratoires=choices)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('--sortie', type=Path)
    args = parser.parse_args()
    result = bilan(args.archive)
    if args.sortie:
        args.sortie.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print('252 essais contrôlés, 84 configurations ; références et sorties cohérentes.')
