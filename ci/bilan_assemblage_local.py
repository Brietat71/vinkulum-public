"""Vérifie et résume les deux campagnes d'assemblage, sans écarter les échecs."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics

import numpy as np


def bilan(avant, apres):
    for key in ('programme_sha256', 'modeles_et_controles_sha256', 'fils'):
        if avant[key] != apres[key]:
            raise ValueError(f'Protocoles différents : {key}')
    if any(v != '1' for v in apres['fils'].values()):
        raise ValueError('Les mesures doivent demander un seul fil')
    a = {(c['famille'], c['corps']): c for c in avant['mesures']}
    b = {(c['famille'], c['corps']): c for c in apres['mesures']}
    if a.keys() != b.keys() or not a or len(a) != len(avant['mesures']) or len(b) != len(apres['mesures']):
        raise ValueError('Configurations de tangentes différentes ou absentes')
    analyses = []
    for (famille, nb), old in a.items():
        new = b[famille, nb]
        errors = {}
        for nom in ('K', 'C'):
            x, y = (np.array(c['projections'][nom]) for c in (old, new))
            errors[nom] = float(np.linalg.norm(x-y)/max(np.linalg.norm(x), 1e-300))
            if errors[nom] > 1e-5:
                raise ValueError(f'Projection de {nom} incompatible : {famille}, {nb}')
        analyses.append(dict(famille=famille, corps=nb, avant_s=old['mediane_s'],
                              apres_s=new['mediane_s'], gain=old['mediane_s']/new['mediane_s'],
                              ecarts_projections=errors))
    tailles = sorted({c['intervalles'] for c in avant['statique']})
    if tailles != sorted({c['intervalles'] for c in apres['statique']}) or not tailles:
        raise ValueError('Configurations statiques différentes ou absentes')
    static = []
    for ne in tailles:
        groupes = [[c for c in r['statique'] if c['intervalles'] == ne] for r in (avant, apres)]
        if any(len(g) < 3 for g in groupes):
            raise ValueError('Au moins trois répétitions statiques requises')
        for g in groupes:
            if len({c['repetition'] for c in g}) != len(g):
                raise ValueError('Répétition dupliquée')
            for c in g:
                if c['erreur'] is not None or c['code_sortie'] != 0:
                    raise ValueError(f'Calcul statique échoué : {ne}')
                rapports = c['statique_rapports']
                if not rapports or not all(r['statut'] == 'tolerance' and r['strict']
                                            and r['tol'] == 1e-8 for r in rapports):
                    raise ValueError(f'Tolérance stricte non atteinte : {ne}')
                if not np.isfinite([c['seconds'], c['rss_max_kib'], c['free_residual'], c['phi'],
                                    *c['position']]).all():
                    raise ValueError(f'Sortie non finie : {ne}')
                if c['free_residual'] > rapports[-1]['tol'] * rapports[-1]['echelle_force'] + 1e-12:
                    raise ValueError(f'Résidu physique final hors tolérance : {ne}')
                if c['phi'] > rapports[-1]['tol']:
                    raise ValueError(f'Contrainte finale hors tolérance : {ne}')
        times = [statistics.median(c['seconds'] for c in g) for g in groupes]
        rss = [statistics.median(c['rss_max_kib'] for c in g)/1024 for g in groupes]
        pos = [np.array([c['position'] for c in g]) for g in groupes]
        ecart = float(np.max(np.abs(pos[0][:, None, :]-pos[1][None, :, :])))
        if ecart > 1e-8:
            raise ValueError(f'Équilibres différents : {ne}, {ecart}')
        static.append(dict(intervalles=ne, repetitions=[len(g) for g in groupes],
                           avant_s=times[0], apres_s=times[1], gain=times[0]/times[1],
                           rss_avant_mio=rss[0], rss_apres_mio=rss[1], ecart_positions_m=ecart,
                           pire_residu_libre_N=max(c['free_residual'] for c in groupes[1]),
                           pire_contrainte=max(c['phi'] for c in groupes[1])))
    return dict(analyses=analyses, statique=static,
                inflows=dict(avant=avant['inflows'], apres=apres['inflows']),
                roulement=dict(avant=avant['roulement'], apres=apres['roulement']),
                precharge=dict(avant=avant['amortissement_sous_precharge'],
                               apres=apres['amortissement_sous_precharge']))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('avant', type=Path)
    p.add_argument('apres', type=Path)
    p.add_argument('sortie', type=Path)
    args = p.parse_args()
    r = bilan(json.loads(args.avant.read_text()), json.loads(args.apres.read_text()))
    r['entrees_sha256'] = {str(f): hashlib.sha256(f.read_bytes()).hexdigest()
                           for f in (args.avant, args.apres)}
    args.sortie.write_text(json.dumps(r, ensure_ascii=False, indent=2)+'\n')
    for c in r['statique']:
        print(c['intervalles'], f"{c['gain']:.2f}×", f"{c['rss_avant_mio']:.1f} → {c['rss_apres_mio']:.1f} Mio")


if __name__ == '__main__':
    main()
