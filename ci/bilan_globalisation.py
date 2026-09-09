"""Relit la campagne de globalisation sans relancer de solveur.

Les deux échecs de précision identifiés pendant l'exploration restent des
échecs numériques. Ce contrôle exige leur présence et leur restauration ;
il ne les convertit pas en équilibres obtenus.
"""
import argparse
import json
import math
from pathlib import Path

from mesure_globalisation import bilan, corpus


def verifie(archive):
    meta = archive['metadata']
    assert meta['repetitions'] >= 3, 'Au moins trois répétitions pour la campagne finale.'
    assert [c['spec'] for c in archive['cas']] == corpus(), 'Corpus incomplet ou modifié.'
    assert meta['binaires']['reference']['sha256'] != meta['binaires']['candidat']['sha256']
    for c in archive['cas']:
        assert set(c['versions']) == {'reference', 'candidat'}
        for samples in c['versions'].values():
            assert len(samples) == meta['repetitions']
            for s in samples:
                assert 'worker_error' not in s, s
                assert s['spec'] == c['spec']
                assert math.isfinite(s['temps_s']) and s['temps_s'] > 0
                assert s['rapports']
                if s['erreur']:
                    assert s['restaure'] is True, 'Échec sans restauration.'
                for r in s['rapports']:
                    assert r['strict'] is True
                    assert r['evaluations'] >= r['tentatives'] >= 1
    rows = bilan(archive['cas'])
    assert rows == archive['bilan'], 'Bilan stocké différent du recalcul.'
    failed = []
    for row in rows:
        spec = row['spec']
        limit = spec['famille'] == 'poutre' and spec.get('unites') == [1e-3, 1.]
        for version in ('reference', 'candidat'):
            assert row['versions'][version]['controles_reussis'] != limit, row
        assert row['etats_equivalents'] is not False, row
        if limit:
            failed.append(spec)
    assert len(failed) == 2
    return dict(configurations=len(rows), repetitions=meta['repetitions'],
        resolutions_chronometrees=len(rows)*2*meta['repetitions'],
        controles_reussis_par_version=len(rows)-len(failed),
        echecs_numeriques_communs=failed,
        regressions_de_convergence_observees=0,
        ecart_position_max=max(e['position'] for r in rows for e in r['ecarts_etats']),
        ecart_rotation_max=max(e['rotation'] for r in rows for e in r['ecarts_etats']),
        resultats=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive', type=Path)
    p.add_argument('--sortie', type=Path)
    a = p.parse_args()
    result = verifie(json.loads(a.archive.read_text()))
    if a.sortie:
        a.sortie.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'resultats'}, ensure_ascii=False, indent=2))
