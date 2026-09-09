"""Compare les verdicts statiques habituel et strict, en conservant les échecs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import warnings

import numpy as np
from vinkulum import Noyau, _vinkulum
from confronte_mbdyn import princeton

ROOT = Path(__file__).resolve().parents[1]


def appelle(n, strict, iters=100, tol=1e-8):
    before = n.etat()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', RuntimeWarning)
        try:
            n.statique(tol=tol, iters=iters, paliers_max=1, strict=strict)
            error = None
        except Exception as e:
            error = str(e)
    return {'rapport': n.statique_info(), 'erreur': error,
            'etat_restaure_verifie': n.etat() == before if error is not None else None,
            'avertissements': [str(w.message) for w in caught]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sortie', type=Path, required=True)
    args = parser.parse_args()
    out = {'extension_sha256': hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
           'sources_sha256': {str(p): hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                              for p in map(Path, ('src/lib.rs','src/tangent.rs',
                                                 'ci/diagnostic_statut_statique.py','ci/confronte_mbdyn.py'))},
           'cas': []}
    for strict in (False, True):
        n = Noyau([0.,0.,0.])
        n.corps('libre',1.,np.eye(3).ravel().tolist(),[0.,0.,0.])
        n.effort(0,[5e-7,0.,0.],[0.,0.,0.])
        out['cas'].append(dict(nom='corps_sans_equilibre',strict=strict,
                               appels=[appelle(n,strict,iters=12,tol=1e-12)]))
    for ne in (10,40):
        for rampe in (False,True):
            for strict in (False,True):
                n,ids=princeton(ne)
                case=dict(nom='princeton',intervalles=ne,rampe=rampe,strict=strict,appels=[])
                previous=0.
                for index,load in enumerate([.5*(1-math.cos(math.pi*k/50)) for k in range(1,51)] if rampe else [1.],1):
                    force=8.896*(load-previous)/math.sqrt(2)
                    n.effort(ids[-1],[0.,force,force],[0.,0.,0.])
                    result=appelle(n,strict)
                    result.update(palier_externe=index,charge=load)
                    case['appels'].append(result)
                    if result['erreur'] is not None:break
                    previous=load
                case['charge_finale_acceptee']=previous
                case['position_finale']=n.pose(ids[-1])[0]
                case['stagnations']=sum(x['rapport']['paliers_stagnation'] for x in case['appels'])
                case['calcul_termine']=all(x['erreur'] is None for x in case['appels'])
                case['tolerance_tous_paliers']=all(x['erreur'] is None and x['rapport']['statut']=='tolerance' for x in case['appels'])
                out['cas'].append(case)
                print(ne,rampe,strict,'terminé',case['calcul_termine'],'stagnations',case['stagnations'],flush=True)
    args.sortie.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')


if __name__ == '__main__':main()
