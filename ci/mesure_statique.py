"""Mesure isolée du correctif statique ; référence MBDyn de la confrontation.

RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv314/bin/python \
    ci/mesure_statique.py --sortie /tmp/statique.json
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import statistics
import time

import numpy as np
from vinkulum import _vinkulum
from confronte_mbdyn import princeton

ROOT = Path(__file__).resolve().parents[1]


def run(ne, rampe=False, formulation="milieu", strict=False):
    n, ids = princeton(ne, formulation=formulation)
    start = time.perf_counter()
    previous = 0.
    rapports = []
    charges = [.5*(1-np.cos(np.pi*k/50)) for k in range(1,51)] if rampe else [1.]
    for load in charges:
        force = 8.896*(load-previous)/np.sqrt(2)
        n.effort(ids[-1], [0.,force,force], [0.,0.,0.])
        options = dict(strict=True) if strict else {}
        residual, iterations = n.statique(tol=1e-8, iters=100, paliers_max=1, **options)
        rapports.append(n.statique_info() if hasattr(n, 'statique_info') else None)
        previous = load
    elapsed = time.perf_counter()-start
    return dict(seconds=elapsed,position=n.pose(ids[-1])[0],
                residual=residual,free_residual=float(np.max(np.abs(n.residu_statique()[6:]))),
                phi=max(map(abs,n.phi())),iterations_last_attempt=iterations,
                statique_rapports=rapports,
                paliers_stagnation=(sum(r['paliers_stagnation'] for r in rapports)
                                    if all(r is not None for r in rapports) else None))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie',type=Path,required=True)
    p.add_argument('--formulation', choices=('milieu', 'integree'), default='milieu')
    p.add_argument('--strict', action='store_true')
    args = p.parse_args()
    archive=json.loads((ROOT/'docs/bancs/confrontation-mbdyn-0.6.2.json').read_text())
    ref=next(c['final_values'] for c in archive['cases']['princeton']['configurations']
             if c['engine']=='mbdyn' and c['step']==160)
    result=dict(metadata=dict(formulation=args.formulation,strict=args.strict,python=platform.python_version(),
                    extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
                    source_sha256=hashlib.sha256((ROOT/'src/lib.rs').read_bytes()).hexdigest(),
                    sources_sha256={str(p):hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                                    for p in (Path('src/lib.rs'), Path('src/tangent.rs'),
                                              Path('ci/mesure_statique.py'), Path('ci/confronte_mbdyn.py'))}),
                reference=ref,direct={},rampe={})
    for ne in (10,20,40,60):
        run(ne, formulation=args.formulation, strict=args.strict)
        samples=[run(ne, formulation=args.formulation, strict=args.strict) for _ in range(3)]
        for sample in samples:
            sample['error_m']=float(np.max(np.abs(np.array(sample['position'])-ref)))
            assert sample['free_residual'] < 1e-7 and sample['phi'] < 1e-9
        result['direct'][ne]=dict(samples=samples,median_seconds=statistics.median(x['seconds'] for x in samples))
        print('direct',ne,result['direct'][ne]['median_seconds'],samples[0]['error_m'],flush=True)
        args.sortie.write_text(json.dumps(result,indent=2)+'\n')
    # Rejouer le chemin exact qui échouait, sans transformer cette validation
    # unique en une mesure statistique de performance.
    for ne in (10,20,40):
        sample=run(ne,rampe=True,formulation=args.formulation,strict=args.strict)
        sample['error_m']=float(np.max(np.abs(np.array(sample['position'])-ref)))
        direct=result['direct'][ne]['samples'][0]
        assert np.max(np.abs(np.array(sample['position'])-direct['position'])) < 1e-8
        assert sample['free_residual'] < 1e-7 and sample['phi'] < 1e-9
        result['rampe'][ne]=sample
        print('rampe',ne,sample['seconds'],sample['free_residual'],flush=True)
        args.sortie.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
