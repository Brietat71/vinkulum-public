"""Confronte le cas fabriqué Rust à sa référence algébrique et à LAPACK gesvd."""
import argparse
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import svd

from bilan_contraintes import sha


def valide(data):
    m, n = data['m'], data['n']
    k = min(m, n)
    a = np.asarray(data['a_col']).reshape(m, n, order='F')
    b, reference, y = (np.asarray(data[key]) for key in ('b', 'reference', 'y'))
    assert np.linalg.norm(a.T@y-reference) < 1e-14
    assert np.linalg.norm(a@reference-b) < 1e-14
    u, sigma, vt = svd(a, full_matrices=False, lapack_driver='gesvd')
    keep = sigma > data['seuil']
    x = vt[keep].T@((u[:, keep].T@b)/sigma[keep])
    checks = {}
    for name in ('nalgebra', 'faer'):
        s = data[name]
        uu = np.asarray(s['u_col']).reshape(m, k, order='F')
        vv = np.asarray(s['vt_col']).reshape(k, n, order='F')
        checks[name] = dict(erreur_solution=float(np.linalg.norm(np.asarray(s['x'])-reference)),
            erreur_sigma=float(np.linalg.norm(np.asarray(s['sigma'])-sigma)),
            erreur_reconstruction_relative=float(np.linalg.norm((uu*np.asarray(s['sigma']))@vv-a)/np.linalg.norm(a)),
            orthogonalite_u=float(np.linalg.norm(uu.T@uu-np.eye(k))),
            orthogonalite_v=float(np.linalg.norm(vv@vv.T-np.eye(k))))
    checks['lapack'] = dict(erreur_solution=float(np.linalg.norm(x-reference)))
    assert checks['lapack']['erreur_solution'] < 1e-13
    assert checks['nalgebra']['erreur_solution'] > 1e-4
    assert all(v < 2e-13 for v in checks['faer'].values())
    return dict(rust=data, lapack=dict(pilote='gesvd', sigma=sigma.tolist(), x=x.tolist()),
                controles=checks, numpy=np.__version__, scipy=scipy.__version__)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('entree', type=Path)
    p.add_argument('sortie', type=Path)
    args = p.parse_args()
    out = valide(json.loads(args.entree.read_text()))
    out.update(entree_sha256=sha(args.entree), programme_sha256=sha(__file__))
    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(json.dumps(out['controles'], ensure_ascii=False, indent=2))
