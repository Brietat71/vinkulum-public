"""Contrôle un spectre fabriqué à cinq échelles contre LAPACK gesvd."""
import argparse
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import svd
from scipy.spatial.transform import Rotation

from bilan_contraintes import sha


def valide(data):
    assert [s['echelle'] for s in data] == [1e-300, 1e-200, 1., 1e200, 1e300]
    q = Rotation.from_rotvec([.6, -.7, 1.8]).as_matrix()
    p = Rotation.from_rotvec([-1.2, .4, .8]).as_matrix()
    reference = q@p.T
    rows = []
    for s in data:
        scale = s['echelle']
        original = np.asarray(s['a_col']).reshape(3, 3, order='F')
        a = original/scale
        assert np.linalg.norm(a-q@np.diag([2., 1., -.5])@p.T) < 3e-15
        u, sigma, vt = svd(original, lapack_driver='gesvd')
        lapack = dict(sigma=(sigma/scale).tolist(),
                      erreur_spectre=float(np.linalg.norm(sigma/scale-[2., 1., .5])),
                      reconstruction_relative=float(np.linalg.norm((u*(sigma/scale))@vt-a)/np.linalg.norm(a)))
        assert max(lapack['erreur_spectre'], lapack['reconstruction_relative']) < 3e-14
        methods = {}
        for name, factors in [('brute', s['brute']), ('normalisee', s)]:
            if factors.get('erreur') is not None:
                methods[name] = dict(erreur=factors['erreur'])
                continue
            uu = np.asarray(factors['u_col']).reshape(3, 3, order='F')
            vv = np.asarray(factors['vt_col']).reshape(3, 3, order='F')
            sigma = np.asarray(factors['sigma'])/scale
            assert np.all(np.diff(sigma) <= 0.) and np.min(sigma) >= 0.
            proper = uu@np.diag([1., 1., np.linalg.det(uu@vv)])@vv
            checks = dict(erreur_spectre=float(np.linalg.norm(sigma-[2., 1., .5])),
                          ecart_lapack=float(np.linalg.norm(sigma-lapack['sigma'])),
                          reconstruction_relative=float(np.linalg.norm((uu*sigma)@vv-a)/np.linalg.norm(a)),
                          orthogonalite_u=float(np.linalg.norm(uu.T@uu-np.eye(3))),
                          orthogonalite_v=float(np.linalg.norm(vv@vv.T-np.eye(3))),
                          erreur_rotation=float(np.linalg.norm(proper-reference)))
            assert all(np.isfinite(v) and v >= 0. for v in checks.values())
            if name == 'normalisee':
                assert max(checks.values()) < 3e-14
            methods[name] = dict(sigma=sigma.tolist(), controles=checks)
        rows.append(dict(echelle=scale, lapack=lapack, **methods))
    assert all(r['brute']['controles']['erreur_rotation'] > 1. for r in rows[:2])
    assert all(r['brute']['erreur'] == 'NoConvergence' for r in rows[3:])
    return dict(numpy=np.__version__, scipy=scipy.__version__, pilote='gesvd', cas=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('entree', type=Path)
    p.add_argument('sortie', type=Path)
    args = p.parse_args()
    out = valide(json.loads(args.entree.read_text()))
    out.update(entree_sha256=sha(args.entree), programme_sha256=sha(__file__))
    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print('Cinq échelles validées ; deux spectres bruts erronés et deux échecs bruts conservés.')
