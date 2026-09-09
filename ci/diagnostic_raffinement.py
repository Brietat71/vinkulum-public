"""Rejoue les essais linéaires sur les matrices de Newton archivées.

Les calculs SciPy servent au diagnostic algébrique. Le rang au seuil absolu
est publié avec le spectre, sans le
confondre avec un certificat de rang des données.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import platform
import tarfile

import numpy as np
import scipy
from scipy import linalg, sparse
from scipy.sparse.linalg import lsmr, splu


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lit(archive, iteration):
    with tarfile.open(archive) as tar:
        def ouvre(name):
            return io.StringIO(tar.extractfile(name).read().decode())
        f = ouvre(f'kkt-{iteration}.txt')
        m, n, nz = map(int, f.readline().split())
        t = np.array([list(map(float, f.readline().split())) for _ in range(nz)])
        b = np.array([float(x) for x in f])
        a = sparse.csc_array((t[:, 2], (t[:, 0].astype(int), t[:, 1].astype(int))), shape=(m, n))
        f = ouvre(f'echelle-{iteration}.txt')
        nq = int(f.readline())
        d = np.array([float(x) for x in f])
        assert len(b) == m == n == len(d) and 0 < nq < n
        return a, b, d, nq


def spectre(a):
    s = linalg.svdvals(a.toarray())
    keep = s > 1e-12
    return dict(coefficient_max=float(abs(a).max()), norme_2=float(s[0]),
                asymetrie_max=float(abs(a-a.T).max()), seuil_absolu=1e-12,
                rang_au_seuil=int(sum(keep)), sigma_min_retenue=float(s[keep][-1]),
                valeurs_singulieres=s.tolist())


def krylov(a, b):
    r = lsmr(a, b, atol=1e-14, btol=1e-14, conlim=1e14, maxiter=8*a.shape[0])
    return dict(code_arret=r[1], iterations=r[2],
                residu_recalcule=float(np.linalg.norm(b-a@r[0])),
                residu_normal_recalcule=float(np.linalg.norm(a.T@(b-a@r[0]))))


def raffine(a, b, nq, shift):
    try:
        auxiliary = a+sparse.diags_array(np.r_[np.zeros(nq), -shift*np.ones(len(b)-nq)])
        lu = splu(auxiliary)
        x = lu.solve(b)
        residuals = []
        for _ in range(12):
            r = b-a@x
            residuals.append(float(np.linalg.norm(r)))
            if residuals[-1] <= 1e-12*np.linalg.norm(b):
                break
            x += lu.solve(r)
        return dict(decalage=shift, residus=residuals, tolerance=1e-12*float(np.linalg.norm(b)),
                    tolerance_atteinte=bool(residuals[-1] <= 1e-12*np.linalg.norm(b)))
    except RuntimeError as e:
        return dict(decalage=shift, erreur=str(e))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jacobi', type=Path, required=True)
    p.add_argument('--blocs', type=Path, required=True)
    p.add_argument('--sortie', type=Path, required=True)
    args = p.parse_args()
    a, b, d, nq = lit(args.jacobi, 0)
    brut = sparse.diags_array(d)@a@sparse.diags_array(d)
    brut_b = d*b
    ab, bb, db, nqb = lit(args.blocs, 0)
    assert nq == nqb and a.shape == ab.shape
    out = dict(programme_sha256=sha(__file__), python=platform.python_version(),
               numpy=np.__version__, scipy=scipy.__version__,
               archives_sha256={str(f): sha(f) for f in (args.jacobi, args.blocs)},
               dimension=len(b), physiques=nq,
               echelles_jacobi=[float(min(d)), float(max(d))],
               echelles_blocs=[float(min(db)), float(max(db))],
               brut=spectre(brut), jacobi=spectre(a), blocs=spectre(ab),
               lsmr_jacobi=krylov(a, b), lsmr_brut=krylov(brut, brut_b),
               raffinement_blocs=[raffine(ab, bb, nq, shift) for shift in (1e-3, 1e-5, 1e-7, 1e-9)])
    args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print('rangs au seuil absolu', {k: out[k]['rang_au_seuil'] for k in ('brut', 'jacobi', 'blocs')})
    print('résidus du raffinement', [(r['decalage'], r.get('residus')) for r in out['raffinement_blocs']])


if __name__ == '__main__':
    main()
