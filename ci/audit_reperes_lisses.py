"""Sépare covariance locale, amplification et raffinement sur la toupie libre.

La sensibilité calculée n'est pas une borne validée de propagation des arrondis.
Elle ne clôt pas l'obligation de certification des trajectoires.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from vinkulum import __version__, _vinkulum
import diagnostic_rotation as d


def skew(w):
    x, y, z = w
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path, required=True)
    a = p.parse_args()
    jinv = np.linalg.inv(d.J)

    def rhs(t, y):
        w = y[:3]
        jac = jinv@(-skew(w)@d.J+skew(d.J@w))
        return np.r_[jinv@np.cross(d.J@w, w), (jac@y[3:].reshape(3, 3)).ravel()]

    sensibilites = []
    for tol in (1e-10, 1e-12, 3e-14):
        ref = solve_ivp(rhs, (0., 20.), np.r_[d.WB0, np.eye(3).ravel()], method='DOP853',
                        rtol=tol, atol=tol/100, t_eval=[.001, 1., 5., 10., 20.])
        if not ref.success or ref.t[-1] != 20. or not np.isfinite(ref.y).all():
            raise RuntimeError('référence variationnelle incomplète')
        sensibilites.append(dict(tolerance=tol, instants=ref.t.tolist(),
                                sensibilites=[s[3:].reshape(3, 3).tolist() for s in ref.y.T],
                                normes=[float(np.linalg.norm(s[3:].reshape(3, 3), 2)) for s in ref.y.T]))
    essais = []
    for t, h in ((.001, .001), (1., .001), (20., .001), (20., .0005), (20., .00025)):
        frames = [d.essai(c, t, h) for c in d.CADRES]
        base = frames[0]['echantillons']
        mesures = []
        for run in frames:
            if len(run['echantillons']) != len(base):
                raise RuntimeError('échantillonnages différents')
            ecarts = []
            for row, target in zip(run['echantillons'], base, strict=True):
                if row[0] != target[0]:
                    raise RuntimeError('instants différents')
                ecarts.append([float(np.linalg.norm(np.asarray(row[k])-target[k])) for k in (1, 2)])
            mesures.append(dict(cadre=run['cadre'], ecarts_max=np.max(ecarts, axis=0).tolist(),
                                ecarts_finaux=ecarts[-1], invariants=run['controles']))
        essais.append(dict(duree=t, pas=h, mesures=mesures))
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    result = dict(schema='vinkulum.diagnostic_reperes.1', version=__version__,
                  extension_sha256=sha(_vinkulum.__file__), programme_sha256=sha(__file__),
                  modele_sha256=sha(d.__file__), sensibilites=sensibilites, essais=essais,
                  amplification_arrondis_certifiee=False, anomalie_long_terme_close=False,
                  interpretation='Diagnostic de sensibilité et de covariance ; pas une attribution exhaustive des écarts.')
    a.sortie.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(__version__, 'diagnostic de repères terminé ; obligation temporelle ouverte.')


if __name__ == '__main__':
    main()
