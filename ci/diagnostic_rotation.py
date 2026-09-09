"""Corps libre : référence Euler–quaternion indépendante, invariants et repères."""
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
from vinkulum import Noyau

J = np.diag([1e-3, 2e-3, 3e-3])
R0 = Rotation.from_rotvec([.43, -1.11, 2.07]).as_matrix()
WB0 = np.array([.05, 20., .03])
W0 = R0@WB0
MONDE = Rotation.from_rotvec([.51, .62, -1.3]).as_matrix()
MATIERE = Rotation.from_rotvec([1.82, -.47, .15]).as_matrix()
CADRES = {'initial': (np.eye(3), np.eye(3)), 'spatial': (MONDE, np.eye(3)),
          'materiel': (np.eye(3), MATIERE), 'combine': (MONDE, MATIERE)}


def reference():
    def rhs(t, y):
        pi, q = y[:3], y[3:]
        w = np.linalg.solve(J, pi)
        dq = .5*np.r_[q[3]*w+np.cross(q[:3], w), -q[:3]@w]
        return np.r_[np.cross(pi, w), dq]
    initial = np.r_[J@WB0, Rotation.from_matrix(R0).as_quat()]
    s = solve_ivp(rhs, (0., 1.), initial, method='DOP853', rtol=2e-13, atol=2e-15)
    assert s.success
    r = Rotation.from_quat(s.y[3:, -1]).as_matrix()
    return dict(methode='Euler matériel et quaternion, DOP853', rtol=2e-13, atol=2e-15,
                evaluations=s.nfev, rotation=r.ravel().tolist(), omega=(r@np.linalg.solve(J, s.y[:3, -1])).tolist())


def essai(cadre, duree, pas):
    monde, matiere = CADRES[cadre]
    n = Noyau([0., 0., 0.])
    n.corps('corps', .7, (matiere.T@J@matiere).ravel().tolist(), [0., 0., 0.],
            rot=(monde@R0@matiere).ravel().tolist(), w=(monde@W0).tolist())
    debut = time.perf_counter()
    traj = n.simule_em(duree, pas, tous=max(1, round(duree/pas)//200))
    elapsed = time.perf_counter()-debut
    state = n.etat_precis()
    if not traj or traj[-1][0] != state[0]:
        traj.append(tuple(state[:6]))
    e0, l0 = .5*WB0@J@WB0, R0@J@WB0
    errors = dict(energie_relative=0., moment_relatif=0., orthogonalite=0., determinant=0.)
    canonical = []
    for t, positions, rotations, vitesses, omegas, _ in traj:
        assert np.max(np.abs(positions)) == np.max(np.abs(vitesses)) == 0.
        r = monde.T@np.asarray(rotations[0]).reshape(3, 3)@matiere.T
        w = monde.T@np.asarray(omegas[0])
        wb = r.T@w
        e, l = .5*wb@J@wb, r@J@wb
        values = dict(energie_relative=abs(e-e0)/e0, moment_relatif=np.linalg.norm(l-l0)/np.linalg.norm(l0),
                      orthogonalite=np.linalg.norm(r.T@r-np.eye(3)), determinant=abs(np.linalg.det(r)-1.))
        for k, v in values.items():
            errors[k] = max(errors[k], float(v))
        canonical.append([t, r.ravel().tolist(), w.tolist()])
    return dict(cadre=cadre, duree=duree, pas=pas, secondes=elapsed, controles=errors,
                echantillons=canonical, etat_final=state)


def worker():
    ref = reference()
    samples = [essai(c, t, h) for c in CADRES for t, h in [(1., .004), (1., .002), (1., .001), (20., .001)]]
    rref, wref = np.asarray(ref['rotation']).reshape(3, 3), np.asarray(ref['omega'])
    for s in samples:
        if s['duree'] == 1.:
            _, r, w = s['echantillons'][-1]
            s['erreur_rotation_reference_rad'] = float(Rotation.from_matrix(rref.T@np.asarray(r).reshape(3, 3)).magnitude())
            s['erreur_omega_reference'] = float(np.linalg.norm(np.asarray(w)-wref))
    rss = next(int(l.split()[1]) for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('VmHWM:'))
    return dict(reference=ref, inertie=J.tolist(), rotation_initiale=R0.tolist(), omega_initial=W0.tolist(),
                cadres={k: [a.tolist(), b.tolist()] for k, (a, b) in CADRES.items()}, essais=samples, rss_kib=rss)


if __name__ == '__main__':
    print(json.dumps(worker(), ensure_ascii=False, allow_nan=False))
