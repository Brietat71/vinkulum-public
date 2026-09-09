"""Parallélogrammes en cascade : graphe connexe, mobilité N, 3N redondances."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from vinkulum import Noyau


def modele(nb, orientation=None):
    q = np.eye(3) if orientation is None else np.asarray(orientation)
    n = Noyau((q @ [0., -9.81, 0.]).tolist())
    longueur, largeur, masse, k = 1./nb, .25, 1./(3*nb), 100./nb
    arm = longueur*np.array([np.cos(np.pi/3), np.sin(np.pi/3), 0.])/2
    base = np.zeros(3)
    joints = []
    for i in range(nb):
        parent = 3*i-2 if i else None
        poses = [base+arm, base+2*arm+[largeur/2, 0., 0.], base+arm+[largeur, 0., 0.]]
        for j, p in enumerate(poses):
            n.corps(f'b{i}_{j}', masse, (masse*np.eye(3)/12).ravel().tolist(),
                    (q @ p).tolist(), rot=q.ravel().tolist())
        specs = [(parent, 3*i, [-largeur/2, 0., 0.] if i else [0., 0., 0.], -arm),
                 (3*i, 3*i+1, arm, [-largeur/2, 0., 0.]),
                 (3*i+1, 3*i+2, [largeur/2, 0., 0.], arm),
                 (parent, 3*i+2, [largeur/2, 0., 0.] if i else [largeur, 0., 0.], -arm)]
        for j, (a, b, pa, pb) in enumerate(specs):
            name = f'p{i}_{j}'
            n.liaison(name, a, b, pa=list(q @ pa if a is None else pa),
                      ra=(q if a is None else np.eye(3)).ravel().tolist(),
                      bloque_t=[0, 1, 2], bloque_r=[0, 1])
            joints.append((name, a, b, np.asarray(pa), np.asarray(pb)))
        axe = q @ [0., 0., 1.] if parent is None else [0., 0., 1.]
        n.couple(f'k{i}', parent, 3*i, list(axe), ('ressort', [k, 0., 0.]))
        base += 2*arm
    return n, joints


def controle(state, reactions, nb, joints, orientation=None):
    from scipy.optimize import brentq
    longueur, largeur, masse, k = 1./nb, .25, 1./(3*nb), 100./nb
    alpha = np.pi/3
    # V = somme_i [mg L (3(N-i)-1) sin(alpha+delta_i) + k delta_i²/2].
    angles = [brentq(lambda d: k*d+masse*9.81*longueur*(3*(nb-i)-1)*np.cos(alpha+d),
                     -.4, 0., xtol=1e-15) for i in range(nb)]
    ref_pos, ref_rot = [], []
    base = np.zeros(3, dtype=np.longdouble)
    for delta in angles:
        arm = longueur*np.array([np.cos(alpha+delta), np.sin(alpha+delta), 0.])/2
        ref_pos.extend([base+arm, base+2*arm+[largeur/2, 0., 0.], base+arm+[largeur, 0., 0.]])
        c, s = np.cos(delta), np.sin(delta)
        r = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
        ref_rot.extend([r, np.eye(3), r])
        base += 2*arm
    pos = np.asarray(state[1], dtype=np.longdouble)+np.asarray(state[6], dtype=np.longdouble)
    rot = np.asarray(state[2]).reshape(3*nb, 3, 3)
    if orientation is not None:
        q = np.asarray(orientation)
        pos = pos @ q.astype(np.longdouble)
        rot = np.einsum('ab,nbc->nac', q.T, rot)
    balance = np.zeros((3*nb, 6), dtype=np.longdouble)
    balance[:, 1] = masse*9.81
    values = dict(reactions)
    gap, off_plane = 0., 0.
    for name, a, b, pa, pb in joints:
        frame = np.eye(3) if a is None else rot[a]
        point_a = pa if a is None else pos[a]+frame@pa
        point_b = pos[b]+rot[b]@pb
        gap = max(gap, float(np.max(np.abs(point_a-point_b))))
        value = np.asarray(values[name])
        off_plane = max(off_plane, float(np.max(np.abs(value[2:]))))
        force = frame@value[:3]
        for body, sign, point in [(a, -1., point_a), (b, 1., point_b)]:
            if body is not None:
                balance[body, :3] += sign*force
                balance[body, 3:] += sign*np.cross(point-pos[body], force)
    for i in range(nb):
        parent = 3*i-2 if i else None
        theta_a = 0. if parent is None else np.arctan2(rot[parent, 1, 0], rot[parent, 0, 0])
        theta_b = np.arctan2(rot[3*i, 1, 0], rot[3*i, 0, 0])
        torque = -k*(theta_b-theta_a)
        balance[3*i, 5] -= torque
        if parent is not None:
            balance[parent, 5] += torque
    return dict(erreur_position_reference_m=float(np.max(np.abs(pos-ref_pos))),
                erreur_rotation_reference=float(np.max(np.abs(rot-ref_rot))),
                erreur_fermeture_m=gap, reactions_hors_plan=off_plane,
                erreur_forces_N=float(np.max(np.abs(balance[:, :3]))),
                erreur_moments_Nm=float(np.max(np.abs(balance[:, 3:]))))


def worker(nb, orientation=None):
    n, joints = modele(nb, orientation)
    initial, error = n.etat_precis(), None
    start = time.perf_counter()
    try:
        n.statique(strict=True, tol=1e-8, iters=100, paliers_max=1)
    except Exception as e:
        error = f'{type(e).__name__}: {e}'
    elapsed = time.perf_counter()-start
    rss = next(int(l.split()[1]) for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('VmHWM:'))
    state, reactions, controls = n.etat_precis(), None, {}
    try:
        reactions = n.reactions()
        controls = controle(state, reactions, nb, joints, orientation)
    except RuntimeError:
        pass
    return dict(famille='cascade', cellules=nb, corps=3*nb, secondes=elapsed, rss_kib=rss,
                erreur=error, rapport=n.statique_info(), etat=state, etat_initial=initial,
                reactions=reactions, controles=controls,
                orientation=None if orientation is None else np.asarray(orientation).tolist())


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cellules', type=int)
    p.add_argument('--repere-tourne', action='store_true')
    a = p.parse_args()
    q = None
    if a.repere_tourne:
        from scipy.spatial.transform import Rotation
        q = Rotation.from_rotvec([.31, -.47, .22]).as_matrix()
    print(json.dumps(worker(a.cellules, q), ensure_ascii=False, allow_nan=False))
