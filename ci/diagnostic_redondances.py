"""Boucles spatiales et chaînes articulées redondantes : contrôles de torseurs."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from vinkulum import Noyau
from diagnostic_contraintes import worker as boucle
from diagnostic_chaine_articulee import modele as chaine_articulee


def rotation(i):
    v = np.array([.2*np.sin(i), .3*np.cos(.7*i), .1*np.sin(.3*i)])
    a = np.linalg.norm(v)
    v /= a
    k = np.array([[0., -v[2], v[1]], [v[2], 0., -v[0]], [-v[1], v[0], 0.]])
    return np.eye(3)+np.sin(a)*k+(1-np.cos(a))*(k@k)


def spatiale(nb):
    gravity = np.array([1.2, -2.4, -9.81])
    n = Noyau(gravity.tolist())
    poses, rotations, masses = [], [], []
    for i in range(nb):
        pos = [.2*i, .11*np.sin(.7*i), .07*np.cos(.2*i)]
        r = rotation(i)
        mass = 1.+.05*np.cos(i)
        n.corps(str(i), mass, np.eye(3).ravel().tolist(), pos, rot=r.ravel().tolist())
        poses.append(pos)
        rotations.append(r)
        masses.append(mass)
        n.liaison('j'+str(i), i-1 if i else None, i,
                  pa=[0., 0., 0.] if i else pos, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    n.liaison('fermeture', None, nb-1, pa=poses[-1], bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    return n, (np.asarray(poses), np.asarray(rotations), np.asarray(masses), gravity)


def controle_spatial(data, reactions):
    pos, rotations, masses, gravity = data
    nb = len(pos)
    points = np.vstack((pos[0], pos[:-1], pos[-1]))
    frames = np.concatenate((np.eye(3)[None], rotations[:-1], np.eye(3)[None]))
    loads = masses[:, None]*gravity
    values = np.asarray([v for _, v in reactions])
    # Une solution particulière est la chaîne encastrée seulement au départ.
    particular = np.zeros((nb+1, 6))
    null = np.zeros((6*(nb+1), 6))
    for i in range(nb):
        particular[i, :3] = frames[i].T@np.sum(loads[i:], axis=0)
        particular[i, 3:] = frames[i].T@np.sum(np.cross(pos[i:]-points[i], loads[i:]), axis=0)
        for k in range(6):
            w = np.eye(6)[k]
            null[6*i:6*i+3, k] = -frames[i].T@w[:3]
            null[6*i+3:6*i+6, k] = frames[i].T@(-w[3:]-np.cross(pos[-1]-points[i], w[:3]))
    null[-6:] = np.eye(6)
    # Seules six inconnues ici, indépendamment du nombre de corps.
    reference = particular.ravel()-null@np.linalg.solve(null.T@null, null.T@particular.ravel())
    balance = np.zeros((nb, 6))
    for i in range(nb+1):
        f, moment = frames[i]@values[i, :3], frames[i]@values[i, 3:]
        b = i if i < nb else nb-1
        balance[b, :3] += f
        balance[b, 3:] += moment+np.cross(points[i]-pos[b], f)
        if 0 < i < nb:
            balance[i-1, :3] -= f
            balance[i-1, 3:] -= moment+np.cross(points[i]-pos[i-1], f)
    balance[:, :3] -= loads
    return dict(erreur_forces_N=float(np.max(np.abs(balance[:, :3]))),
                erreur_moments_Nm=float(np.max(np.abs(balance[:, 3:]))),
                erreur_reactions_reference_SI=float(np.max(np.abs(values.ravel()-reference))),
                orthogonalite_autocontraintes=float(np.linalg.norm(null.T@values.ravel()) /
                                                   (np.linalg.norm(null)*np.linalg.norm(values))))


def double(nb):
    n = chaine_articulee(nb)
    for i in range(nb):
        n.liaison('bis'+str(i), i-1 if i else None, i,
                  pa=[.5/nb, 0., 0.] if i else [0., 0., 0.],
                  bloque_t=[0, 1, 2], bloque_r=[0, 1])
    return n


def controle_double(state, reactions, nb):
    pos = np.asarray(state[1], dtype=np.longdouble)+np.asarray(state[6], dtype=np.longdouble)
    rot = np.asarray(state[2]).reshape(nb, 3, 3)
    angles = np.arctan2(rot[:, 1, 0], rot[:, 0, 0])
    rel = angles-np.r_[0., angles[:-1]]
    points = np.zeros((nb, 3), dtype=np.longdouble)
    points[1:] = pos[:-1]+rot[:-1, :, 0]/(2*nb)
    values = dict(reactions)
    errors, differences = [], []
    for i in range(nb):
        frame = np.eye(3) if i == 0 else rot[i-1]
        ref = np.r_[frame.T@np.array([0., -9.81*(nb-i)/nb, 0.]), 0., 0.]
        a, b = np.asarray(values['j'+str(i)]), np.asarray(values['bis'+str(i)])
        errors.append(float(np.max(np.abs(a+b-ref))))
        differences.append(float(np.max(np.abs(a-b))))
    torques = np.array([9.81/nb*np.sum(pos[i:, 0]-points[i, 0]) for i in range(nb)])
    return dict(erreur_reactions_N=max(errors), ecart_doublons=max(differences),
                erreur_moments_Nm=float(np.max(np.abs(-10.*nb*rel-torques))),
                erreur_fermeture_m=float(np.max(np.abs(pos-rot[:, :, 0]/(2*nb)-points))))


def worker(family, nb):
    if family in ('boucle', 'chaine', 'rotors'):
        return boucle(family, nb)
    n, data = spatiale(nb) if family == 'spatiale' else (double(nb), None)
    initial = n.etat_precis()
    error = None
    start = time.perf_counter()
    try:
        n.statique(strict=True, tol=1e-8, iters=100, paliers_max=1)
    except Exception as e:
        error = f'{type(e).__name__}: {e}'
    elapsed = time.perf_counter()-start
    rss = next(int(l.split()[1]) for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('VmHWM:'))
    state = n.etat_precis()
    controls, reactions = {}, None
    try:
        reactions = n.reactions()
        controls = controle_spatial(data, reactions) if data is not None else controle_double(state, reactions, nb)
    except RuntimeError:
        pass
    return dict(famille=family, corps=nb, erreur=error, secondes=elapsed, rss_kib=rss,
                rapport=n.statique_info(), etat=state, etat_initial=initial,
                reactions=reactions, controles=controls)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('famille', choices=('boucle', 'chaine', 'rotors', 'spatiale', 'double'))
    p.add_argument('corps', type=int)
    a = p.parse_args()
    print(json.dumps(worker(a.famille, a.corps), ensure_ascii=False, allow_nan=False))
