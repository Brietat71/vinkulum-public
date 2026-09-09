"""Mécanismes mobiles redondants, avec références de pose et de torseurs indépendantes."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from vinkulum import Noyau
from diagnostic_chaine_articulee import modele as chaine, run as articulee
from diagnostic_redondances import worker as redondances


def mixte(nb):
    n = chaine(nb)
    for i in range(nb):
        for j in range(1, 1+i % 3):
            s = -1. if j % 2 else 1.
            n.liaison(f'copie{i}_{j}', i-1 if i else None, i,
                      pa=[.5/nb, 0., 0.] if i else [0., 0., 0.],
                      ra=np.diag([s, s, 1.]).ravel().tolist(),
                      bloque_t=[0, 1, 2], bloque_r=[0, 1])
    return n


def controle_mixte(state, reactions, nb):
    pos = np.asarray(state[1], dtype=np.longdouble)+np.asarray(state[6], dtype=np.longdouble)
    rot = np.asarray(state[2]).reshape(nb, 3, 3)
    angles = np.arctan2(rot[:, 1, 0], rot[:, 0, 0])
    points = np.zeros((nb, 3), dtype=np.longdouble)
    points[1:] = pos[:-1]+rot[:-1, :, 0]/(2*nb)
    values = dict(reactions)
    errors, sharing = [], []
    for i in range(nb):
        frame = np.eye(3) if i == 0 else rot[i-1]
        force = frame.T @ [0., -9.81*(nb-i)/nb, 0.]
        copies = []
        count = 1+i % 3
        for j in range(count):
            s = -1. if j % 2 else 1.
            name = f'copie{i}_{j}' if j else f'j{i}'
            value = np.asarray(values[name])*[s, s, 1., s, s]
            errors.append(float(np.max(np.abs(value-np.r_[force/count, 0., 0.]))))
            copies.append(value)
        sharing.append(float(np.max(np.abs(np.asarray(copies)-copies[0]))))
    torque = np.array([9.81/nb*np.sum(pos[i:, 0]-points[i, 0]) for i in range(nb)])
    return dict(erreur_reactions_N=max(errors), ecart_copies=max(sharing),
                erreur_moments_Nm=float(np.max(np.abs(-10.*nb*(angles-np.r_[0., angles[:-1]])-torque))),
                erreur_fermeture_m=float(np.max(np.abs(pos-rot[:, :, 0]/(2*nb)-points))))


def parallelogrammes(nb):
    # Quatre pivots spatiaux par mécanisme : 20 équations, 18 coordonnées,
    # rang 17, une mobilité. Les trois dépendances ne sont pas des doublons.
    n = Noyau([0., -9.81, 0.])
    alpha = np.pi/3
    x, y = np.cos(alpha), np.sin(alpha)
    for i in range(nb):
        for j, p in enumerate(([x/2, y/2, 0.], [.5+x, y, 0.], [1+x/2, y/2, 0.])):
            p[0] += 3*i
            n.corps(f'b{i}_{j}', 1., np.eye(3).ravel().tolist(), p)
        specs = [(None, 3*i, [3.*i, 0., 0.]),
                 (3*i, 3*i+1, [x/2, y/2, 0.]),
                 (3*i+1, 3*i+2, [.5, 0., 0.]),
                 (None, 3*i+2, [3.*i+1, 0., 0.])]
        for j, (a, b, p) in enumerate(specs):
            n.liaison(f'p{i}_{j}', a, b, pa=p, bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.couple(f'k{i}', None, 3*i, [0., 0., 1.], ('ressort', [100., 0., 0.]))
    return n


def princeton_double(nb, formulation):
    from confronte_mbdyn import princeton
    n, ids = princeton(nb, formulation)
    n.liaison('copie_enc', None, ids[0], bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    n.effort(ids[-1], [0., 8.896/np.sqrt(2), 8.896/np.sqrt(2)], [0., 0., 0.])
    return n


def controle_princeton(state, reactions, nb):
    pos = np.asarray(state[1], dtype=np.longdouble)+np.asarray(state[6], dtype=np.longdouble)
    force = np.array([0., 8.896/np.sqrt(2), 8.896/np.sqrt(2)])
    reference = np.r_[force, np.cross(pos[-1], force)]/2
    values = np.asarray([v for _, v in reactions])
    return dict(erreur_reactions_N=float(np.max(np.abs(values[:, :3]-reference[:3]))),
                erreur_moments_Nm=float(np.max(np.abs(values[:, 3:]-reference[3:]))),
                ecart_copies=float(np.max(np.abs(values[0]-values[1]))))


def controle_parallelogrammes(state, reactions, nb):
    from scipy.optimize import brentq
    alpha = np.pi/3
    # Énergie : V = 2 m g L sin(alpha+delta) + k delta²/2.
    delta = brentq(lambda d: 100*d+2*9.81*np.cos(alpha+d), -.4, 0., xtol=1e-15)
    x, y = np.cos(alpha+delta), np.sin(alpha+delta)
    c, s = np.cos(delta), np.sin(delta)
    r = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    reference_rot = np.array([r, np.eye(3), r])
    pos = np.asarray(state[1], dtype=np.longdouble)+np.asarray(state[6], dtype=np.longdouble)
    rot = np.asarray(state[2]).reshape(3*nb, 3, 3)
    values = dict(reactions)
    controls = dict(erreur_pose_reference_m=0., erreur_rotation_reference=0.,
                    erreur_forces_N=0., erreur_moments_Nm=0., erreur_reactions_N=0.,
                    reactions_hors_plan=0., erreur_fermeture_m=0.)
    for i in range(nb):
        p = pos[3*i:3*i+3]
        rr = rot[3*i:3*i+3]
        ref = np.array([[3*i+x/2, y/2, 0.], [3*i+.5+x, y, 0.], [3*i+1+x/2, y/2, 0.]])
        controls['erreur_pose_reference_m'] = max(controls['erreur_pose_reference_m'], float(np.max(np.abs(p-ref))))
        controls['erreur_rotation_reference'] = max(controls['erreur_rotation_reference'], float(np.max(np.abs(rr-reference_rot))))
        arm = np.array([np.cos(alpha)/2, np.sin(alpha)/2, 0.])
        points = np.array([[3.*i, 0., 0.], p[0]+rr[0]@arm, p[1]+rr[1]@[.5, 0., 0.], [3.*i+1., 0., 0.]])
        homologues = np.array([p[0]-rr[0]@arm, p[1]-rr[1]@[.5, 0., 0.], p[2]+rr[2]@arm, p[2]-rr[2]@arm])
        controls['erreur_fermeture_m'] = max(controls['erreur_fermeture_m'], float(np.max(np.abs(points-homologues))))
        frames = [np.eye(3), rr[0], rr[1], np.eye(3)]
        # Équilibre plan de chaque corps : forces et moments aux points de
        # liaison, sans utiliser les gradients ou tangentes du noyau.
        a = np.zeros((9, 8))
        forces = []
        for j, (ba, bb) in enumerate([(None, 0), (0, 1), (1, 2), (None, 2)]):
            value = np.asarray(values[f'p{i}_{j}'])
            controls['reactions_hors_plan'] = max(controls['reactions_hors_plan'], float(np.max(np.abs(value[2:]))))
            forces.extend((frames[j]@value[:3])[:2])
            for body, sign in [(ba, -1.), (bb, 1.)]:
                if body is not None:
                    dx, dy = (points[j]-p[body])[:2]
                    a[3*body:3*body+3, 2*j:2*j+2] = sign*np.array([[1., 0.], [0., 1.], [-dy, dx]])
        f = np.tile([0., -9.81, 0.], 3)
        f[2] = -100*np.arctan2(rr[0, 1, 0], rr[0, 0, 0])
        residual = (a@forces-f).reshape(3, 3)
        controls['erreur_forces_N'] = max(controls['erreur_forces_N'], float(np.max(np.abs(residual[:, :2]))))
        controls['erreur_moments_Nm'] = max(controls['erreur_moments_Nm'], float(np.max(np.abs(residual[:, 2]))))
        reference_forces = np.linalg.lstsq(a, f, rcond=None)[0]
        controls['erreur_reactions_N'] = max(controls['erreur_reactions_N'], float(np.max(np.abs(reference_forces-forces))))
    return controls


def worker(family, nb):
    if family in ('double', 'boucle', 'spatiale', 'chaine', 'rotors'):
        return redondances(family, nb)
    if family == 'articulee':
        d = articulee(nb)
        d['controles'] = {k: v for k, v in d.items() if k.startswith('erreur_') and v is not None}
        return d
    if family.startswith('princeton_'):
        n, controle = princeton_double(nb, family.split('_')[1]), controle_princeton
    elif family == 'mixte':
        n, controle = mixte(nb), controle_mixte
    else:
        n, controle = parallelogrammes(nb), controle_parallelogrammes
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
        controls = controle(state, reactions, nb)
    except RuntimeError:
        pass
    return dict(famille=family, taille=nb, corps=len(state[1]), secondes=elapsed, rss_kib=rss,
                erreur=error, rapport=n.statique_info(), etat=state, etat_initial=initial,
                reactions=reactions, controles=controls)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('famille', choices=('double', 'mixte', 'parallelogrammes', 'boucle', 'spatiale', 'chaine', 'rotors', 'articulee', 'princeton_milieu', 'princeton_integree'))
    p.add_argument('taille', type=int)
    a = p.parse_args()
    print(json.dumps(worker(a.famille, a.taille), ensure_ascii=False, allow_nan=False))
