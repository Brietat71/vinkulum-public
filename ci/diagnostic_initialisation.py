"""Initialisation : références analytiques, analyse redondante, décollement.

Chaque worker fonctionne aussi avec la roue témoin 0.7.1. Le temps de
l'appel public exclut le contre-calcul NumPy ; son pic RSS est relevé avant.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys
import time

import numpy as np
from vinkulum import Noyau
import vinkulum._vinkulum as extension

from diagnostic_cascade import modele


def cascade(nb, tourne, appel):
    from scipy.spatial.transform import Rotation
    q = Rotation.from_rotvec([.31, -.47, .22]).as_matrix() if tourne else np.eye(3)
    n, _ = modele(nb, q)
    before = n.etat_precis()
    if appel == 'initialisation':
        n.enregistre_schema()
    result, error = None, None
    start = time.perf_counter()
    try:
        if appel == 'analyse':
            result = n.k_c_m_z()
        else:
            n.simule(0., .001)
            result = n.schema()[0][0]
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
    elapsed = time.perf_counter()-start
    rss = next(int(l.split()[1]) for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('VmHWM:'))
    controls = dict(etat_inchange=n.etat_precis() == before)
    if result is not None:
        alpha, length, mass = np.pi/3, 1./nb, 1./(3*nb)
        S = np.zeros((18*nb, nb))
        d = q @ (length*np.array([-np.sin(alpha), np.cos(alpha), 0.])/2)
        for i in range(nb):
            for j in range(i+1):
                for body in range(3):
                    S[18*i+6*body:18*i+6*body+3, j] = d*(1. if j == i and body != 1 else 2.)
                if j == i:
                    S[18*i+3:18*i+6, j] = q[:, 2]
                    S[18*i+15:18*i+18, j] = q[:, 2]
        masses = np.tile([mass]*3+[mass/12]*3, 3*nb)
        if appel == 'analyse':
            K, C, M, Z, G = map(np.asarray, result)
            kr = 100./nb-mass*9.81*length*(3*(nb-np.arange(nb))-1)*np.sin(alpha)
            controls.update(mobilites=Z.shape[1], erreur_masse=float(abs(M-np.diag(masses)).max()),
                            erreur_tangente_reduite=float(abs(S.T @ K @ S-np.diag(kr)).max()),
                            erreur_amortissement=float(abs(C).max()),
                            erreur_noyau=float(abs(G @ Z).max()),
                            erreur_parametrage=float(abs(G @ S).max()),
                            erreur_orthogonalite=float(abs(Z.T @ Z-np.eye(nb)).max()))
        else:
            f = np.zeros(18*nb)
            for i in range(3*nb):
                f[6*i:6*i+3] = mass*(q @ [0., -9.81, 0.])
            ref = S @ np.linalg.solve(S.T @ (masses[:, None]*S), S.T @ f)
            controls.update(erreur_acceleration=float(abs(np.asarray(result)-ref).max()))
    return dict(cas=appel, cellules=nb, corps=3*nb, tourne=tourne, secondes=elapsed,
                rss_kib=rss, erreur=error, controles=controls)


def courbure():
    results = []
    for length, speed in [(1., 1.), (1e-6, 10.), (1., 1e8), (1e6, 1e8)]:
        n = Noyau([0., 0., 0.])
        n.corps('masse', 1., np.eye(3).ravel().tolist(), [length, 0., 0.], v=[0., speed, 0.])
        n.distance('bielle', None, 0, [0., 0., 0.], [0., 0., 0.])
        n.enregistre_schema()
        n.simule(0., .001)
        a, ref = n.schema()[0][0][0], -speed*speed/length
        results.append(dict(longueur_m=length, vitesse_m_s=speed, acceleration_m_s2=a,
                            reference_m_s2=ref, erreur_relative=abs(a/ref-1)))
    return dict(cas='courbure', resultats=results)


def commande():
    n = Noyau([0., 0., 0.])
    J = np.eye(3).ravel().tolist()
    n.corps('porteur', 1., J, [0., 0., 0.], w=[2., 0., 0.])
    n.corps('arbre', 1., J, [0., 0., 0.], w=[2., 0., 3.])
    n.liaison('x', None, 0, cible_r=([1., 0., 0.], ('lineaire', [0., 2.])))
    n.liaison('z', 0, 1, cible_r=([0., 0., 1.], ('lineaire', [0., 3.])))
    n.enregistre_schema()
    n.simule(0., .001)
    a = n.schema()[0][0]
    ref = np.zeros(12); ref[10] = -6.
    return dict(cas='commande', acceleration=a, reference=ref.tolist(), erreur=float(abs(np.asarray(a)-ref).max()))


def decollement():
    results = []
    J = np.eye(3).ravel().tolist()
    for rebond in (False, True):
        n = Noyau([0., 0., -9.81 if rebond else 0.])
        if rebond:
            n.corps('bille', 1., J, [0., 0., .5], v=[0., 0., -1.])
            n.contact('sol', 0, [0., 0., 0.], .5, nonlisse=True, restitution=.8)
            ref = [0., 0., -9.81, 0., 0., 0.]
        else:
            n.corps('a', 1., J, [0., 0., 0.])
            n.corps('b', 1., J, [1., 0., 0.], v=[0., 1., 0.])
            n.contact('tangence', 0, [0., 0., 0.], .5, b=1, pb=[0., 0., 0.], rayon_b=.5, nonlisse=True)
            ref = [0.]*12
        n.enregistre_schema()
        n.simule(.0001, .0001)
        a = n.schema()[-1][0]
        results.append(dict(rebond=rebond, acceleration=a, reference=ref,
                            erreur=float(abs(np.asarray(a)-ref).max()), vitesses=n.etat()[3]))
    return dict(cas='decollement', resultats=results)


def libre(nb):
    n = Noyau([0., 0., -9.81])
    J = np.diag([.5, 1., 2.]).ravel().tolist()
    for i in range(nb):
        n.corps(str(i), 1., J, [float(i), 0., 0.])
    n.enregistre_schema()
    start = time.perf_counter()
    n.simule(0., .001)
    elapsed = time.perf_counter()-start
    rss = next(int(l.split()[1]) for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('VmHWM:'))
    a = np.asarray(n.schema()[0][0]).reshape(nb, 6)
    reference = np.zeros_like(a); reference[:, 2] = -9.81
    return dict(cas='libre', corps=nb, secondes=elapsed, rss_kib=rss, erreur=None,
                controles=dict(erreur_acceleration=float(abs(a-reference).max())))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cas', choices=['analyse', 'initialisation', 'courbure', 'commande', 'decollement', 'libre'])
    parser.add_argument('--cellules', type=int, default=8)
    parser.add_argument('--corps', type=int, default=128)
    parser.add_argument('--tourne', action='store_true')
    args = parser.parse_args()
    if args.cas in ('analyse', 'initialisation'):
        out = cascade(args.cellules, args.tourne, args.cas)
    elif args.cas == 'libre':
        out = libre(args.corps)
    else:
        out = dict(courbure=courbure, commande=commande, decollement=decollement)[args.cas]()
    out.update(version=importlib.metadata.version('vinkulum'),
               python=sys.version, numpy=np.__version__, scipy=importlib.metadata.version('scipy'),
               chemin_binaire=extension.__file__,
               binaire_sha256=hashlib.sha256(Path(extension.__file__).read_bytes()).hexdigest())
    print(json.dumps(out, ensure_ascii=False, allow_nan=False))
