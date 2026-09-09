"""Équilibre non linéaire d'une chaîne à pivots et ressorts, longueur et masse fixes."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from vinkulum import Noyau, _vinkulum


def modele(nb):
    n = Noyau([0., -9.81, 0.])
    length, mass = 1./nb, 1./nb
    for i in range(nb):
        j = mass*length**2/12
        n.corps(str(i), mass, [j, 0., 0., 0., j, 0., 0., 0., j], [(i+.5)*length, 0., 0.])
        a = i-1 if i else None
        n.liaison('j'+str(i), a, i, pa=[length/2, 0., 0.] if i else [0., 0., 0.],
                  bloque_t=[0, 1, 2], bloque_r=[0, 1])
        n.couple('k'+str(i), a, i, [0., 0., 1.], ('ressort', [10.*nb, 0., 0.]))
    return n


def run(nb):
    n = modele(nb)
    start = time.perf_counter()
    error = None
    try:
        n.statique(strict=True, tol=1e-8, iters=100, paliers_max=1)
    except Exception as e:
        error = f'{type(e).__name__}: {e}'
    elapsed = time.perf_counter()-start
    rss = next(int(l.split()[1]) for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('VmHWM:'))
    state = n.etat_precis()
    pos = np.asarray(state[1], dtype=np.longdouble)+np.asarray(state[6], dtype=np.longdouble)
    rot = np.asarray(state[2]).reshape(nb, 3, 3)
    theta = np.arctan2(rot[:, 1, 0], rot[:, 0, 0])
    theta_rel = theta-np.r_[0., theta[:-1]]
    points = np.zeros((nb, 3), dtype=np.longdouble)
    points[1:] = pos[:-1]+rot[:-1, :, 0]/(2*nb)
    gaps = pos-rot[:, :, 0]/(2*nb)-points
    # Le ressort à chaque pivot équilibre le moment de tous les poids en aval.
    torques = np.array([9.81/nb*np.sum(pos[i:, 0]-points[i, 0]) for i in range(nb)])
    err_moment = float(np.max(np.abs(-10.*nb*theta_rel-torques)))
    reactions = None
    err_reactions = None
    try:
        reactions = n.reactions()
        refs = []
        for i in range(nb):
            r = np.eye(3) if i == 0 else rot[i-1]
            refs.append([*(r.T@np.array([0., -9.81*(nb-i)/nb, 0.])), 0., 0.])
        err_reactions = float(np.max(np.abs(np.array([v for _, v in reactions])-refs)))
    except RuntimeError:
        pass
    return dict(corps=nb, secondes=elapsed, rss_kib=rss, erreur=error,
                rapport=n.statique_info(), etat=state, reactions=reactions,
                erreur_moment_independant_Nm=err_moment, erreur_reactions_N=err_reactions,
                erreur_fermeture_m=float(np.max(np.abs(gaps))))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path)
    p.add_argument('--tailles', nargs='+', type=int, default=[8, 16, 32, 64])
    p.add_argument('--repetitions', type=int, default=3)
    p.add_argument('--worker', type=int)
    args = p.parse_args()
    if args.worker is not None:
        print(json.dumps(run(args.worker)))
        return
    if args.sortie is None:
        p.error('--sortie requis')
    cpu = min(os.sched_getaffinity(0))
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    out = dict(cpu=cpu, longueur_m=1., masse_kg=1., raideur_articulation_Nm_rad='10*N',
               extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
               programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), mesures=[])
    for nb in args.tailles:
        for rep in range(args.repetitions):
            cmd = ['taskset', '-c', str(cpu), sys.executable, str(Path(__file__).resolve()), '--worker', str(nb)]
            try:
                c = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=90)
                r = json.loads(c.stdout) if c.returncode == 0 else dict(erreur=c.stderr)
                r.update(code_sortie=c.returncode, stderr=c.stderr)
            except subprocess.TimeoutExpired:
                r = dict(erreur='délai de 90 s dépassé')
            r.update(corps=nb, repetition=rep)
            out['mesures'].append(r)
            args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2)+'\n')
            print(nb, rep, r.get('secondes'), r.get('erreur_moment_independant_Nm'), r['erreur'], flush=True)


if __name__ == '__main__':
    main()
