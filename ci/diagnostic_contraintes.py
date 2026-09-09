"""Coût et précision des contraintes statiques, avec références analytiques."""
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

J = [1., 0., 0., 0., 1., 0., 0., 0., 1.]


def modele(famille, nb):
    n = Noyau([0., 0., -9.81])
    for i in range(nb):
        n.corps(str(i), 1., J, [.2*i, 0., 0.])
        if famille == 'rotors':
            n.liaison('j'+str(i), None, i, pa=[.2*i, 0., 0.],
                      bloque_t=[0, 1, 2], bloque_r=[0, 1])
            n.couple('k'+str(i), None, i, [0., 0., 1.], ('ressort', [10., 0., .1]))
        else:
            n.liaison('j'+str(i), i-1 if i else None, i,
                      pa=[.2, 0., 0.] if i else [0., 0., 0.],
                      bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    if famille == 'boucle':
        n.liaison('fermeture', None, nb-1, pa=[.2*(nb-1), 0., 0.],
                  bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
    return n


def worker(famille, nb):
    n = modele(famille, nb)
    start = time.perf_counter()
    error = None
    try:
        n.statique(strict=True, tol=1e-8, iters=10, paliers_max=1)
    except Exception as e:
        error = f'{type(e).__name__}: {e}'
    elapsed = time.perf_counter()-start
    rss = next(int(l.split()[1]) for l in Path('/proc/self/status').read_text().splitlines()
               if l.startswith('VmHWM:'))
    state = n.etat_precis()
    try:
        reactions = n.reactions()
    except RuntimeError:
        reactions = []
    if famille == 'chaine':
        reference = [[0., 0., -9.81*(nb-i), 0.,
                      9.81*.2*(nb-i)*(nb-i-1)/2, 0.] for i in range(nb)]
        values = np.array([v for _, v in reactions])
        reaction_error = float(np.max(np.abs(values-reference))) if values.size else None
    else:
        reaction_error = None
    balance_error = selfstress_error = None
    if reactions:
        if famille == 'rotors':
            values = np.array([v for _, v in reactions])
            values[:, 2] += 9.81
            angles = np.arctan2(np.array(state[2])[:, 3], np.array(state[2])[:, 0])
            balance_error = max(float(np.max(np.abs(values))), float(np.max(np.abs(10.*(angles-.1)))))
        else:
            # Bilan des forces et moments au CdM, indépendant du G du noyau.
            values = np.array([v for _, v in reactions])
            efforts = np.zeros((nb, 6))
            points = [0.] + [.2*(i-1)+.2 for i in range(1, nb)]
            for i in range(nb):
                force, moment = values[i, :3], values[i, 3:]
                efforts[i, :3] += force
                efforts[i, 3:] += moment + np.cross([points[i]-.2*i, 0., 0.], force)
                if i:
                    efforts[i-1, :3] -= force
                    efforts[i-1, 3:] -= moment + np.cross([.2, 0., 0.], force)
            if famille == 'boucle':
                efforts[-1] += values[-1]
                # La solution de norme minimale est orthogonale aux six
                # auto-contraintes : torseur de fermeture transporté aux joints.
                null = np.zeros((6*(nb+1), 6))
                for k in range(6):
                    w = np.eye(6)[k]
                    for i in range(nb):
                        null[6*i:6*i+3, k] = -w[:3]
                        null[6*i+3:6*i+6, k] = -w[3:] - np.cross([.2*(nb-1)-points[i], 0., 0.], w[:3])
                    null[-6:, k] = w
                selfstress_error = float(np.linalg.norm(null.T@values.ravel()) /
                                         max(np.linalg.norm(null)*np.linalg.norm(values), 1e-300))
            efforts[:, 2] += 9.81
            balance_error = float(np.max(np.abs(efforts)))
    rotref = np.array([[np.cos(.1), -np.sin(.1), 0.],
                       [np.sin(.1), np.cos(.1), 0.], [0., 0., 1.]]) if famille == 'rotors' else np.eye(3)
    pose_error = max(float(np.max(np.abs(np.array(state[1])-[[.2*i, 0., 0.] for i in range(nb)]))),
                     float(np.max(np.abs(np.array(state[2]).reshape(nb, 3, 3)-rotref))))
    return dict(famille=famille, corps=nb, secondes=elapsed, rss_kib=rss,
                erreur=error, rapport=n.statique_info(), etat=state, reactions=reactions,
                erreur_reactions_absolue=reaction_error, erreur_pose_max=pose_error,
                erreur_equilibre_independant=balance_error, orthogonalite_autocontraintes=selfstress_error)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sortie', type=Path)
    p.add_argument('--tailles', type=int, nargs='+', default=[32, 64, 128, 256])
    p.add_argument('--familles', nargs='+', choices=('chaine', 'rotors', 'boucle'), default=['chaine', 'rotors'])
    p.add_argument('--worker', nargs=2, metavar=('FAMILLE', 'CORPS'))
    p.add_argument('--repetitions', type=int, default=3)
    args = p.parse_args()
    if args.worker:
        print(json.dumps(worker(args.worker[0], int(args.worker[1]))))
        return
    if args.sortie is None:
        p.error('--sortie requis')
    cpu = min(os.sched_getaffinity(0))
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    result = dict(cpu=cpu, extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
                  programme_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), mesures=[])
    for family in args.familles:
        for nb in args.tailles:
            for rep in range(args.repetitions):
                cmd = ['taskset', '-c', str(cpu), sys.executable, str(Path(__file__).resolve()), '--worker', family, str(nb)]
                try:
                    c = subprocess.run(cmd, env=env, text=True, capture_output=True, timeout=90)
                    m = json.loads(c.stdout) if c.returncode == 0 else dict(erreur=c.stderr)
                    m.update(code_sortie=c.returncode, stderr=c.stderr)
                except subprocess.TimeoutExpired:
                    m = dict(erreur='délai de 90 s dépassé')
                m.update(famille=family, corps=nb, repetition=rep)
                result['mesures'].append(m)
                args.sortie.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
                print(family, nb, rep, m.get('secondes'), m.get('rss_kib'), m['erreur'], flush=True)


if __name__ == '__main__':
    main()
