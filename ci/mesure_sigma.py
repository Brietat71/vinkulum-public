"""Expérience sigma : mêmes équations, binaire, sorties et tolérance Newton.

Références indépendantes Euler + Rdot=R[omega_b] (DOP853) ou rotations
commandées exactes. La branche flexible utilise une référence interne affinée,
explicitement identifiée. Aucun code de solveur concurrent.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

import vinkulum
from vinkulum import Noyau, _vinkulum


def skew(w):
    x, y, z = w
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def rot(w):
    return Rotation.from_rotvec(w).as_matrix()


def params(name):
    if name == 'affine':
        return dict(j=np.eye(3), r0=rot([.1, .3, -.2]), c=np.zeros(3),
                    w0=np.array([1.1, -.4, .7]), moment=np.array([-.3, .8, .2]),
                    g=np.zeros(3), m=1., pivot=False)
    if name == 'libre':
        return dict(j=np.diag([.7, 1.1, 1.6]), r0=rot([.1, .3, -.2]), c=np.zeros(3),
                    w0=np.array([2., -1., 3.]), moment=np.zeros(3),
                    g=np.zeros(3), m=1., pivot=False)
    if name == 'plan':
        return dict(j=1e-3*np.eye(3), r0=np.eye(3), c=np.array([.5, 0., -np.sqrt(.75)]),
                    w0=np.zeros(3), moment=np.zeros(3), g=np.array([0., 0., -9.80665]),
                    m=.2, pivot=True)
    if name == 'toupie':
        r0 = rot([np.pi/6, 0., 0.])
        return dict(j=np.diag([3e-4, 3e-4, 1e-4]), r0=r0, c=np.array([0., 0., .05]),
                    w0=r0 @ [0., 0., 300.], moment=np.zeros(3),
                    g=np.array([0., 0., -9.80665]), m=.1, pivot=True)
    raise ValueError(name)


def model(name):
    if name == 'boucle_spatiale':
        from modeles_confrontation import spatial
        return spatial()[0]
    if name == 'axes_commandes':
        n = Noyau([0., 0., 0.])
        n.corps('porteur', 1., np.eye(3).ravel().tolist(), [0., 0., 0.], w=[2., 0., 0.])
        n.corps('arbre', 1., np.eye(3).ravel().tolist(), [0., 0., 0.], w=[2., 0., 3.])
        n.liaison('x', None, 0, cible_r=([1., 0., 0.], ('lineaire', [0., 2.])))
        n.liaison('z', 0, 1, cible_r=([0., 0., 1.], ('lineaire', [0., 3.])))
        return n
    if name == 'branche_flexible':
        n = Noyau([0., 0., 0.])
        points = [[0., 0., 0.], [.5, 0., 0.], [1., .3, 0.], [1., -.3, .2]]
        for i, p in enumerate(points):
            n.corps(str(i), .2, np.diag([.002, .003, .004]).ravel().tolist(), p,
                    v=[0., .1, -.1] if i > 1 else [0., 0., 0.],
                    w=[.4, -.2, .3] if i > 1 else [0., 0., 0.])
        n.liaison('base', None, 0)
        for i, j in [(0, 1), (1, 2), (1, 3)]:
            n.poutre(f'{i}_{j}', i, j, 200., 100., .3, .5)
        return n
    p = params(name)
    n = Noyau(p['g'].tolist())
    n.corps(name, p['m'], p['j'].ravel().tolist(), (p['r0'] @ p['c']).tolist(),
            rot=p['r0'].ravel().tolist(), w=p['w0'].tolist(),
            v=np.cross(p['w0'], p['r0'] @ p['c']).tolist())
    if p['pivot']:
        n.liaison('pointe', None, 0, bloque_r=[])
    n.effort(0, [0., 0., 0.], p['moment'].tolist())
    return n


def reference(name, end, rtol=2e-12):
    if name == 'axes_commandes':
        def exact(t):
            t = np.asarray(t)
            r = np.zeros((len(t), 2, 3))
            mats = np.array([[rot([2*s, 0., 0.]), rot([2*s, 0., 0.]) @ rot([0., 0., 3*s])]
                             for s in t])
            w = np.array([[[2., 0., 0.], [2., -3*np.sin(2*s), 3*np.cos(2*s)]] for s in t])
            return r, mats, r.copy(), w
        return exact, {'type': 'rotations_commandées_exactes'}
    p = params(name)
    j, c, mass = p['j'], p['c'], p['m']
    jp = j + mass*((c @ c)*np.eye(3)-np.outer(c, c)) if p['pivot'] else j

    def rhs(t, y):
        r = y[:9].reshape(3, 3)
        wb = y[9:]
        torque = r.T @ p['moment']
        if p['pivot']:
            torque += np.cross(c, r.T @ (mass*p['g']))
        return np.r_[(r @ skew(wb)).ravel(), np.linalg.solve(jp, torque-np.cross(wb, jp @ wb))]

    y0 = np.r_[p['r0'].ravel(), p['r0'].T @ p['w0']]
    sol = solve_ivp(rhs, (0., end), y0, method='DOP853', rtol=rtol, atol=rtol/100,
                    dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)

    def evaluate(t):
        y = sol.sol(t).T
        mats = y[:, :9].reshape(-1, 3, 3)
        wb = y[:, 9:]
        # Pas de projection de R dans l'EDO : l'écart d'orthogonalité est contrôlé.
        r = np.einsum('tij,j->ti', mats, c)
        w = np.einsum('tij,tj->ti', mats, wb)
        v = np.cross(w, r)
        return r[:, None], mats[:, None], v[:, None], w[:, None]

    if p['pivot']:
        def reaction_norm(t):
            ys = sol.sol(t).T
            values = []
            for ti, y in zip(t, ys):
                r, wb = y[:9].reshape(3, 3), y[9:]
                ab = rhs(ti, y)[9:]
                acceleration = r @ (np.cross(ab, c)+np.cross(wb, np.cross(wb, c)))
                values.append(mass*np.linalg.norm(p['g']-acceleration))
            return np.array(values)
        evaluate.reaction_norm = reaction_norm

    mats = sol.sol(np.linspace(0., end, 501))[:9].T.reshape(-1, 3, 3)
    orth = np.max(np.abs(mats.transpose(0, 2, 1) @ mats-np.eye(3)))
    return evaluate, {'type': 'Euler_repère_corps_DOP853', 'rtol': rtol,
                      'atol': rtol/100, 'nfev': sol.nfev, 'orthogonalite_max': float(orth)}


def unpack(tr):
    return tuple(np.asarray([row[k] for row in tr]).reshape(len(tr), -1, 3, 3)
                 if k == 2 else np.asarray([row[k] for row in tr]) for k in range(1, 5))


def errors(actual, expected):
    r, a, v, w = actual
    rr, b, vv, ww = expected
    relative = a @ b.transpose(0, 1, 3, 2)
    angular = Rotation.from_matrix(relative.reshape(-1, 3, 3)).magnitude()
    return {'position_m': float(np.max(np.linalg.norm(r-rr, axis=-1))),
            'orientation_rad': float(np.max(angular)),
            'vitesse_m_s': float(np.max(np.linalg.norm(v-vv, axis=-1))),
            'rotation_rad_s': float(np.max(np.linalg.norm(w-ww, axis=-1)))}


def run(name, end, h, sigma, rho=.9, output_step=None):
    n = model(name)
    every = max(1, round((output_step or h)/h))
    tic = time.perf_counter()
    tr = n.simule(end, h, rho=rho, tol=1e-12, newton_max=25, tous=every, sigma_lie=sigma)
    elapsed = time.perf_counter()-tic
    return tr, {'temps_s': elapsed, 'stats': list(n.stats()), 'profil_s': list(n.chronos()),
                'adapt': list(n.adapt_stats()), 'phi_max_final': float(np.max(np.abs(n.phi()), initial=0.)),
                'phi_dot_max_final': float(np.max(np.abs(n.phi_dot()), initial=0.))}


CASES = {'affine': (1., .04), 'libre': (2., .04), 'plan': (1., .02),
         'toupie': (.2, .001), 'axes_commandes': (1., .02), 'branche_flexible': (.2, .002),
         'boucle_spatiale': (1., .004)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sortie', type=Path, required=True)
    parser.add_argument('--cas', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--repetitions', type=int, default=3)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error('repetitions >= 1')
    cpu = min(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {cpu})
    rho = .9
    gamma = .5 + (1-rho)/(1+rho)
    beta = .25*(gamma+.5)**2
    sigmas = [0., gamma/(3*beta), 1.]
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    out = dict(version=vinkulum.__version__, binaire_sha256=sha(_vinkulum.__file__),
               script_sha256=sha(__file__), cpu=cpu, machine=platform.platform(),
               python=platform.python_version(), rho=rho, sigma=sigmas,
               repetitions=args.repetitions, tol_newton=1e-12,
               threads={k: os.environ.get(k) for k in ('RAYON_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS')},
               references={}, mesures=[])
    args.sortie.parent.mkdir(parents=True, exist_ok=True)

    def save():
        args.sortie.write_text(json.dumps(out, ensure_ascii=False, indent=2, allow_nan=False)+'\n')

    for name in args.cas:
        end, coarse = CASES[name]
        if name in ('branche_flexible', 'boucle_spatiale'):
            # Référence interne : convergence en h ET concordance entre sigma.
            grids = []
            for div, sigma in [(64, 0.), (128, 0.), (128, sigmas[1])]:
                tr, _ = run(name, end, coarse/div, sigma, output_step=coarse/8)
                grids.append(unpack(tr))
            out['references'][name] = {'type': 'Vinkulum_affiné_non_indépendant',
                'h': coarse/128, 'ecart_h': errors(grids[0], grids[1]),
                'ecart_sigma': errors(grids[1], grids[2])}
            # Toutes les sorties mesurées appartiennent à cette grille.
            times = np.asarray([row[0] for row in tr])
            def ref(t):
                indices = np.rint(np.asarray(t)/(coarse/8)).astype(int)-1
                if not np.allclose(times[indices], t, rtol=0., atol=2e-12):
                    raise RuntimeError('grille de référence incompatible')
                return tuple(a[indices] for a in grids[1])
        else:
            ref, meta = reference(name, end)
            tighter, _ = reference(name, end, rtol=2e-13)
            grid = np.linspace(0., end, 501)
            meta['ecart_reference_resserree'] = errors(ref(grid), tighter(grid))
            out['references'][name] = meta
        # Une chauffe complète pour chaque variante, hors chrono conservé.
        for sigma in sigmas:
            run(name, end, coarse, sigma)
        for rep in range(args.repetitions):
            for div in [1, 2, 4, 8]:
                for sigma in sigmas[::1 if rep % 2 == 0 else -1]:
                    h = coarse/div
                    record = {'cas': name, 'repetition': rep, 'h': h, 'sigma': sigma}
                    try:
                        tr, stats = run(name, end, h, sigma, output_step=coarse)
                        record.update(stats)
                        record['erreurs_grille_commune'] = errors(unpack(tr), ref(np.array([r[0] for r in tr])))
                        # Rejeu de validation hors chrono conservé : contrôler
                        # chaque pas, y compris le transitoire initial, sans
                        # faire payer plus de sorties aux candidats à petit h.
                        dense, _ = run(name, end, h, sigma)
                        if tr[-1] != dense[-1]:
                            raise RuntimeError('la cadence de sortie a modifié le résultat')
                        times_trial = np.array([r[0] for r in dense])
                        record['erreurs'] = errors(unpack(dense), ref(times_trial))
                        if hasattr(ref, 'reaction_norm'):
                            actual = np.linalg.norm(np.array([r[5] for r in dense]), axis=1)
                            record['erreurs']['norme_reaction_N'] = float(np.max(np.abs(actual-ref.reaction_norm(times_trial))))
                        record['nombre_sorties'] = len(tr)
                        record['nombre_etats_controles'] = len(dense)
                    except (RuntimeError, ValueError) as e:
                        record['echec'] = str(e)
                    out['mesures'].append(record)
                    save()
                    print(name, rep, h, sigma, record.get('erreurs', record.get('echec')), flush=True)
    save()


if __name__ == '__main__':
    main()
