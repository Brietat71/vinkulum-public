"""Références analytiques et limite câble des formulations de poutre.

Les mesures d'erreur sont indépendantes de MBDyn. Le câble réutilise le
modèle de la campagne interne, avec un délai borné par processus.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from vinkulum import Noyau, _vinkulum

ROOT = Path(__file__).resolve().parents[1]
J = np.eye(3).ravel().tolist()


def beam(ne, formulation):
    n = Noyau([0., 0., 0.])
    for i in range(ne+1):
        n.corps(str(i), 1., J, [i/ne, 0., 0.])
    for i in range(ne):
        n.poutre(str(i), i, i+1, 1e6, 1e5, 50., 100., formulation=formulation)
    return n


def modele_factory(formulation):
    # Adaptateur de constructeur : paramètres et géométrie du corpus
    # inchangés ; seule l'option explicite de formulation diffère.
    class Modele:
        def __init__(self, *args, **kwargs):
            self.n = Noyau(*args, **kwargs)

        def __getattr__(self, key):
            return getattr(self.n, key)

        def poutre(self, *args, **kwargs):
            return self.n.poutre(*args, **kwargs, formulation=formulation)
    return Modele


def cable_worker(formulation):
    from vinkulum import campagne
    campagne.Noyau = modele_factory(formulation)
    n, ids = campagne._chaine(20, True)
    n.statique()
    print(json.dumps({'fleche_m': campagne._fleche(n, ids)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--formulation', choices=('milieu', 'integree'), required=True)
    parser.add_argument('--sortie', type=Path)
    parser.add_argument('--cable-worker', action='store_true')
    parser.add_argument('--adjoint-worker', action='store_true')
    args = parser.parse_args()
    if args.adjoint_worker:
        import vinkulum
        from vinkulum.adjoint_temps import demo_pont
        original = vinkulum.Noyau
        vinkulum.Noyau = modele_factory(args.formulation)
        try:
            demo_pont(rapide=True)
        finally:
            vinkulum.Noyau = original
        return
    if args.cable_worker:
        cable_worker(args.formulation)
        return
    if args.sortie is None:
        parser.error('--sortie requis')
    data = {'formulation': args.formulation,
            'extension_sha256': hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
            'sources_sha256': {str(p): hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                              for p in map(Path, ('src/lib.rs', 'src/tangent.rs',
                                                 'ci/diagnostic_poutres.py', 'python/vinkulum/campagne.py'))},
            'parameters': {'L': 1., 'EA': 1e6, 'GA': 1e5, 'GJ': 50., 'EI': 100.},
            'lineaire': {}, 'moment': {}}
    for ne in (1, 2, 4, 8):
        n = beam(ne, args.formulation)
        k = np.asarray(n.k_c_m_z()[0])[6:, 6:]
        force = np.zeros(6*ne)
        force[-5] = 1.
        tip = float(np.linalg.solve(k, force)[-5])
        exact = 1/(3*100.)+1/1e5
        data['lineaire'][ne] = {'fleche_m': tip, 'reference_m': exact,
                               'erreur_relative': abs(tip/exact-1)}
        n.liaison('enc', None, 0, bloque_t=[0, 1, 2], bloque_r=[0, 1, 2])
        n.effort(ne, [0., 0., 0.], [0., 0., 100.])
        n.statique(tol=1e-9, iters=100)
        exact_tip = [math.sin(1.), 1-math.cos(1.), 0.]
        tip = n.pose(ne)[0]
        data['moment'][ne] = {'position': tip, 'reference': exact_tip,
                             'erreur_m': float(np.max(np.abs(np.array(tip)-exact_tip))),
                             'residu_libre': float(np.max(np.abs(n.residu_statique()[6:])))}
    env = dict(os.environ, RAYON_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    try:
        p = subprocess.run([sys.executable, str(Path(__file__).resolve()),
                            '--formulation', args.formulation, '--cable-worker'],
                           capture_output=True, text=True, env=env, timeout=20)
        data['cable_20'] = {'statut': 'ok' if p.returncode == 0 else 'echec',
                            'code': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}
    except subprocess.TimeoutExpired:
        data['cable_20'] = {'statut': 'delai_depasse', 'limite_secondes': 20}
    args.sortie.write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n')
    print(args.sortie)


if __name__ == '__main__':
    main()
