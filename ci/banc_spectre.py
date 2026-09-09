"""Mesure isolée des modes complexes historiques sur une console de 30 corps.

Comparer deux paquets dans des processus neufs, en ordres AB et BA, sans
compilation concurrente. --paquet permet d'utiliser une extension de référence.
Le temps inclut K, C, M, base admissible et Schur ; la construction est exclue.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--paquet', type=Path)
    parser.add_argument('--repetitions', type=int, default=7)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error('repetitions >= 1')
    if args.paquet:
        sys.path.insert(0, str(args.paquet.resolve()))
    import numpy as np
    from vinkulum import Noyau, _vinkulum

    n = Noyau(g=[0, 0, 0])
    for i in range(30):
        n.corps(str(i), 1, np.eye(3).ravel().tolist(), [float(i), 0, 0])
        if i:
            n.poutre(str(i), i-1, i, 1000, 100, 10, 10)
    n.liaison('fixe', None, 0)
    mesures = []
    for i in range(args.repetitions + 1):
        debut = time.perf_counter()
        modes = n.modes_complexes(combien=12)
        duree = time.perf_counter() - debut
        if i:
            mesures.append(dict(secondes=duree, modes=modes))
    temps = [r['secondes'] for r in mesures]
    mediane = statistics.median(temps)
    print(json.dumps(dict(
        extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
        python=platform.python_version(), rayon=os.environ.get('RAYON_NUM_THREADS'),
        mediane_s=mediane,
        mad_s=statistics.median(abs(t-mediane) for t in temps), mesures=mesures,
    )))


if __name__ == '__main__':
    main()
