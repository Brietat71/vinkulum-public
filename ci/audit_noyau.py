"""Banc isolé pour l'audit : JSON sur stdout, sans écriture implicite.

Exemple : RAYON_NUM_THREADS=1 .venv/bin/python ci/audit_noyau.py chaine900
--paquet /tmp/reference sélectionne une copie du paquet avant reconstruction.
Comparer des processus neufs, dans les deux ordres ; ne pas compiler en même temps.
Les compteurs de chronos sont des mesures internes, le RSS un pic du processus.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import sys
import sysconfig
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cas', choices=['pendule', 'chaine17', 'chaine18', 'chaine30',
                                 'chaine221', 'chaine223', 'chaine900', 'poutre8',
                                 'super8', 'contact', 'modal30', 'sensibilite30'])
    p.add_argument('--paquet', type=Path)
    p.add_argument('--repetitions', type=int, default=7)
    args = p.parse_args()
    if args.repetitions < 1:
        p.error('repetitions ≥ 1')
    if args.paquet:
        sys.path.insert(0, str(args.paquet.resolve()))
    import numpy as np
    from vinkulum import Noyau, bancs, _vinkulum

    def modele(nb):
        n = Noyau(g=[0, 0, 0])
        for i in range(nb):
            n.corps(str(i), 1, np.eye(3).ravel().tolist(), [i, 0, 0])
        return n

    def ponctuel():
        if args.cas == 'super8':
            n = modele(8)
            b = np.diff(np.eye(8), axis=0)
            k = np.kron(b.T @ b, np.eye(6)) * 100
            n.superelement('s', list(range(8)), k.ravel().tolist(), beta=.01)
            etat = list(n.etat())
            etat[1][-1][0] += .01
            n.pose_etat(*etat)
        elif args.cas == 'contact':
            n = modele(1)
            n.contact('sol', 0, [0, 0, 0], .01, k=1e4, c=.1, mu=.3)
        else:
            n = modele(30)
            for i in range(29):
                if args.cas == 'modal30':
                    n.liaison(str(i), i, i+1, bloque_r=[])
                else:
                    n.poutre(str(i), i, i+1, 1000, 100, 10, 10)
        debut = time.perf_counter()
        if args.cas == 'modal30':
            _, _, z = n.k_m_z()
            physique = dict(dimension=len(z[0]), norme_base=float(np.linalg.norm(z)))
        elif args.cas == 'sensibilite30':
            dk = n.d_raideur_poutre(10, 'ea')
            physique = dict(norme=float(np.linalg.norm(dk)))
        else:
            n.simule(.05, .001, tous=10**9)
            physique = dict(temps_final=n.t(), energie=n.energie(), etat=n.etat(),
                            phi=float(np.linalg.norm(n.phi())),
                            stats=n.stats(), adapt=n.adapt_stats(), chronos=n.chronos())
        return dict(temps=time.perf_counter()-debut, **physique)

    if args.cas.startswith('chaine'):
        nb = int(args.cas[6:])
        def f():
            n = Noyau()
            longueur, masse = .05, .02
            inertie = np.diag([1e-8, masse*longueur**2/12, masse*longueur**2/12])
            for i in range(nb):
                n.corps(str(i), masse, inertie.ravel().tolist(),
                        [longueur*(i+.5), 0, 0])
                if i == 0:
                    n.liaison('r0', None, 0, bloque_r=[])
                else:
                    n.liaison(str(i), i-1, i, pa=[longueur/2, 0, 0], bloque_r=[])
            debut = time.perf_counter()
            n.simule(.02, .001, tous=10**9)
            duree = time.perf_counter()-debut
            return dict(temps=duree, ms_pas=duree*1000/20, chronos=n.chronos(),
                        stats=n.stats(), adapt=n.adapt_stats(),
                        phi=float(np.linalg.norm(n.phi())), bout=n.pose(nb-1),
                        energie=n.energie(), temps_final=n.t())
    elif args.cas == 'pendule':
        f = lambda: bancs.pendule(.001)
    elif args.cas == 'poutre8':
        f = lambda: bancs.console(8)
    else:
        f = ponctuel
    f()
    sorties = [f() for _ in range(args.repetitions)]
    secondes = [r['temps'] for r in sorties]
    med = statistics.median(secondes)
    print(json.dumps(dict(cas=args.cas, python=sys.version, machine=platform.machine(),
                          gil_desactive=bool(sysconfig.get_config_var('Py_GIL_DISABLED')),
                          threads=os.environ.get('RAYON_NUM_THREADS'),
                          extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
                          mediane_s=med, mad_s=statistics.median(abs(t-med) for t in secondes),
                          min_s=min(secondes), max_s=max(secondes),
                          rss_pic_kio=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          sorties=sorties), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
