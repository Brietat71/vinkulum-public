"""Comparaison isolée ; PYTHONPATH sélectionne la roue, le modèle est commun.

RAYON_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python ci/mesure_multirythme.py
Un échauffement puis cinq répétitions, aucun autre banc en parallèle.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import statistics
import time

import numpy as np
import vinkulum
from vinkulum import _vinkulum

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / "python/vinkulum/_banc_multirythme.py"
spec = importlib.util.spec_from_file_location("modele_mesure", MODEL)
modele = importlib.util.module_from_spec(spec)
spec.loader.exec_module(modele)


def mesure(n, mode, T):
    noyau, a, c, b = modele.monte(n)
    start = time.perf_counter()
    if mode == "multi":
        tr = noyau.simule_multirythme(T, .001, [c, b], 100)
    else:
        tr = noyau.simule(T, 1e-5 if mode == "fin" else .001,
                          tous=100 if mode == "fin" else 1)
    duree = time.perf_counter() - start
    fin = [f for f in tr if f[0] > T - .05]
    t = np.array([f[0] for f in fin])
    d = np.array([np.linalg.norm(np.array(f[1][b]) - f[1][c]) - .05 for f in fin])
    d -= np.polyval(np.polyfit(t - t[0], d, 2), t - t[0])
    return dict(temps=duree, amplitude=float((d.max() - d.min()) / 2),
                position=tr[-1][1][a], etat=tr[-1],
                contraintes=noyau.phi(), stats=noyau.stats(),
                chronos=noyau.chronos())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tailles", nargs="+", type=int, default=[8, 16, 32])
    p.add_argument("--modes", nargs="+", choices=["fin", "gros", "multi"], default=["fin", "gros", "multi"])
    p.add_argument("--repetitions", type=int, default=5)
    p.add_argument("--duree", type=float, default=.3)
    p.add_argument("--sortie", type=Path, required=True)
    args = p.parse_args()
    rapport = dict(version=vinkulum.__version__, python=platform.python_version(),
                   plateforme=platform.platform(),
                   extension_sha256=hashlib.sha256(Path(_vinkulum.__file__).read_bytes()).hexdigest(),
                   modele_sha256=hashlib.sha256(MODEL.read_bytes()).hexdigest(),
                   affinite=sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
                   threads={k: os.environ.get(k) for k in ["RAYON_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"]},
                   duree=args.duree, h=.001, k=100, echauffements=1, mesures={})
    for n in args.tailles:
        for mode in args.modes:
            mesure(n, mode, args.duree)
            sorties = [mesure(n, mode, args.duree) for _ in range(args.repetitions)]
            temps = [s["temps"] for s in sorties]
            rapport["mesures"][f"{n}/{mode}"] = dict(mediane=statistics.median(temps),
                minimum=min(temps), maximum=max(temps), sorties=sorties)
            args.sortie.parent.mkdir(parents=True, exist_ok=True)
            args.sortie.write_text(json.dumps(rapport, indent=2) + "\n")
            print(f"{n}/{mode}: {statistics.median(temps):.3f} s", flush=True)


if __name__ == "__main__":
    main()
