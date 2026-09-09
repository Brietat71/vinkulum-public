"""Mesures sans écriture de rapport : un échauffement, sept répétitions.

RAYON_NUM_THREADS=1 .venv/bin/python ci/mesure_noyau.py
Pour comparer deux extensions, choisir leur paquet avec PYTHONPATH.
Les pics RSS sont ceux du processus, cumulés ; comparer des processus neufs.
"""
import json
import os
import resource
import statistics

from vinkulum import bancs


def main():
    mesures = {}
    cas = {
        "pendule": lambda: bancs.pendule(.001),
        "chaine30": lambda: bancs.chaine(30),
        "poutre8": lambda: bancs.console(8),
        "creux900": lambda: bancs.chaine(900, t_end=.005),
    }
    for nom, f in cas.items():
        f()
        sorties = [f() for _ in range(7)]
        temps = [r["ms_pas"] if "ms_pas" in r else r["temps"] * 1e3 / r["pas"] for r in sorties]
        mesures[nom] = dict(ms_pas_mediane=statistics.median(temps),
                            ms_pas_min=min(temps), ms_pas_max=max(temps),
                            sorties=sorties, rss_pic_kio=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print(json.dumps(dict(threads=os.environ.get("RAYON_NUM_THREADS"), mesures=mesures), indent=2))


if __name__ == "__main__":
    main()
