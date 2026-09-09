"""Contre-calcul direct aux points témoins des configurations HCB écartées.

Diagnostic après campagne, sans classement de vitesse. Vérifie la même
paire réduite par une résolution directe distincte du préconditionneur.
Une réussite aux points choisis ne qualifie pas les 257 fréquences.
"""
import argparse
import json
import os
from pathlib import Path

from confronte_ports_exudyn import (THREADS, analyse, distribution_figee,
                                    ecrire, juger, lire, sha)


def controler(dossier, sortie, cpu):
    import sys
    import numpy as np
    import exudyn
    from reference_hcb_exudyn import prepare, reponse
    dossier = Path(dossier)
    rapport = json.loads((dossier/"rapport.json").read_text())
    bilan = analyse(rapport)
    ecartes = {(c["n"], c["max_hz"], c["variante"]) for c in bilan["configurations"]
               if c["variante"].startswith("hcb") and not c["eligible"]}
    module = sys.modules[exudyn.SystemContainer.__module__]
    resultats = dict(schema=1, campagne_sha256=sha(dossier/"rapport.json"),
        source_sha256=sha(__file__), exudyn=exudyn.__version__,
        extension_sha256=sha(module.__file__), distribution=distribution_figee("exudyn"),
        cpu=sorted(os.sched_getaffinity(0)), fils={k: os.environ.get(k) for k in THREADS},
        portee="Diagnostic sur les pires points mesurés, zéro classement de vitesse ; "
               "une réussite ponctuelle ne qualifie pas toute la grille",
        sources_sha256={nom: sha(Path(__file__).parent/nom)
                        for nom in ("reference_hcb_exudyn.py", "confronte_ports_exudyn.py")},
        configurations=[])
    if resultats["cpu"] != [cpu] or resultats["fils"] != THREADS:
        raise ValueError("budget de calcul incohérent")
    for cas in rapport["cas"]:
        nom = f"n{cas['n']}-f{cas['max_hz']:g}"
        for config in cas["configurations"]:
            if (cas["n"], cas["max_hz"], config["variante"]) not in ecartes:
                continue
            ids = {0, cas["nombre_frequences"]-1}
            for e in config["essais"]:
                if e["statut"] == "termine":
                    for courbe in e["juge"]["operateurs"].values():
                        ids.add(int(np.argmax(courbe)))
                    for courbes in e["juge"]["erreurs"].values():
                        ids.add(int(np.argmax(np.max(courbes, axis=1))))
            ids = sorted(ids)
            entree, oracle = dossier/(nom+".npz"), dossier/(nom+".ref.npy")
            if sha(entree) != cas["entree_sha256"] or sha(oracle) != cas["reference"]["fichier_sha256"]:
                raise ValueError("entrée ou oracle différent de la campagne")
            d, m, metrique, forces, omega = lire(entree)
            ref = np.load(oracle, mmap_mode="r", allow_pickle=False)[ids]
            np.random.seed(42123)
            p = prepare(d, m, config["modes"], projection=config["variante"].split("_")[1],
                        diagnostics_complets=False)
            lignes = []
            for methode in ("spectrale_corrigee", "directe"):
                x = np.asarray([reponse(p, omega[j], forces, resolution=methode) for j in ids])
                lignes.append(dict(methode=methode, juge=juger(d, m, metrique, x, ref)))
            resultats["configurations"].append(dict(n=cas["n"], max_hz=cas["max_hz"],
                variante=config["variante"], modes=config["modes"], indices=ids,
                frequences_hz=(omega[ids]/(2*np.pi)).tolist(), methodes=lignes))
            ecrire(sortie, resultats)
            print(nom, config["variante"], {r["methode"]: r["juge"]["maxima_operateurs"]
                                           for r in lignes}, flush=True)
    ecrire(sortie, resultats)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("dossier", type=Path)
    p.add_argument("sortie", type=Path)
    p.add_argument("--cpu", type=int, default=8)
    a = p.parse_args()
    os.sched_setaffinity(0, {a.cpu})
    controler(a.dossier, a.sortie, a.cpu)
