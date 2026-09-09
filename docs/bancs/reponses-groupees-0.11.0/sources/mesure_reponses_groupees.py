"""Mesure des réponses contrôlées, roues 0.10.0/0.11.0 et témoins LU/HCB.

Entrées et références de la confrontation des ports, sans nouveau réglage
HCB. Champs complets, mêmes six charges et juge commun. Sous Linux, VmHWM
mesure le pic de l'image exécutée, avant chargement de l'oracle.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time

from confronte_ports_exudyn import (CIBLE, THREADS, distribution_figee, ecrire,
                                    juger, lire, worker as temoin)

ROOT = Path(__file__).resolve().parents[1]
VARIANTES = ("reference_010", "api_011", "lot_011", "lu_corrigee", "hcb_standard", "hcb_energie")
SOURCES = ("mesure_reponses_groupees.py", "confronte_ports_exudyn.py",
           "reference_hcb_exudyn.py", "oracle_champs_ports.py",
           "reference_ports_precision.py", "experience_reduction_native.py")


def sha(path):
    # Ne pas charger le grand oracle pour calculer son empreinte avant le
    # chronomètre et la mesure de VmHWM.
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for bloc in iter(lambda: f.read(1024*1024), b""):
            h.update(bloc)
    return h.hexdigest()


def worker(entree, reference, variante, modes, sortie):
    import numpy as np
    import scipy
    if variante not in VARIANTES:
        raise ValueError("variante inconnue")
    if variante in ("lu_corrigee", "hcb_standard", "hcb_energie"):
        temoin(entree, reference, variante, modes, sortie)
        r = json.loads(Path(sortie).read_text())
        r["memoire_comparable"] = False  # ru_maxrss de l'ancien pilote, laissé brut
        ecrire(sortie, r)
        return
    import vinkulum
    from vinkulum.reduction_ports import ReductionMaterielle
    from vinkulum._ports import certificat_spectral, inverse_selectionnee
    attendu = "0.10.0" if variante == "reference_010" else "0.11.0"
    if vinkulum.__version__ != attendu:
        raise ValueError("roue de mauvaise version")
    d, m, metrique, forces, omegas = lire(entree)
    paquet = Path(vinkulum.__file__).parent
    r = dict(variante=variante, modes_demandes=modes, moteur=vinkulum.__version__,
        extension_sha256=sha(next(paquet.glob("*.so"))), distribution=distribution_figee("vinkulum"),
        sources_python_sha256={str(p.relative_to(paquet)): sha(p) for p in sorted(paquet.rglob("*.py"))},
        entree_sha256=sha(entree), reference_sha256=sha(reference), python=sys.version,
        numpy=np.__version__, scipy=scipy.__version__, cpu=sorted(os.sched_getaffinity(0)),
        fils={k: os.environ.get(k) for k in THREADS}, nombre_frequences=len(omegas), charges=6,
        date_utc=datetime.now(timezone.utc).isoformat(), statut="en_cours", cible=CIBLE,
        memoire_comparable=True)
    n = d.shape[1]
    champs = np.empty((len(omegas), n, 6))
    bornes = {k: [] for k in ("masse", "deformation")}
    absolues = {k: [] for k in bornes}
    normes = {k: [] for k in bornes}
    debut, phase = time.perf_counter(), "preparation"
    try:
        reduction = ReductionMaterielle(d, m, np.arange(n-6), np.arange(n-6, n),
            metrique, omegas[-1], tolerance=1e-10, max_blocs=12, max_directions=128)
        r.update(directions=reduction.taille_interieure,
                 certificat_spectral=reduction.certificat_spectral,
                 borne_schur=reduction.borne_uniforme, statut_schur=reduction.statut)
        prepare = time.perf_counter()
        phase = "reponses"
        for j, w in enumerate(omegas):
            reps = (reduction.reponses(w, forces) if variante == "lot_011" else
                    [reduction.reponse(w, forces[:, col]) for col in range(6)])
            for col, rep in enumerate(reps):
                champs[j, :, col] = rep["champ_physique"]
            for nom in bornes:
                bornes[nom].append([rep[nom]["borne_relative"] for rep in reps])
                absolues[nom].append([rep[nom]["borne_absolue"] for rep in reps])
                normes[nom].append([rep[nom]["norme_candidate"] for rep in reps])
        termine = time.perf_counter()
        memoire = {ligne.split(":")[0]: int(ligne.split()[1])
                   for ligne in Path("/proc/self/status").read_text().splitlines()
                   if ligne.startswith(("VmHWM:", "VmRSS:"))}
        r.update(statut="termine", preparation_s=prepare-debut,
                 reponses_s=termine-prepare, total_s=termine-debut,
                 pic_image_kib=memoire["VmHWM"], rss_image_kib=memoire["VmRSS"])
        phase = "juge"
        ref = np.load(reference, mmap_mode="r", allow_pickle=False)
        r["juge"] = juger(d, m, metrique, champs, ref)
        r.update(bornes_api=bornes, bornes_absolues=absolues, normes_candidates=normes)
        # Diagnostic des majorants flottants par colonne, séparé du juge.
        controles = {nom: [] for nom in bornes}
        for j, (x, y) in enumerate(zip(champs, ref)):
            delta = x-y
            for nom, erreur, norme in (
                ("masse", np.sqrt(np.maximum(np.sum(delta*(m@delta), axis=0), 0.)),
                 np.sqrt(np.sum(y*(m@y), axis=0))),
                ("deformation", np.linalg.norm(d@delta, axis=0), np.linalg.norm(d@y, axis=0))):
                controles[nom].append([b is not None and e <= b+64*np.finfo(float).eps*ny
                                      for b, e, ny in zip(absolues[nom][j], erreur, norme)])
        r["controle_bornes"] = {n: np.asarray(v, dtype=bool).tolist() for n, v in controles.items()}
        r["marge_arrondi_controle_bornes"] = "64 epsilon * norme de référence ; pas un certificat machine"
    except Exception as exc:
        r.update(statut="refus_ou_echec", phase=phase,
                 erreur=dict(type=type(exc).__name__, message=str(exc)),
                 temps_jusquau_refus_s=time.perf_counter()-debut)
    ecrire(sortie, r)


def analyse(rapport):
    attendus = {(n, float(hz)) for n in (32, 128, 512) for hz in (2, 20, 40)}
    cas = rapport["cas"]
    if len(cas) != 9 or {(c["n"], c["max_hz"]) for c in cas} != attendus:
        raise ValueError("familles de cas incomplètes")
    identites, environnements, lignes = {}, set(), []
    for c in cas:
        ref = c["reference"]
        bornes_ref = [ref[cle][n] for cle in ("maxima", "maxima_operateurs")
                      for n in ("masse", "deformation", "port")]
        if (ref["precision"] != [70, 90] or max(bornes_ref) != ref["ecart_relatif_max"]
                or not all(math.isfinite(v) and 0 <= v <= CIBLE/100 for v in bornes_ref)):
            raise ValueError("référence non qualifiée")
        if len(c["configurations"]) != 6 or {v["variante"] for v in c["configurations"]} != set(VARIANTES):
            raise ValueError("configurations incomplètes")
        for config in c["configurations"]:
            variante = config["variante"]
            es = config["essais"]
            roles = [(e["role"], e["repetition"]) for e in es]
            if len(set(roles)) != len(roles):
                raise ValueError("répétition dupliquée")
            if not set(roles) <= {("chauffe", -1), ("mesure", 0), ("mesure", 1), ("mesure", 2)}:
                raise ValueError("rôle de répétition inconnu")
            complet = set(roles) == {("chauffe", -1), ("mesure", 0), ("mesure", 1), ("mesure", 2)}
            eligible, bornes_ok = complet, complet
            for e in es:
                if e["cpu"] != [rapport["cpu"]] or e["fils"] != THREADS:
                    raise ValueError("budget de calcul différent")
                if e["variante"] != variante or e["modes_demandes"] != config["modes"]:
                    raise ValueError("configuration incohérente")
                if e["entree_sha256"] != c["entree_sha256"] or e["reference_sha256"] != c["reference"]["fichier_sha256"]:
                    raise ValueError("entrées ou référence différentes")
                environnements.add((e["python"], e["numpy"], e["scipy"]))
                if "moteur" in e:
                    cle = ("v010" if variante == "reference_010" else "v011"
                           if variante in ("api_011", "lot_011") else "hcb")
                    identite = json.dumps({k: e.get(k) for k in ("moteur", "extension_sha256",
                        "distribution", "sources_python_sha256")}, sort_keys=True)
                    if cle in identites and identites[cle] != identite:
                        raise ValueError("installation différente dans une famille")
                    identites[cle] = identite
                    if cle in ("v010", "v011") and e["moteur"] != ("0.10.0" if cle == "v010" else "0.11.0"):
                        raise ValueError("version incorrecte")
                if e["statut"] != "termine":
                    eligible = bornes_ok = False
                    continue
                if e["nombre_frequences"] != 257 or e["charges"] != 6:
                    raise ValueError("couverture incomplète")
                bon = True
                for nom in ("masse", "deformation", "port"):
                    v, op = e["juge"]["erreurs"][nom], e["juge"]["operateurs"][nom]
                    if len(v) != 257 or any(len(a) != 6 for a in v) or len(op) != 257:
                        raise ValueError("couverture du juge incomplète")
                    vals = [x for a in v for x in a]
                    if not all(math.isfinite(x) and x >= 0 for x in vals+op):
                        raise ValueError("erreur invalide")
                    if max(vals) != e["juge"]["maxima"][nom] or max(op) != e["juge"]["maxima_operateurs"][nom]:
                        raise ValueError("maximum incohérent")
                    bon &= max(vals+op) <= CIBLE
                if bool(bon) != e["juge"]["accepte"]:
                    raise ValueError("qualification incohérente")
                eligible &= bon
                for cle in ("preparation_s", "reponses_s", "total_s"):
                    if not math.isfinite(e[cle]) or e[cle] <= 0:
                        raise ValueError("temps invalide")
                if abs(e["total_s"]-e["preparation_s"]-e["reponses_s"]) > 1e-9:
                    raise ValueError("temps total incohérent")
                if variante in VARIANTES[:3]:
                    if (not e["memoire_comparable"] or e["pic_image_kib"] <= 0
                            or e["rss_image_kib"] <= 0 or e["pic_image_kib"] < e["rss_image_kib"]):
                        raise ValueError("compteur de mémoire invalide")
                    for nom in ("masse", "deformation"):
                        controles = e["controle_bornes"][nom]
                        if len(controles) != 257 or any(len(a) != 6 for a in controles):
                            raise ValueError("couverture des bornes incomplète")
                        if not all(isinstance(v, bool) for a in controles for v in a):
                            raise ValueError("contrôle des bornes invalide")
                        bornes_ok &= all(v is True for a in controles for v in a)
            ligne = dict(n=c["n"], max_hz=c["max_hz"], variante=variante,
                         modes=config["modes"], complet=complet, eligible=bool(eligible))
            mesures = [e for e in es if e["role"] == "mesure"]
            if complet and all(e["statut"] == "termine" for e in mesures):
                for cle in ("preparation_s", "reponses_s", "total_s"):
                    ligne[cle] = statistics.median(e[cle] for e in mesures)
                    ligne[cle+"_plage"] = [min(e[cle] for e in mesures), max(e[cle] for e in mesures)]
                ligne["erreurs_operateur_max"] = {n: max(e["juge"]["maxima_operateurs"][n] for e in mesures)
                                                 for n in ("masse", "deformation", "port")}
                if variante in VARIANTES[:3]:
                    ligne["bornes_confrontees"] = bool(bornes_ok)
                    ligne["pic_image_kib"] = statistics.median(e["pic_image_kib"] for e in mesures)
            lignes.append(ligne)
    if len(environnements) != 1:
        raise ValueError("versions numériques différentes")
    return dict(cible=CIBLE, configurations=lignes)


def campagne(dossier, entrees, py010, py011, pyhcb, cpu):
    dossier, entrees = Path(dossier), Path(entrees)
    dossier.mkdir(parents=True, exist_ok=True)
    precedent = ROOT/"docs/bancs/confrontation-ports-exudyn-0.10.0/bilan.json"
    rangs = {(v["n"], v["max_hz"], v["variante"]): v["modes"]
             for v in json.loads(precedent.read_text())["configurations"]}
    r = dict(schema=1, date_utc=datetime.now(timezone.utc).isoformat(), cpu=cpu, cas=[],
             sources_sha256={n: sha(ROOT/"ci"/n) for n in SOURCES},
             choix_hcb_sha256=sha(precedent), protocole=dict(
                 frequences=257, charges=6, echauffements=1, repetitions=3,
                 ordre_inverse_un_passage_sur_deux=True, processus_frais=True,
                 cout="préparation et champs avec normes ; les trois API Vinkulum ajoutent les bornes",
                 choix_hcb="rangs de la campagne corrigée 0.10.0, aucun réglage supplémentaire",
                 memoire="VmHWM avant oracle pour Vinkulum uniquement ; ru_maxrss des témoins exclu"))
    env = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE="1")
    for n in (32, 128, 512):
        for hz in (2., 20., 40.):
            nom = f"n{n}-f{hz:g}"
            c = json.loads((entrees/(nom+".npz.json")).read_text())
            c.update(reference=json.loads((entrees/(nom+".ref.npy.json")).read_text()), configurations=[])
            for variante in VARIANTES:
                c["configurations"].append(dict(variante=variante, modes=rangs.get((n, hz, variante), 0), essais=[]))
            r["cas"].append(c)
            for passage in range(4):
                ordre = c["configurations"] if passage % 2 == 0 else c["configurations"][::-1]
                for config in ordre:
                    variante, modes = config["variante"], config["modes"]
                    py = py010 if variante == "reference_010" else pyhcb if variante.startswith("hcb") else py011
                    role, rep = ("chauffe", -1) if passage == 0 else ("mesure", passage-1)
                    stem = f"{nom}-{variante}-{role}{rep}"
                    sortie = dossier/(stem+".json")
                    print(stem, flush=True)
                    with (dossier/(stem+".log")).open("w") as log:
                        commande = [str(py), str(Path(__file__).resolve()), "--worker",
                            str(entrees/(nom+".npz")), str(entrees/(nom+".ref.npy")),
                            variante, str(modes), str(sortie), str(cpu)]
                        subprocess.run(commande, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
                    essai = json.loads(sortie.read_text())
                    essai.update(role=role, repetition=rep, journal_sha256=sha(dossier/(stem+".log")))
                    config["essais"].append(essai)
    r["bilan"] = analyse(r)
    ecrire(dossier/"rapport.json", r)
    print("Campagne terminée : "+str(dossier/"rapport.json"), flush=True)


def archiver(dossier, entrees, destination):
    dossier, entrees, destination = map(Path, (dossier, entrees, destination))
    r = json.loads((dossier/"rapport.json").read_text())
    bilan = analyse(r)
    destination.mkdir(parents=True, exist_ok=True)
    (destination/"essais.json.gz").write_bytes(gzip.compress(
        json.dumps(r, ensure_ascii=False, allow_nan=False).encode(), mtime=0))
    ecrire(destination/"bilan.json", bilan)
    (destination/"sources").mkdir(exist_ok=True)
    for n, h in r["sources_sha256"].items():
        if sha(ROOT/"ci"/n) != h:
            raise ValueError("source modifiée depuis les mesures : "+n)
        shutil.copy2(ROOT/"ci"/n, destination/"sources"/n)
    (destination/"entrees").mkdir(exist_ok=True)
    for c in r["cas"]:
        n = f"n{c['n']}-f{c['max_hz']:g}.npz"
        if sha(entrees/n) != c["entree_sha256"]:
            raise ValueError("entrée modifiée")
        shutil.copy2(entrees/n, destination/"entrees"/n)
    logs = {p.name: p.read_text() for p in dossier.glob("*.log")}
    (destination/"journaux.json.gz").write_bytes(gzip.compress(json.dumps(logs).encode(), mtime=0))
    ecrire(destination/"manifest.json", dict(version="0.11.0",
        fichiers_sha256={str(p.relative_to(destination)): sha(p)
                         for p in sorted(destination.rglob("*")) if p.is_file() and p.name != "manifest.json"}))


def verifier(destination):
    destination = Path(destination)
    m = json.loads((destination/"manifest.json").read_text())["fichiers_sha256"]
    presents = {str(p.relative_to(destination)) for p in destination.rglob("*") if p.is_file() and p.name != "manifest.json"}
    indispensables = {"essais.json.gz", "bilan.json", "journaux.json.gz"}
    indispensables.update("sources/"+n for n in SOURCES)
    indispensables.update(f"entrees/n{n}-f{hz}.npz" for n in (32, 128, 512) for hz in (2, 20, 40))
    if set(m) != presents or not indispensables <= presents:
        raise ValueError("archive incomplète")
    if any(sha(destination/n) != h for n, h in m.items()):
        raise ValueError("archive altérée")
    r = json.loads(gzip.decompress((destination/"essais.json.gz").read_bytes()))
    if analyse(r) != json.loads((destination/"bilan.json").read_text()):
        raise ValueError("bilan incohérent")
    if set(r["sources_sha256"]) != set(SOURCES):
        raise ValueError("sources incomplètes")
    for n, h in r["sources_sha256"].items():
        if sha(destination/"sources"/n) != h:
            raise ValueError("provenance de source incohérente")
    logs = json.loads(gzip.decompress((destination/"journaux.json.gz").read_bytes()))
    import hashlib
    for c in r["cas"]:
        nom = f"n{c['n']}-f{c['max_hz']:g}"
        if sha(destination/"entrees"/(nom+".npz")) != c["entree_sha256"]:
            raise ValueError("provenance d'entrée incohérente")
        for config in c["configurations"]:
            for e in config["essais"]:
                log = f"{nom}-{e['variante']}-{e['role']}{e['repetition']}.log"
                if log not in logs or hashlib.sha256(logs[log].encode()).hexdigest() != e["journal_sha256"]:
                    raise ValueError("provenance de journal incohérente")
    print("Archive des réponses groupées conforme : "+str(destination))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--worker", nargs=6)
    p.add_argument("--campagne", nargs=5, metavar=("SORTIE", "ENTREES", "PY010", "PY011", "PYHCB"))
    p.add_argument("--cpu", type=int, default=8)
    p.add_argument("--archiver", nargs=3)
    p.add_argument("--verifier", type=Path)
    a = p.parse_args()
    if a.worker:
        entree, ref, v, modes, sortie, cpu = a.worker
        os.sched_setaffinity(0, {int(cpu)})
        worker(entree, ref, v, int(modes), sortie)
    elif a.campagne:
        campagne(*a.campagne, a.cpu)
    elif a.archiver:
        archiver(*a.archiver)
    elif a.verifier:
        verifier(a.verifier)
    else:
        p.error("choisir une opération")
