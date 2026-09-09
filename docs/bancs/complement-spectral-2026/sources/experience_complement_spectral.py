"""Qualification du complément certifié et témoin KKT sur la bande 0–40 Hz.

Une exécution par cas : les durées sont indicatives, sans classement.
Les trois jeux natifs, directions choisies et oracles préexistants sont
identifiés. Le certificat spectral et le juge du champ restent distincts.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time

import numpy as np
import scipy
from scipy.sparse import csr_matrix
import vinkulum
from vinkulum._ports.condensation_energie import CondensationEnergie
from trace_complement_dirigee import certifier_complement
from retention_complement import TemoinRetention
from confronte_ports_exudyn import lire, juger, THREADS, distribution_figee

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ("trace_complement_dirigee.py", "retention_complement.py",
           "experience_complement_spectral.py", "confronte_ports_exudyn.py")


def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for bloc in iter(lambda: f.read(1024*1024), b""):
            h.update(bloc)
    return h.hexdigest()


def ecrire(p, r):
    Path(p).write_text(json.dumps(r, indent=2, ensure_ascii=False, allow_nan=False)+"\n")


def mesurer(entrees, directions, sortie, cpu):
    entrees, directions, sortie = map(Path, (entrees, directions, sortie))
    sortie.mkdir(parents=True, exist_ok=True)
    if {k: os.environ.get(k) for k in THREADS} != THREADS:
        raise ValueError("un fil par bibliothèque requis")
    os.sched_setaffinity(0, {cpu})
    paquet = Path(vinkulum.__file__).parent
    rapport = dict(date_utc=datetime.now(timezone.utc).isoformat(),
        python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
        vinkulum=vinkulum.__version__, distribution=distribution_figee("vinkulum"),
        extension_sha256=sha(next(paquet.glob("*.so"))),
        sources_python_sha256={str(p.relative_to(paquet)): sha(p) for p in sorted(paquet.rglob("*.py"))},
        sources_sha256={n: sha(ROOT/"ci"/n) for n in SOURCES},
        cpu=sorted(os.sched_getaffinity(0)), fils=THREADS,
        protocole="une exécution, durées indicatives ; aucune comparaison de vitesse",
        cas=[])
    for n in (32, 128, 512):
        fichier = entrees/f"n{n}-f40.npz"
        ref = entrees/f"n{n}-f40.ref.npy"
        dirfile = directions/f"n{n}-directions.npz"
        d, m, metrique, forces, omegas = lire(fichier)
        i, s = np.arange(d.shape[1]-6), np.arange(d.shape[1]-6, d.shape[1])
        c = dict(n=n, entree_sha256=sha(fichier), reference_sha256=sha(ref),
            directions_sha256=sha(dirfile), reference=json.loads(ref.with_suffix(ref.suffix+".json").read_text()),
            frequences=257, charges=6, bande_hz=[0, 40], certificats=[])
        with np.load(dirfile, allow_pickle=False) as a:
            r = csr_matrix((a["r_data"], a["r_indices"], a["r_indptr"]), shape=tuple(a["r_shape"]))
            e, b, phi = a["echelles"].copy(), a["b"].copy(), a["phi"].copy()
        for k in (1, 2):
            cert = certifier_complement(d[:, i], r, e, m[i][:, i], b[:, :k],
                                        budget_operations=10_000_000)
            cert["bande_demandee_binary64_couverte"] = Fraction(cert["lambda_min"]) > Fraction(float(omegas[-1]))**2
            cert["omega_max_binary64"] = float(omegas[-1])
            c["certificats"].append(cert)
        print(n, "certificats", [v["lambda_min"] for v in c["certificats"]], flush=True)
        if not all(v["bande_demandee_binary64_couverte"] for v in c["certificats"]):
            raise ValueError("bande non couverte par le certificat")
        debut = time.perf_counter()
        qr = CondensationEnergie(d, i, s, metrique)
        temoin = TemoinRetention(d, m, i, s, b[:, :1], phi[:, :1],
            c["certificats"][0]["lambda_min"], float(omegas[-1]),
            psi=qr.psi, normalisation=qr.w)
        c["preparation_temoin_indicative_s"] = time.perf_counter()-debut
        champs = np.empty((len(omegas), d.shape[1], 6))
        diagnostics = {k: [] for k in ("sigma_min_bloc_conserve", "defaut_contrainte",
                                       "residu_kkt", "ecart_fonctionnel")}
        debut = time.perf_counter()
        for j, w in enumerate(omegas):
            rep = temoin.reponses(w, forces)
            champs[j] = rep["champ"]
            for k in diagnostics:
                diagnostics[k].append(rep[k])
        c.update(reponses_temoin_indicative_s=time.perf_counter()-debut,
                 taille_conservee=rep["taille_conservee"], taille_kkt=rep["taille_kkt"],
                 diagnostics=diagnostics, certification_champ_machine=False)
        c["juge"] = juger(d, m, metrique, champs, np.load(ref, mmap_mode="r", allow_pickle=False))
        rapport["cas"].append(c)
        ecrire(sortie/"rapport.json", rapport)
        print(n, "juge", c["juge"]["accepte"], c["juge"]["maxima_operateurs"], flush=True)


def archiver(entrees, directions, sortie, destination):
    entrees, directions, sortie, destination = map(Path, (entrees, directions, sortie, destination))
    r = json.loads((sortie/"rapport.json").read_text())
    if [c["n"] for c in r["cas"]] != [32, 128, 512]:
        raise ValueError("expérience incomplète")
    destination.mkdir(parents=True, exist_ok=True)
    (destination/"rapport.json.gz").write_bytes(gzip.compress(json.dumps(r, ensure_ascii=False, allow_nan=False).encode(), mtime=0))
    for folder in ("sources", "entrees", "directions"):
        (destination/folder).mkdir(exist_ok=True)
    for nom, h in r["sources_sha256"].items():
        if sha(ROOT/"ci"/nom) != h:
            raise ValueError("source modifiée depuis l'expérience")
        shutil.copy2(ROOT/"ci"/nom, destination/"sources"/nom)
    for c in r["cas"]:
        for source, folder, nom, key in (
            (entrees, "entrees", f"n{c['n']}-f40.npz", "entree_sha256"),
            (directions, "directions", f"n{c['n']}-directions.npz", "directions_sha256")):
            if sha(source/nom) != c[key]:
                raise ValueError("donnée modifiée depuis l'expérience")
            shutil.copy2(source/nom, destination/folder/nom)
    ecrire(destination/"manifest.json", dict(
        version_support="0.11.0", prototype="complément spectral, hors API de production",
        fichiers_sha256={str(p.relative_to(destination)): sha(p) for p in sorted(destination.rglob("*"))
                         if p.is_file() and p.name != "manifest.json"}))


def verifier(destination):
    destination = Path(destination)
    m = json.loads((destination/"manifest.json").read_text())
    hashes = m["fichiers_sha256"]
    attendus = {"rapport.json.gz"} | {"sources/"+n for n in SOURCES}
    attendus |= {f"entrees/n{n}-f40.npz" for n in (32, 128, 512)}
    attendus |= {f"directions/n{n}-directions.npz" for n in (32, 128, 512)}
    presents = {str(p.relative_to(destination)) for p in destination.rglob("*")
                if p.is_file() and p.name != "manifest.json"}
    if set(hashes) != presents or not attendus <= presents or any(sha(destination/n) != h for n,h in hashes.items()):
        raise ValueError("archive incomplète ou altérée")
    r = json.loads(gzip.decompress((destination/"rapport.json.gz").read_bytes()))
    if len(r["cas"]) != 3 or [c["n"] for c in r["cas"]] != [32, 128, 512]:
        raise ValueError("couverture de maillages incomplète")
    if set(r["sources_sha256"]) != set(SOURCES):
        raise ValueError("sources incomplètes")
    for n, h in r["sources_sha256"].items():
        if hashes["sources/"+n] != h:
            raise ValueError("provenance source différente")
    for c in r["cas"]:
        if (hashes[f"entrees/n{c['n']}-f40.npz"] != c["entree_sha256"]
                or hashes[f"directions/n{c['n']}-directions.npz"] != c["directions_sha256"]):
            raise ValueError("provenance entrée différente")
        if len(c["certificats"]) != 2 or [v["rang_contraintes"] for v in c["certificats"]] != [1, 2]:
            raise ValueError("couverture des rangs incomplète")
        for v in c["certificats"]:
            if not v["certification_machine"] or not v["bande_demandee_binary64_couverte"]:
                raise ValueError("certificat absent")
            if Fraction(v["lambda_min"]) <= Fraction(v["omega_max_binary64"])**2:
                raise ValueError("bande non couverte")
        juge = c["juge"]
        for nom in ("masse", "deformation", "port"):
            vals, ops = juge["erreurs"][nom], juge["operateurs"][nom]
            if len(vals) != 257 or any(len(x) != 6 for x in vals) or len(ops) != 257:
                raise ValueError("couverture du juge incomplète")
            if max(x for v in vals for x in v) != juge["maxima"][nom] or max(ops) != juge["maxima_operateurs"][nom]:
                raise ValueError("maxima incohérents")
        bon = max(*juge["maxima"].values(), *juge["maxima_operateurs"].values()) <= 1e-6
        if bool(juge["accepte"]) != bon:
            raise ValueError("qualification incohérente")
    print("Archive du complément conforme : 6 certificats et 3 bandes complètes")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesurer", nargs=3)
    p.add_argument("--cpu", type=int, default=8)
    p.add_argument("--archiver", nargs=4)
    p.add_argument("--verifier", type=Path)
    a = p.parse_args()
    if a.mesurer: mesurer(*a.mesurer, a.cpu)
    elif a.archiver: archiver(*a.archiver)
    elif a.verifier: verifier(a.verifier)
    else: p.error("choisir une opération")
