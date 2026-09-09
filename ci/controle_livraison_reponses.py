"""Lie la roue livrée au code chronométré et vérifie les exemples publics.

Les roues chronométrées restent intactes. Seules les métadonnées de la
roue livrée peuvent différer : le README final décrit les mesures obtenues.
Usage : python ci/controle_livraison_reponses.py DEPOT ROUE_MESUREE ROUE_LIVREE SORTIE
Exécuter dans l'environnement isolé contenant la roue livrée.
"""
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import zipfile

from controle_roue_reduction import main as controle_roue


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(depot, mesuree, livree, sortie):
    depot, mesuree, livree, sortie = map(Path, (depot, mesuree, livree, sortie))
    controle_roue(depot, livree, sortie)
    r = json.loads(sortie.read_text())
    assert r["version"] == "0.11.0"
    with zipfile.ZipFile(mesuree) as a, zipfile.ZipFile(livree) as b:
        noms = {n for n in a.namelist() if n.startswith("vinkulum/")}
        assert noms == {n for n in b.namelist() if n.startswith("vinkulum/")}
        assert all(a.read(n) == b.read(n) for n in noms), "code mesuré et livré différent"
        empreintes = {n: hashlib.sha256(a.read(n)).hexdigest() for n in sorted(noms)}
    archive = depot/"docs/bancs/reponses-groupees-0.11.0/essais.json.gz"
    campagne = json.loads(gzip.decompress(archive.read_bytes()))
    essais = [e for c in campagne["cas"] for v in c["configurations"]
              if v["variante"] in ("api_011", "lot_011") for e in v["essais"]]
    assert len(essais) == 72
    for e in essais:
        assert e["moteur"] == r["version"]
        assert e["extension_sha256"] == r["roue"]["extension_sha256"]
        assert {"vinkulum/"+k: v for k, v in e["sources_python_sha256"].items()} == r["sources_python_sha256"]
    blocs = re.findall(r"```python\n(.*?)\n```", (depot/"docs/REDUCTION_PORTS.md").read_text(), re.S)
    contexte = {}
    for i in (0, 2):
        exec(compile(blocs[i], "docs/REDUCTION_PORTS.md", "exec"), contexte)
    assert len(contexte["reponses"]) == 3
    r["exemple_charges_groupees"] = {"charges": 3, "superposition_verifiee": True}
    r["concordance_mesure_livraison"] = {
        "roue_mesuree": {"nom": mesuree.name, "sha256": sha(mesuree)},
        "archive_essais_sha256": sha(archive), "essais_identiques": len(essais),
        "fichiers_paquet_identiques_sha256": empreintes,
        "note": "Tous les fichiers vinkulum/ sont identiques ; métadonnées distinctes pour le README final."}
    sortie.write_text(json.dumps(r, ensure_ascii=False, indent=2)+"\n")
    print(f"Code livré identique à la roue mesurée et aux {len(essais)} essais 0.11 ; exemple groupé vérifié.")


if __name__ == "__main__":
    main(*sys.argv[1:])
