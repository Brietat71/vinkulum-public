"""Vérifie une roue installée hors dépôt et rejoue les exemples documentés.

Usage : python ci/controle_roue_reduction.py DEPOT ROUE SORTIE_JSON
L'interpréteur doit désigner l'environnement neuf recevant cette roue.
"""
import base64
import csv
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path
import re
import sys
import tomllib
import zipfile
from email.parser import BytesParser

import vinkulum
from vinkulum import _vinkulum
from vinkulum import reduction_ports


def sha(contenu):
    return hashlib.sha256(contenu).hexdigest()


def main(depot, roue, sortie):
    depot, roue, sortie = map(Path, (depot, roue, sortie))
    assert vinkulum.depot_docs() is None, "exécuter avec une installation hors dépôt"
    version = importlib.metadata.version("vinkulum")
    assert version == vinkulum.__version__
    assert version == tomllib.loads((depot/"Cargo.toml").read_text())["package"]["version"]
    assert version == tomllib.loads((depot/"pyproject.toml").read_text())["project"]["version"]
    assert hasattr(vinkulum.Noyau, "facteurs_materiels_poutres")
    with zipfile.ZipFile(roue) as archive:
        metadata = BytesParser().parsebytes(archive.read(f"vinkulum-{version}.dist-info/METADATA"))
        assert metadata["Version"] == version
        assert "reduction" in metadata.get_all("Provides-Extra", [])
        readme = (depot/"README.md").read_text()
        assert metadata.get_payload(decode=True).decode("utf8").rstrip("\n") == readme.rstrip("\n")
        record = archive.read(f"vinkulum-{version}.dist-info/RECORD").decode()
        entrees = 0
        for nom, empreinte, taille in csv.reader(io.StringIO(record)):
            if not empreinte:
                continue
            contenu = archive.read(nom)
            algo, attendu = empreinte.split("=", 1)
            assert algo == "sha256"
            actuel = base64.urlsafe_b64encode(hashlib.sha256(contenu).digest()).decode().rstrip("=")
            assert actuel == attendu and len(contenu) == int(taille), nom
            entrees += 1
        sources = {}
        for nom in archive.namelist():
            if nom.startswith("vinkulum/") and nom.endswith(".py"):
                contenu = archive.read(nom)
                assert contenu == (depot/"python"/nom).read_bytes(), nom
                sources[nom] = sha(contenu)
        extension = Path(_vinkulum.__file__)
        assert archive.read("vinkulum/"+extension.name) == extension.read_bytes()
    resultats = {}
    for nom in ("README.md", "docs/REDUCTION_PORTS.md"):
        texte = (depot/nom).read_text()
        exemple = re.findall(r"```python\n(.*?)\n```", texte, re.S)[0]
        contexte = {}
        exec(compile(exemple, nom, "exec"), contexte)
        if nom == "README.md":
            resultats[nom] = {"angle_rad": contexte["angle"]}
        else:
            r, rep, bilan = (contexte[k] for k in ("r", "rep", "bilan"))
            assert bilan["statut"] == "tolerance_aux_frequences_demandees"
            assert r.certificat_spectral["certification_machine"] is True
            assert rep["certification_machine"] is False
            assert rep["deformation"]["borne_relative"] <= 1e-8
            resultats[nom] = {
                "statut": bilan["statut"], "uy_m": float(rep["deplacement_ports_physique"][1]),
                "borne_relative_deformation": rep["deformation"]["borne_relative"],
                "lambda_min_s_2": r.certificat_spectral["lambda_min"],
                "spectre_certifie": True, "champ_certifie": False}
    resultat = {"version": version, "module": str(extension),
                "module_reduction": reduction_ports.__file__, "installation_hors_depot": True,
                "roue": {"nom": roue.name, "sha256": sha(roue.read_bytes()),
                         "octets": roue.stat().st_size, "entrees_record_verifiees": entrees,
                         "readme_conforme": True, "extension_sha256": sha(extension.read_bytes())},
                "sources_python_sha256": sources, "exemples": resultats}
    sortie.write_text(json.dumps(resultat, indent=2, ensure_ascii=False)+"\n")
    print(f"Roue {version} conforme ; {entrees} entrées RECORD ; exemples exécutés.")


if __name__ == "__main__":
    main(*sys.argv[1:])
