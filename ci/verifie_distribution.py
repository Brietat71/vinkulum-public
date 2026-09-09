"""Vérifie les modèles et notices distribués, sans réseau ni dépendance tierce."""

import hashlib
import json
from pathlib import Path


def verifier(root):
    root = Path(root).resolve()
    base = root / "python/vinkulum/licences"
    modeles = json.loads((base / "modeles.json").read_text())
    dependances = json.loads((base / "rust.json").read_text())

    def empreinte(chemin, attendu):
        path = (root / chemin).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f"Chemin hors distribution : {chemin}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != attendu:
            raise ValueError(f"Empreinte différente : {chemin}")

    attendus = {str(p.relative_to(root)) for p in
                (root / "python/vinkulum/donnees").rglob("*")
                if p.suffix in (".xml", ".urdf")}
    declares = {m["local"] for m in modeles}
    if declares != attendus or len(declares) != len(modeles):
        raise ValueError("Inventaire XML/URDF incomplet ou dupliqué")
    for modele in modeles:
        empreinte(modele["local"], modele["sha256"])
        if not modele["license_sha256"]:
            raise ValueError(f"Notice absente : {modele['local']}")
        for chemin, sha in modele["license_sha256"].items():
            empreinte(chemin, sha)
    for dep in dependances:
        if not dep["license"] or not dep["notice_sha256"]:
            raise ValueError(f"Licence ou notice absente : {dep['name']}")
        for chemin, sha in dep["notice_sha256"].items():
            empreinte(chemin, sha)
    print(f"Distribution : {len(modeles)} modèles et "
          f"{len(dependances)} dépendances, notices et empreintes vérifiées.")


if __name__ == "__main__":
    verifier(Path(__file__).resolve().parents[1])
