"""Prepare the small, editable examples shipped with the desktop application."""

import hashlib
import json
import shutil
import sys
import zipfile
from pathlib import Path


def prepare(destination):
    root = Path(__file__).resolve().parents[1]
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    for name, target in (
        ("platine-parametrique.vinkulum.json", "parametric-plate.vinkulum.json"),
        ("platine-percee.vinkulum.json", "machined-plate.vinkulum.json"),
    ):
        shutil.copyfile(root / "examples/studio" / name, destination / target)
    shutil.copytree(
        root / "examples/studio/articulated/double-pendulum",
        destination / "pinocchio",
    )
    shutil.copytree(
        root / "docs/bancs/studio-tetrahedra-060/calculation",
        destination / "tetra-bending",
    )
    bench = root / "docs/bancs/studio-cad-meshing-060"
    archive = bench / "mesh-and-static-example.zip"
    expected = json.loads((bench / "qualification.json").read_text())["example_sha256"]
    if hashlib.sha256(archive.read_bytes()).hexdigest() != expected:
        raise ValueError(
            "The packaged CAD/mesh example does not match its qualification."
        )
    with zipfile.ZipFile(archive) as stream:
        if any(
            Path(n).is_absolute() or ".." in Path(n).parts for n in stream.namelist()
        ):
            raise ValueError("Invalid archive path.")
        stream.extractall(destination / "cad-statics")
    shutil.copyfile(
        root / "apps/studio/packaging/EXAMPLES.md", destination / "README.md"
    )
    print(f"Prepared delivered examples: {destination}")


if __name__ == "__main__":
    prepare(sys.argv[1])
