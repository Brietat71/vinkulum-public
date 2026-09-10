"""Build the self-contained FreeCAD Mod archive. SPDX-License-Identifier: Apache-2.0"""

import argparse
import hashlib
import json
import zipfile
from pathlib import Path


def package(destination):
    root = Path(__file__).resolve().parent
    names = (
        "InitGui.py",
        "vinkulum_freecad.py",
        "bridge.py",
        "worker.py",
        "vinkulum.svg",
    )
    files = {name: (root / name).read_bytes() for name in names}
    files["LICENSE"] = (root.parent.parent / "LICENSE").read_bytes()
    files["INSTALL.txt"] = (
        "Vinkulum for FreeCAD 0.1.0a1\n\n"
        "Extract the Vinkulum directory into FreeCAD's user Mod directory and restart FreeCAD.\n"
        "Select the Vinkulum workbench. Configure its external Vinkulum Python engine, then\n"
        "create the parametric pendulum example and select Analyse.\n\n"
        "Instructions and scope: https://github.com/Brietat71/vinkulum-public/tree/main/apps/freecad\n"
        "Source modelling stays in FreeCAD. This first workbench supports one rigid solid\n"
        "with one explicit world-frame revolute joint. Results retain their capture directory.\n"
    ).encode()
    manifest = {
        "version": "0.1.0a1",
        "sha256": {
            name: hashlib.sha256(data).hexdigest() for name, data in files.items()
        },
    }
    files["build-info.json"] = (json.dumps(manifest, indent=2) + "\n").encode()
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(destination, "x", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(files.items()):
            info = zipfile.ZipInfo("Vinkulum/" + name, (2026, 9, 10, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, data)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    print(package(args.destination))
