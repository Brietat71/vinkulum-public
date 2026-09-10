"""Build the FreeCAD extension ZIP without bundling FreeCAD or a solver runtime."""

import argparse
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path

VERSION = "0.1.0a1"


def build(destination, require_clean=False):
    source = Path(__file__).resolve().parent
    root = source.parents[1]
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=root))
    if require_clean and dirty:
        raise ValueError("Commit the reviewed extension before a release build")
    names = ("InitGui.py", "host.py", "bridge.py", "worker.py", "INSTALLATION.txt")
    report = {
        "extension_version": VERSION,
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "source_dirty": dirty,
        "files_sha256": {
            name: hashlib.sha256((source / name).read_bytes()).hexdigest()
            for name in names
        },
        "runtime_scope": "FreeCAD 1.1.3 / Qt 6 on Linux; separate Vinkulum 0.20 / OCCT 8 Python environment",
    }
    with zipfile.ZipFile(destination, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("Vinkulum/Init.py", '"""Vinkulum FreeCAD extension."""\n')
        archive.write(source / "InitGui.py", "Vinkulum/InitGui.py")
        archive.write(source / "INSTALLATION.txt", "Vinkulum/INSTALLATION.txt")
        archive.writestr(
            "Vinkulum/extension-info.json", json.dumps(report, indent=2) + "\n"
        )
        archive.writestr(
            "Vinkulum/vinkulum_freecad/__init__.py",
            f'"""FreeCAD host for Vinkulum."""\n__version__ = "{VERSION}"\n',
        )
        for name in ("host.py", "bridge.py", "worker.py"):
            archive.write(source / name, "Vinkulum/vinkulum_freecad/" + name)
        for name in ("LICENSE", "NOTICE"):
            archive.write(root / name, "Vinkulum/" + name)
        with zipfile.ZipFile(
            root / "docs/bancs/freecad-bridge-2026/record.zip"
        ) as record:
            archive.writestr(
                "Vinkulum/vinkulum_freecad/Examples/Pendulum.FCStd",
                record.read("record/length-800-step-0.005/source.FCStd"),
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--require-clean", action="store_true")
    args = parser.parse_args()
    build(args.output, args.require_clean)
