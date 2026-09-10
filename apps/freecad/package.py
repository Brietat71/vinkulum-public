"""Build the FreeCAD extension ZIP without bundling FreeCAD or a solver runtime."""

import argparse
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path

VERSION = "0.1.0a3.dev6"


def build(destination, require_clean=False):
    source = Path(__file__).resolve().parent
    root = source.parents[1]
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=root))
    if require_clean and dirty:
        raise ValueError("Commit the reviewed extension before a release build")
    modules = (
        "host.py",
        "analysis.py",
        "bridge.py",
        "worker.py",
        "assembly_capture.py",
        "assembly_worker.py",
        "static_analysis.py",
        "static_examples.py",
        "static_bridge.py",
        "static_host.py",
        "static_job.py",
        "static_guard.py",
        "static_worker.py",
    )
    names = ("InitGui.py", *modules, "INSTALLATION.txt")
    examples = {}
    for archive_name, member, name in (
        (
            "freecad-bridge-2026",
            "record/length-800-step-0.005/source.FCStd",
            "Pendulum.FCStd",
        ),
        ("freecad-assembly-2026", "record/assembly-2.FCStd", "DoublePendulum.FCStd"),
    ):
        with zipfile.ZipFile(
            root / "docs/bancs" / archive_name / "record.zip"
        ) as record:
            examples["Examples/" + name] = record.read(member)
    report = {
        "extension_version": VERSION,
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "source_dirty": dirty,
        "files_sha256": {
            name: hashlib.sha256((source / name).read_bytes()).hexdigest()
            for name in names
        }
        | {name: hashlib.sha256(data).hexdigest() for name, data in examples.items()},
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
        for name in modules:
            archive.write(source / name, "Vinkulum/vinkulum_freecad/" + name)
        for name in ("LICENSE", "NOTICE"):
            archive.write(root / name, "Vinkulum/" + name)
        for name, data in examples.items():
            archive.writestr("Vinkulum/vinkulum_freecad/" + name, data)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--require-clean", action="store_true")
    args = parser.parse_args()
    build(args.output, args.require_clean)
