"""Run the isolated FreeCAD/OCCT8/Gmsh/CalculiX transfer experiment."""

import argparse
import hashlib
import json
import os
import signal
import subprocess
import zipfile
from pathlib import Path


def qualify(args):
    if args.native_inputs and args.recipe != "static":
        raise ValueError("--native-inputs requires --recipe static.")
    root = Path(__file__).resolve().parents[1]
    output = args.output.absolute()
    output.mkdir(parents=True, exist_ok=False)
    sources = [
        "ci/freecad_static.py",
        "apps/freecad/static_bridge.py",
        "apps/freecad/native_static.py",
        "apps/freecad/static_worker.py",
        "apps/freecad/qualify_static.FCMacro",
        "apps/freecad/bridge.py",
    ]
    if args.recipe in ("static-task", "static-fingerprint"):
        sources += [
            "apps/freecad/qualify_static_task.FCMacro",
            "apps/freecad/static_host.py",
            "apps/freecad/static_job.py",
            "apps/freecad/static_guard.py",
            "apps/freecad/static_analysis.py",
            "apps/freecad/static_examples.py",
            "apps/freecad/host.py",
            "apps/freecad/package.py",
        ]
        if args.recipe == "static-fingerprint":
            sources.append("apps/freecad/qualify_static_fingerprint.FCMacro")
        payload = output / "Vinkulum-FreeCAD.zip"
        subprocess.run(
            ["python3", str(root / "apps/freecad/package.py"), str(payload)], check=True
        )
        with zipfile.ZipFile(payload) as archive:
            archive.extractall(output / "data/FreeCAD/Mod")
    provenance = {
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "source_status": subprocess.check_output(
            ["git", "status", "--short"], cwd=root, text=True
        ),
        "source_sha256": {
            name: hashlib.sha256((root / name).read_bytes()).hexdigest()
            for name in sources
        },
        "executables": {
            name: str(getattr(args, name).absolute())
            for name in ("freecad", "engine_python", "gmsh", "ccx")
        },
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    environment = os.environ.copy()
    for name in (
        "PYTHONHOME",
        "PYTHONPATH",
        "VIRTUAL_ENV",
        "QT_PLUGIN_PATH",
        "QT_QPA_PLATFORM_PLUGIN_PATH",
    ):
        environment.pop(name, None)
    environment.update(
        {
            "XDG_CONFIG_HOME": str(output / "config"),
            "XDG_DATA_HOME": str(output / "data"),
            "XDG_CACHE_HOME": str(output / "cache"),
            "VINKULUM_STATIC_NATIVE_INPUTS": "1" if args.native_inputs else "0",
            "VINKULUM_STATIC_SOURCE": str(root),
            "VINKULUM_STATIC_CHECK": str(output),
            "VINKULUM_STATIC_PYTHON": str(args.engine_python.absolute()),
            "VINKULUM_STATIC_GMSH": str(args.gmsh.absolute()),
            "VINKULUM_STATIC_CCX": str(args.ccx.absolute()),
        }
    )
    command = [
        str(args.freecad.absolute()),
        str(
            root
            / "apps/freecad"
            / ("qualify_" + args.recipe.replace("-", "_") + ".FCMacro")
        ),
    ]
    if not environment.get("DISPLAY"):
        command = ["xvfb-run", "-a", "-s", "-screen 0 1600x1100x24", *command]
    with (output / "console.log").open("wb") as log:
        process = subprocess.Popen(
            command,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            code = process.wait(timeout=270)
        finally:
            # Kill the dedicated host group. Static-job guardians detect its death
            # and kill their separately owned worker/mesher/solver groups.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait(timeout=5)
    if code:
        raise RuntimeError(f"FreeCAD exited with {code}; inspect {output}")
    if "Traceback (most recent call last)" in (output / "console.log").read_text(
        errors="replace"
    ):
        raise RuntimeError(
            f"Native Python callback failed; inspect {output / 'console.log'}"
        )
    report = json.loads((output / "static-check.json").read_text())
    if report["status"] != "passed":
        raise RuntimeError(report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("freecad", "engine-python", "gmsh", "ccx", "output"):
        parser.add_argument("--" + option, type=Path, required=True)
    parser.add_argument("--native-inputs", action="store_true")
    parser.add_argument(
        "--recipe",
        choices=("static", "static-task", "static-fingerprint"),
        default="static",
    )
    qualify(parser.parse_args())
