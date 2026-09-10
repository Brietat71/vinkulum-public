"""Qualify an installed extension in an isolated real FreeCAD Linux process."""

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import zipfile
from pathlib import Path


def qualify(freecad, engine_python, output, archive=None, recipe="extension"):
    root = Path(__file__).resolve().parents[1]
    output = output.absolute()
    output.mkdir(parents=True, exist_ok=False)
    payload = output / "Vinkulum-FreeCAD.zip"
    if archive is None:
        subprocess.run(
            [sys.executable, str(root / "apps/freecad/package.py"), str(payload)],
            check=True,
        )
    else:
        shutil.copyfile(archive, payload)
    module_root = output / "data/FreeCAD/Mod"
    module_root.mkdir(parents=True)
    with zipfile.ZipFile(payload) as contents:
        for entry in contents.infolist():
            path = Path(entry.filename)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("Extension archive has an unsafe path")
        contents.extractall(module_root)
    barrier = output / "barrier-python"
    barrier.write_text(
        '#!/bin/sh\n: > "$VINKULUM_FREECAD_EXTENSION_CHECK/barrier-entered"\n'
        "while :; do sleep 0.05; done\n"
    )
    barrier.chmod(0o755)
    environment = os.environ.copy()
    for key in (
        "PYTHONHOME",
        "PYTHONPATH",
        "VIRTUAL_ENV",
        "QT_PLUGIN_PATH",
        "QT_QPA_PLATFORM_PLUGIN_PATH",
    ):
        environment.pop(key, None)
    environment.update(
        {
            "XDG_CONFIG_HOME": str(output / "config"),
            "XDG_DATA_HOME": str(output / "data"),
            "XDG_CACHE_HOME": str(output / "cache"),
            "VINKULUM_FREECAD_EXTENSION_CHECK": str(output),
            "VINKULUM_CLOSE_CAPTURE": str(output / "example-close"),
            "VINKULUM_FREECAD_PYTHON": str(engine_python.absolute()),
            "VINKULUM_FREECAD_BARRIER_PYTHON": str(barrier),
            "VINKULUM_FREECAD_OUTPUT_ROOT": str(output / "runs"),
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
        }
    )
    command = [
        str(freecad.absolute()),
        str(root / "apps/freecad" / f"qualify_{recipe.replace('-', '_')}.FCMacro"),
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
            code = process.wait(timeout=130)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait(timeout=5)
            raise
    if code:
        raise RuntimeError(f"FreeCAD exited with {code}; see {output / 'console.log'}")
    report_path = (
        output / "extension-check.json"
        if recipe == "extension"
        else output / "example-close/report.json"
    )
    report = json.loads(report_path.read_text())
    if report["status"] != "passed":
        raise RuntimeError(f"FreeCAD qualification failed: {report}")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freecad", type=Path, required=True)
    parser.add_argument("--engine-python", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument(
        "--recipe", choices=("extension", "example-close"), default="extension"
    )
    args = parser.parse_args()
    qualify(args.freecad, args.engine_python, args.output, args.archive, args.recipe)
