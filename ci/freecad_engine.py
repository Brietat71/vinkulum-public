"""Install and verify a separate FreeCAD engine environment on Linux x86-64.

SPDX-License-Identifier: Apache-2.0
"""

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path


def install(destination, python="3.14"):
    root = Path(__file__).resolve().parents[1]
    destination = destination.expanduser().absolute()
    if destination.exists() or destination.is_symlink():
        raise ValueError(
            f"Choose a new directory; existing files are retained: {destination}"
        )
    if sys.platform != "linux" or platform.machine() != "x86_64":
        raise ValueError("This installation recipe is qualified only on Linux x86-64.")
    commands = {}
    for name in ("uv", "cargo", "patch", "git"):
        executable = shutil.which(name)
        if executable is None:
            raise ValueError(f"Install {name} before creating the engine environment.")
        commands[name] = executable
    environment = os.environ.copy()
    for name in (
        "PYTHONHOME",
        "PYTHONPATH",
        "VIRTUAL_ENV",
        "LD_LIBRARY_PATH",
        "QT_PLUGIN_PATH",
        "QT_QPA_PLATFORM_PLUGIN_PATH",
    ):
        environment.pop(name, None)
    environment.setdefault("UV_CACHE_DIR", str(destination / "cache"))
    report = {
        "status": "installing",
        "source_commit": subprocess.check_output(
            [commands["git"], "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "source_dirty": bool(
            subprocess.check_output(
                [commands["git"], "status", "--porcelain"], cwd=root
            )
        ),
        "platform": platform.platform(),
        "requested_python": python,
        "installer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "steps": [],
    }
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(__file__, destination / "installer.py")
    report_path = destination / "installation.json"

    def save():
        report_path.write_text(json.dumps(report, indent=2) + "\n")

    def run(label, args):
        report["steps"].append({"label": label, "command": [str(a) for a in args]})
        save()
        print(label, flush=True)
        with (destination / "installation.log").open("ab") as log:
            log.write(("\n" + label + "\n").encode())
            log.flush()
            subprocess.run(
                args,
                cwd=root,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        report["steps"][-1]["status"] = "passed"
        save()

    interpreter = destination / "venv/bin/python"
    sources = destination / "cad-sources"
    try:
        run(
            "Create isolated Python environment",
            [commands["uv"], "venv", "--python", python, str(destination / "venv")],
        )
        run(
            "Verify Python requirement",
            [
                str(interpreter),
                "-c",
                "import sys; assert sys.version_info >= (3,14), 'Python 3.14 or newer is required'",
            ],
        )
        run(
            "Check upstream source hashes and apply reviewed OCCT 8 patches",
            [str(interpreter), str(root / "ci/prepare_cad.py"), str(sources)],
        )
        run(
            "Build and install the native kernel and CAD/mechanical adapters",
            [
                commands["uv"],
                "pip",
                "install",
                "--python",
                str(interpreter),
                str(root),
                str(sources / "build123d-0.11.1"),
                str(sources / "ocpsvg-0.6.0"),
                str(root / "apps/studio") + "[cad,test]",
            ],
        )
        run(
            "Check installed dependency constraints",
            [commands["uv"], "pip", "check", "--python", str(interpreter)],
        )
        run(
            "Verify OCCT 8 and installed runtime versions",
            [
                str(interpreter),
                "-c",
                (
                    "import importlib.metadata as m,json,sys,OCP,build123d,vinkulum; "
                    "from pathlib import Path; "
                    "assert OCP.__version__.split('.')[0] == '8'; "
                    "Path(sys.argv[1]).write_text(json.dumps({'python':sys.version,"
                    "'ocp_binding':OCP.__version__,'packages':"
                    "{d.metadata['Name']:d.version for d in m.distributions()}},sort_keys=True,indent=2))"
                ),
                str(destination / "runtime.json"),
            ],
        )
        report["runtime"] = json.loads((destination / "runtime.json").read_text())
        run(
            "Verify installed engine adapters without Qt or VTK",
            [
                str(interpreter),
                str(root / "ci/verify_engine_headless.py"),
                str(destination / "headless.json"),
            ],
        )
        run(
            "Check actual mechanics against the independent pendulum and reject corrupted mass",
            [
                str(interpreter),
                "-m",
                "unittest",
                "discover",
                "-s",
                str(root / "apps/studio/tests"),
                "-p",
                "test_freecad_bridge.py",
                "-v",
            ],
        )
        report["engine_python"] = str(interpreter)
        report["status"] = "passed"
        save()
        print(f"\nFreeCAD → Motion analysis → Engine Python:\n{interpreter}")
        print(f"Installation record: {report_path}")
    except Exception as error:
        report["status"] = "failed"
        report["error"] = str(error)
        save()
        raise RuntimeError(
            f"Installation failed; inspect {destination / 'installation.log'}. "
            "The directory was retained for diagnosis."
        ) from error


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "destination",
        type=Path,
        help="New directory; never modifies an existing environment",
    )
    parser.add_argument(
        "--python", default="3.14", help="Python version or executable understood by uv"
    )
    args = parser.parse_args()
    try:
        install(args.destination, args.python)
    except (ValueError, RuntimeError) as error:
        parser.exit(1, str(error) + "\n")
