"""Keep the static worker and its children tied to the Linux FreeCAD lifetime.

SPDX-License-Identifier: Apache-2.0
Launched as a new session by QProcess. Imports no GUI or CAD library.
"""

import os
import signal
import subprocess
import sys
from pathlib import Path


def run(parent, arguments):
    if os.getpgrp() != os.getpid():
        raise RuntimeError("The static guardian requires its own process group.")
    if os.getppid() != parent:
        raise RuntimeError("FreeCAD closed before the static guardian started.")
    worker = subprocess.Popen(
        [sys.executable, str(Path(__file__).with_name("static_worker.py")), *arguments]
    )
    while True:
        if os.getppid() != parent:
            # This also terminates this guardian; all external children inherit
            # its session and group, even when FreeCAD itself crashes.
            os.killpg(os.getpgrp(), signal.SIGKILL)
        try:
            return worker.wait(timeout=0.1)
        except subprocess.TimeoutExpired:
            pass


if __name__ == "__main__":
    raise SystemExit(run(int(sys.argv[1]), sys.argv[2:]))
