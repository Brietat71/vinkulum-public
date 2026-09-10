"""Deterministic process barriers for FreeCAD lifecycle qualification only.

The released workbench never invokes this helper. The real worker executes in
this same process when the before-calculation barrier is released.
SPDX-License-Identifier: Apache-2.0
"""

import json
import os
import runpy
import sys
import time
from pathlib import Path

worker, root = Path(sys.argv[1]), Path(sys.argv[2])
stage = json.loads((root / "gate.json").read_text())["stage"]


def barrier(name):
    if stage != name:
        return
    marker = root / "entered.tmp"
    marker.write_text(json.dumps({"stage": name, "pid": os.getpid()}))
    marker.replace(root / "entered.json")
    deadline = time.monotonic() + 15
    while not (root / "release").exists():
        if time.monotonic() > deadline:
            raise RuntimeError("Qualification barrier expired")
        time.sleep(0.01)


barrier("before")
sys.argv = [str(worker), str(root)]
runpy.run_path(str(worker), run_name="__main__")
barrier("after")
