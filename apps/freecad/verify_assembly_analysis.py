"""Check actual Motion-task Assembly results against the independent reference.

SPDX-License-Identifier: Apache-2.0
Requires NumPy only; uses the original finite-section double-pendulum equations.
"""

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

from verify_assembly import reference


def verify(root):
    report = json.loads((root / "report.json").read_text())
    assert report["status"] == "passed", report
    yaw = math.radians(17)
    frame = np.array(
        (
            (math.cos(yaw), -math.sin(yaw), 0),
            (math.sin(yaw), math.cos(yaw), 0),
            (0, 0, 1),
        )
    )
    identities, runs, rows = [], [], []
    for key in ("first_capture", "second_capture"):
        directory = root / "runs" / Path(report[key]).name
        raw = (directory / "request.json").read_bytes()
        request = json.loads(raw)
        result = json.loads((directory / "result/playback.json").read_text())
        metadata = json.loads((directory / "result/result.json").read_text())
        assert hashlib.sha256(raw).hexdigest() == result["request_sha256"]
        assert request["run_id"] == result["run_id"] == metadata["run_id"]
        assert (
            request["density_kg_m3"] == 2700
            and request["duration_s"] == 0.5
            and request["step_s"] == 0.005
        )
        for body in request["bodies"]:
            assert (
                hashlib.sha256((directory / body["step_file"]).read_bytes()).hexdigest()
                == body["step_sha256"]
            )
        ids = [body["body_id"] for body in request["bodies"]]
        assert ids == result["body_ids"]
        identities.append((request["project_id"], ids))
        runs.append(request["run_id"])
        names = [body["name"] for body in request["bodies"]]
        assert set(names) == {"Upper", "Lower"}
        path = directory / "result/trajectory.npz"
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest() == metadata["archive_sha256"]
        )
        with np.load(path, allow_pickle=False) as archive:
            for external, native in (
                ("time_s", "time"),
                ("position_m", "position"),
                ("rotation", "rotation"),
            ):
                np.testing.assert_array_equal(result[external], archive[native])
        positions = (np.array(result["position_m"]) - (0.125, -0.075, 0.210)) @ frame
        upper = positions[:, names.index("Upper")]
        lower = positions[:, names.index("Lower")] - 2 * upper
        angles = np.column_stack(
            [np.arctan2(p[:, 0], -p[:, 2]) for p in (upper, lower)]
        )
        coarse, fine = (reference(result["time_s"], n) for n in (16, 32))
        gap = float(np.max(abs(coarse - fine)))
        error = float(np.max(abs(angles - fine)))
        assert gap < 1e-10 and error < 2e-3, (gap, error)
        rows.append(
            {
                "capture": directory.name,
                "max_angle_error_rad": error,
                "reference_refinement_gap_rad": gap,
            }
        )
    assert identities[0] == identities[1] and runs[0] != runs[1]
    assert max(report["display_errors_m"]) < 1e-9
    result = {
        "status": "passed",
        "trials": rows,
        "stable_project_and_body_ids": True,
        "distinct_run_ids": True,
        "max_display_error_m": max(report["display_errors_m"]),
        "reference": "Independent finite-section double-pendulum Lagrange equations, RK4 with 16/32 subdivisions",
        "imports_absent": [
            name
            for name in (
                "FreeCAD",
                "vinkulum",
                "vinkulum_studio",
                "OCP",
                "build123d",
                "PySide6",
            )
            if importlib.util.find_spec(name) is None
        ],
    }
    (root / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    print(json.dumps(verify(parser.parse_args().root), indent=2))
