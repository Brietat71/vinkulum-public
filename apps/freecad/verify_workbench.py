"""Independent checks of installed-workbench captures; requires only NumPy.

SPDX-License-Identifier: Apache-2.0
"""

import argparse
import hashlib
import importlib.util
import json
import math
import zipfile
from pathlib import Path

import numpy as np

from verify import reference


def verify(root):
    report = json.loads((root / "report.json").read_text())
    assert report["status"] == "passed", report
    assert [trial["axis"] for trial in report["trials"]] == ["Y", "X", "Z"]
    rows = []
    for trial in report["trials"]:
        directory = root / trial["directory"]
        request_bytes = (directory / "request.json").read_bytes()
        request = json.loads(request_bytes)
        result = json.loads((directory / "result/playback.json").read_text())
        metadata = json.loads((directory / "result/result.json").read_text())
        assert hashlib.sha256(request_bytes).hexdigest() == result["request_sha256"]
        assert (
            hashlib.sha256((directory / "part.step").read_bytes()).hexdigest()
            == request["step_sha256"]
        )
        assert (
            hashlib.sha256(
                (directory / "result/trajectory.npz").read_bytes()
            ).hexdigest()
            == metadata["archive_sha256"]
        )
        assert request["run_id"] == result["run_id"] == metadata["run_id"]
        length = trial["length_mm"] * 0.001
        mass = 0.020 * 0.030 * length * 7800
        pivot = np.array(trial["pivot_mm"]) * 0.001
        np.testing.assert_array_equal(request["pivot_m"], pivot)
        axis = np.eye(3)["XYZ".index(trial["axis"])]
        np.testing.assert_array_equal(request["axis_world"], axis)
        a, yaw = math.radians(160), math.radians(trial["yaw_deg"])
        ry = np.array(
            ((math.cos(a), 0, math.sin(a)), (0, 1, 0), (-math.sin(a), 0, math.cos(a)))
        )
        rz = np.array(
            (
                (math.cos(yaw), -math.sin(yaw), 0),
                (math.sin(yaw), math.cos(yaw), 0),
                (0, 0, 1),
            )
        )
        orientation = rz @ ry
        centre = pivot + orientation @ (0, 0, length / 2)
        inertia = (
            orientation
            @ np.diag(
                mass
                * np.array(
                    (0.030**2 + length**2, 0.020**2 + length**2, 0.020**2 + 0.030**2)
                )
                / 12
            )
            @ orientation.T
        )
        for properties in (request["properties_si"], result["properties_si"]):
            for name, expected, scale in (
                ("volume_m3", mass / 7800, mass / 7800),
                ("mass_kg", mass, mass),
                ("centre_m", centre, length),
                ("inertia_kg_m2", inertia.ravel(), mass * length**2),
            ):
                np.testing.assert_allclose(
                    properties[name], expected, atol=1e-8 * scale, rtol=0
                )
        times, positions, rotations = (
            np.array(result[key]) for key in ("time_s", "position_m", "rotation")
        )
        with np.load(directory / "result/trajectory.npz", allow_pickle=False) as native:
            np.testing.assert_array_equal(times, native["time"])
            np.testing.assert_array_equal(positions, native["position"][:, 0])
            np.testing.assert_array_equal(rotations, native["rotation"][:, 0])
        rotations = rotations.reshape(-1, 3, 3)
        np.testing.assert_allclose(
            positions,
            pivot + np.einsum("tij,j->ti", rotations, centre - pivot),
            atol=1e-9,
            rtol=0,
        )
        np.testing.assert_allclose(
            np.einsum("tij,j->ti", rotations, axis),
            np.broadcast_to(axis, positions.shape),
            atol=1e-9,
            rtol=0,
        )
        if trial["axis"] == "Z":
            # Gravity is parallel to the fixed axis, so its torque about Z is zero.
            np.testing.assert_allclose(
                positions, np.broadcast_to(centre, positions.shape), atol=1e-9, rtol=0
            )
            np.testing.assert_allclose(
                rotations,
                np.broadcast_to(np.eye(3), rotations.shape),
                atol=1e-9,
                rtol=0,
            )
            rows.append(
                {
                    "axis": "Z",
                    "max_centre_drift_m": float(np.max(abs(positions - centre))),
                }
            )
        else:
            local_positions = (positions - pivot) @ rz
            theta = np.arctan2(local_positions[:, 0], -local_positions[:, 2])
            coarse, fine = reference(times, length, 16), reference(times, length, 32)
            gap = float(np.max(abs(coarse - fine)))
            error = float(np.max(abs(theta - fine)))
            assert gap < 1e-10, gap
            assert error < 3e-4, error
            rows.append(
                {
                    "axis": trial["axis"],
                    "length_m": length,
                    "mass_kg": mass,
                    "max_angle_error_rad": error,
                    "reference_refinement_gap_rad": gap,
                }
            )
    checked = {
        "status": "passed",
        "trials": rows,
        "reference": "Independent box integrals and scalar physical-pendulum RK4 (16/32 subdivisions); zero gravitational torque about vertical Z",
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
    (root / "verification.json").write_text(json.dumps(checked, indent=2) + "\n")
    return checked


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--installed-archive", type=Path)
    parser.add_argument("--lifecycle-root", type=Path)
    args = parser.parse_args()
    if args.installed_archive:
        report = json.loads((args.root / "report.json").read_text())
        with zipfile.ZipFile(args.installed_archive) as archive:
            manifest = json.loads(archive.read("Vinkulum/build-info.json"))
            assert manifest == report["installed_manifest"]
            for name, digest in manifest["sha256"].items():
                assert (
                    hashlib.sha256(archive.read("Vinkulum/" + name)).hexdigest()
                    == digest
                )
        if args.lifecycle_root:
            lifecycle = json.loads((args.lifecycle_root / "report.json").read_text())
            assert lifecycle["installed_manifest"] == manifest
            assert lifecycle["status"] == "passed"
            expected = {
                "failed_start": "failed",
                "cancel_starting": "cancelled",
                "cancel_running": "cancelled",
                "delete_source": "cancelled",
                "close_document": "cancelled",
                "late_cancel": "cancelled",
                "timeout": "failed",
                "changed_source": "failed",
                "success_after_cancellations": "completed",
                "close_window": "cancelled",
            }
            assert {case["name"]: case["signals"] for case in lifecycle["cases"]} == {
                k: [v] for k, v in expected.items()
            }
            assert all(
                case["child_reaped_before_notification"] for case in lifecycle["cases"]
            )
            assert (
                lifecycle["close_veto_kept_document"]
                and lifecycle["native_save_chosen"]
            )
            assert lifecycle["normal_freecad_close_after_worker_reaped"]
    print(json.dumps(verify(args.root), indent=2))
