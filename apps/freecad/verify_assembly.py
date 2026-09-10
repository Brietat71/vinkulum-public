"""Independent finite-section double-pendulum checks for native Assembly captures.

SPDX-License-Identifier: Apache-2.0
Requires NumPy only; imports no CAD or dynamics engine.
"""

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


def reference(times, subdivisions):
    lengths = np.array((0.8, 0.6))
    masses = 2700 * 0.02 * 0.03 * lengths
    centres = lengths / 2
    inertia = masses * (lengths**2 + 0.02**2) / 12
    a = inertia[0] + masses[0] * centres[0] ** 2 + masses[1] * lengths[0] ** 2
    d = inertia[1] + masses[1] * centres[1] ** 2
    b = masses[1] * lengths[0] * centres[1]
    gravity = 9.80665 * np.array(
        (masses[0] * centres[0] + masses[1] * lengths[0], masses[1] * centres[1])
    )

    def rhs(y):
        alpha, beta, wa, wb = y
        delta = alpha - beta
        off = b * math.cos(delta)
        fa = -b * math.sin(delta) * wb**2 - gravity[0] * math.sin(alpha)
        fb = b * math.sin(delta) * wa**2 - gravity[1] * math.sin(beta)
        determinant = a * d - off**2
        return np.array(
            (
                wa,
                wb,
                (d * fa - off * fb) / determinant,
                (a * fb - off * fa) / determinant,
            )
        )

    y = np.array((math.radians(20), math.radians(-25), 0.0, 0.0))
    rows = [y.copy()]
    for start, end in zip(times, times[1:]):
        h = (end - start) / subdivisions
        for _ in range(subdivisions):
            k1 = rhs(y)
            k2 = rhs(y + h * k1 / 2)
            k3 = rhs(y + h * k2 / 2)
            k4 = rhs(y + h * k3)
            y += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        rows.append(y.copy())
    return np.array(rows)[:, :2]


def verify(root):
    report = json.loads((root / "report.json").read_text())
    assert report["status"] == "passed", report
    assert {(t["yaw_deg"], t["step_s"]) for t in report["trials"]} == {
        (yaw, h) for yaw in (0, 17) for h in (0.005, 0.0025)
    }
    expected_rejections = {
        "displaced_joint_connector",
        "unsupported_joint",
        "unmapped_joint_limit",
        "missing_component_reference",
        "moved_parent_assembly",
        "changed_grounding",
        "improper_rotation",
        "swapped_body_identity",
        "boolean_native_time",
        "replaced_request_and_matching_result_hash",
        "changed_parametric_shape",
    }
    assert set(report["rejections"]) == expected_rejections
    assert all(
        isinstance(error, str) and error for error in report["rejections"].values()
    )
    rows = []
    errors = {}
    series = {}
    for trial in report["trials"]:
        directory = root / trial["directory"]
        raw = (directory / "request.json").read_bytes()
        request = json.loads(raw)
        result = json.loads((directory / "result/playback.json").read_text())
        metadata = json.loads((directory / "result/result.json").read_text())
        assert hashlib.sha256(raw).hexdigest() == result["request_sha256"]
        assert request["run_id"] == result["run_id"] == metadata["run_id"]
        assert (
            hashlib.sha256(
                (directory / "result/trajectory.npz").read_bytes()
            ).hexdigest()
            == metadata["archive_sha256"]
        )
        names = [body["name"] for body in request["bodies"]]
        assert set(names) == {"Upper", "Lower"}
        assert [body["body_id"] for body in request["bodies"]] == result["body_ids"]
        yaw = math.radians(trial["yaw_deg"])
        rz = np.array(
            (
                (math.cos(yaw), -math.sin(yaw), 0),
                (math.sin(yaw), math.cos(yaw), 0),
                (0, 0, 1),
            )
        )
        origin = np.array(trial["translation_mm"]) * 0.001
        shoulder = np.array(
            (0.8 * math.sin(math.radians(20)), 0, -0.8 * math.cos(math.radians(20)))
        )
        for i, body in enumerate(request["bodies"]):
            assert (
                hashlib.sha256((directory / body["step_file"]).read_bytes()).hexdigest()
                == body["step_sha256"]
            )
            length, angle = (0.8, 20) if body["name"] == "Upper" else (0.6, -25)
            mass = 2700 * 0.02 * 0.03 * length
            q = math.radians(180 - angle)
            orientation = rz @ np.array(
                (
                    (math.cos(q), 0, math.sin(q)),
                    (0, 1, 0),
                    (-math.sin(q), 0, math.cos(q)),
                )
            )
            centre = origin + orientation @ np.array((0, 0, length / 2))
            if body["name"] == "Lower":
                centre += rz @ shoulder
            inertia = (
                orientation
                @ np.diag(
                    mass
                    * np.array(
                        (0.03**2 + length**2, 0.02**2 + length**2, 0.02**2 + 0.03**2)
                    )
                    / 12
                )
                @ orientation.T
            )
            for actual in (body["properties_si"], result["properties_si"][i]):
                for key, expected, budget in (
                    ("volume_m3", mass / 2700, 1e-8 * mass / 2700),
                    ("mass_kg", mass, 1e-8 * mass),
                    ("centre_m", centre, 1e-8 * length),
                    ("inertia_kg_m2", inertia.ravel(), 1e-8 * mass * length**2),
                ):
                    np.testing.assert_allclose(
                        actual[key], expected, atol=budget, rtol=0
                    )
        times, positions, rotations = (
            np.array(result[key]) for key in ("time_s", "position_m", "rotation")
        )
        with np.load(directory / "result/trajectory.npz", allow_pickle=False) as native:
            np.testing.assert_array_equal(times, native["time"])
            np.testing.assert_array_equal(positions, native["position"])
            np.testing.assert_array_equal(rotations, native["rotation"])
        local = (positions - origin) @ rz
        upper = local[:, names.index("Upper")]
        lower = local[:, names.index("Lower")] - 2 * upper
        angles = np.column_stack(
            (
                np.arctan2(upper[:, 0], -upper[:, 2]),
                np.arctan2(lower[:, 0], -lower[:, 2]),
            )
        )
        for p, length in ((upper, 0.4), (lower, 0.3)):
            np.testing.assert_allclose(
                np.linalg.norm(p, axis=1), length, atol=1e-9, rtol=0
            )
            np.testing.assert_allclose(p[:, 1], 0, atol=1e-9, rtol=0)
        coarse, fine = reference(times, 16), reference(times, 32)
        gap = float(np.max(abs(coarse - fine)))
        error = float(np.max(abs(angles - fine)))
        assert gap < 1e-10, gap
        assert error < 2e-3, error
        if trial["step_s"] == 0.0025:
            assert error < 5e-4, error
        assert trial["heartbeat_ticks"] > 2
        assert len(trial["replayed_samples"]) == 41
        replay = max(row["centre_error_m"] for row in trial["replayed_samples"])
        assert replay < 1e-9, replay
        errors[trial["yaw_deg"], trial["step_s"]] = error
        series[trial["yaw_deg"], trial["step_s"]] = (
            positions,
            rotations.reshape(-1, 2, 3, 3),
            rz,
            origin,
        )
        rows.append(
            {
                "yaw_deg": trial["yaw_deg"],
                "step_s": trial["step_s"],
                "max_angle_error_rad": error,
                "reference_refinement_gap_rad": gap,
                "replayed_centre_error_m": replay,
            }
        )
    ratios = {str(yaw): errors[yaw, 0.0025] / errors[yaw, 0.005] for yaw in (0, 17)}
    assert all(value < 0.4 for value in ratios.values()), ratios
    transport = []
    for h in (0.005, 0.0025):
        pos0, rot0, _, _ = series[0, h]
        pos, rot, rz, origin = series[17, h]
        position_error = float(np.max(abs(pos - (pos0 @ rz.T + origin))))
        rotation_error = float(
            np.max(abs(rot - np.einsum("ij,tbjk,kl->tbil", rz, rot0, rz.T)))
        )
        assert position_error < 1e-8 and rotation_error < 1e-8, (
            position_error,
            rotation_error,
        )
        transport.append(
            {
                "step_s": h,
                "position_error_m": position_error,
                "rotation_error": rotation_error,
            }
        )
    checked = {
        "status": "passed",
        "trials": rows,
        "rejections_checked": sorted(expected_rejections),
        "fine_over_coarse_error_ratio": ratios,
        "world_frame_transport": transport,
        "reference": "Independent 2x2 Lagrange equations for two uniform finite-section rods, RK4 with 16/32 subdivisions",
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
    print(json.dumps(verify(parser.parse_args().root), indent=2))
