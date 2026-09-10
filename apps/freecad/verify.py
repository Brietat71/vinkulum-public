"""Independent rectangular physical-pendulum reference for the FreeCAD probe.

SPDX-License-Identifier: Apache-2.0
Needs numpy; imports neither FreeCAD, Studio nor a geometry/physics engine.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def reference(times, length, subdivisions):
    # Rectangle 20 x 30 mm, revolute axis world Y at the end of its length.
    coefficient = 9.80665 * (length / 2) / (length**2 / 3 + 0.020**2 / 12)
    angle, speed = math.radians(20), 0.0
    result = [angle]

    def derivative(a, v):
        return v, -coefficient * math.sin(a)

    for start, end in zip(times, times[1:]):
        h = (end - start) / subdivisions
        for _ in range(subdivisions):
            a1, v1 = derivative(angle, speed)
            a2, v2 = derivative(angle + h * a1 / 2, speed + h * v1 / 2)
            a3, v3 = derivative(angle + h * a2 / 2, speed + h * v2 / 2)
            a4, v4 = derivative(angle + h * a3, speed + h * v3)
            angle += h * (a1 + 2 * a2 + 2 * a3 + a4) / 6
            speed += h * (v1 + 2 * v2 + 2 * v3 + v4) / 6
        result.append(angle)
    return np.array(result)


def verify(root, plot=False):
    report = json.loads((root / "report.json").read_text())
    assert report["status"] == "completed", report
    assert report["stale_geometry_rejected"] and report["corrupt_rotation_rejected"]
    assert report["moved_parent_rejected"]
    assert report["changed_request_rejected"]
    assert {(r["length_mm"], r["step_s"]) for r in report["trials"]} == {
        (l, h) for l in (800, 1200) for h in (0.005, 0.0025)
    }
    errors, rows, series, identities = {}, [], {}, set()
    for trial in report["trials"]:
        directory = root / trial["directory"]
        request_path = directory / "request.json"
        request = json.loads(request_path.read_text())
        data = json.loads((directory / "result/playback.json").read_text())
        metadata = json.loads((directory / "result/result.json").read_text())
        assert (
            hashlib.sha256(request_path.read_bytes()).hexdigest()
            == data["request_sha256"]
        )
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
        assert data["run_id"] == request["run_id"] == metadata["run_id"]
        assert data["run_id"] not in identities
        identities.add(data["run_id"])
        length, step = trial["length_mm"] * 1e-3, trial["step_s"]
        mass = 0.020 * 0.030 * length * 7800
        angle = math.radians(160)
        c, s = math.cos(angle), math.sin(angle)
        q = np.array(((c, 0, s), (0, 1, 0), (-s, 0, c)))
        centre = q @ (0, 0, length / 2)
        inertia = (
            q
            @ np.diag(
                mass
                * np.array(
                    (0.030**2 + length**2, 0.020**2 + length**2, 0.020**2 + 0.030**2)
                )
                / 12
            )
            @ q.T
        )
        for source in (request["properties_si"], data["properties_si"]):
            np.testing.assert_allclose(
                source["volume_m3"], mass / 7800, atol=1e-8 * mass / 7800, rtol=0
            )
            np.testing.assert_allclose(
                source["mass_kg"], mass, atol=1e-8 * mass, rtol=0
            )
            np.testing.assert_allclose(
                source["centre_m"], centre, atol=1e-8 * length, rtol=0
            )
            np.testing.assert_allclose(
                source["inertia_kg_m2"],
                inertia.ravel(),
                atol=1e-8 * mass * length**2,
                rtol=0,
            )
        times, position, rotation = (
            np.array(data[key]) for key in ("time_s", "position_m", "rotation")
        )
        with np.load(directory / "result/trajectory.npz", allow_pickle=False) as native:
            np.testing.assert_array_equal(times, native["time"])
            np.testing.assert_array_equal(position, native["position"][:, 0])
            np.testing.assert_array_equal(rotation, native["rotation"][:, 0])
        np.testing.assert_allclose(position[:, 1], 0, atol=1e-10, rtol=0)
        np.testing.assert_allclose(
            position,
            np.einsum("tij,j->ti", rotation.reshape(-1, 3, 3), centre),
            atol=1e-9,
            rtol=0,
        )
        theta = np.arctan2(position[:, 0], -position[:, 2])
        expected = reference(times, length, 16)
        refined = reference(times, length, 32)
        refinement = float(np.max(abs(refined - expected)))
        assert refinement < 1e-10, refinement
        error = float(np.max(abs(theta - refined)))
        errors[length, step] = error
        assert error < 3e-4, error
        assert trial["heartbeat_ticks_during_calculation"] > 2
        assert len(trial["replayed_samples"]) >= 40
        replay_error = max(row["centre_error_m"] for row in trial["replayed_samples"])
        assert replay_error < 1e-9
        rows.append(
            {
                "length_m": length,
                "step_s": step,
                "mass_kg": mass,
                "max_angle_error_rad": error,
                "reference_refinement_gap_rad": refinement,
                "max_replayed_centre_error_m": replay_error,
            }
        )
        series[length, step] = (times, theta, refined)
    ratios = {
        str(length): errors[length, 0.0025] / errors[length, 0.005]
        for length in (0.8, 1.2)
    }
    assert all(ratio < 0.4 for ratio in ratios.values()), ratios
    verified = {
        "status": "passed",
        "reference": "theta'' = -g(L/2)/(L^2/3 + width^2/12) sin(theta); theta(0)=20deg, theta'(0)=0",
        "trials": rows,
        "fine_over_coarse_error_ratio": ratios,
        "scope": "Two uniform rectangular rigid pendula, explicit revolute joint; no general trajectory error bound",
    }
    (root / "verification.json").write_text(json.dumps(verified, indent=2) + "\n")
    if plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots(2, 2, figsize=(10, 6), layout="constrained")
        for column, length in enumerate((0.8, 1.2)):
            for h, label in ((0.005, "Vinkulum h=5 ms"), (0.0025, "Vinkulum h=2.5 ms")):
                times, theta, expected = series[length, h]
                axes[0, column].plot(times, np.degrees(theta), label=label)
                axes[1, column].plot(times, (theta - expected) * 1e3, label=label)
            axes[0, column].plot(
                times, np.degrees(expected), "k--", linewidth=1, label="Independent ODE"
            )
            axes[0, column].set_title(f"FreeCAD parametric length: {length:g} m")
            axes[0, column].set_ylabel("Angle from downward vertical (degrees)")
            axes[1, column].set_ylabel("Native − reference (mrad)")
            axes[1, column].set_xlabel("Time (s)")
            axes[0, column].legend(fontsize=8)
            for row in range(2):
                axes[row, column].grid(alpha=0.25)
        plot_path = root / "pendulum-reference.svg"
        figure.savefig(plot_path, metadata={"Date": None})
        plot_path.write_text(
            "\n".join(line.rstrip() for line in plot_path.read_text().splitlines())
            + "\n"
        )
        plt.close(figure)
    print(json.dumps(verified, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    verify(args.root.resolve(), args.plot)
