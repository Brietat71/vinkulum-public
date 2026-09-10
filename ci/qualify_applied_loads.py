#!/usr/bin/env python3
"""Retain analytic comparisons, Hessian binding observations and separated timings."""

import argparse
import hashlib
import json
import os
import platform
import runpy
import statistics
import subprocess
import sys
import time
from dataclasses import replace
from importlib.metadata import version
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np

from vinkulum_studio.applied_loads import world_load_derivative
from vinkulum_studio.articulated_result import load_operators
from vinkulum_studio.document import Body, Law, Load, Project, joint_at, save_project
from vinkulum_studio.pinocchio_backend import ArticulatedModel, run_analysis

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = runpy.run_path(
    ROOT / "apps/studio/tests/fixtures/applied_loads/reference.py"
)


def write(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def samples(call, count):
    values = []
    for _ in range(count):
        start = time.perf_counter()
        call()
        values.append(time.perf_counter() - start)
    return {"seconds": values, "median_s": statistics.median(values)}


def differences(directory):
    rows = []
    for kind in ("spatial-force", "spatial-moment", "spatial-wrench", "mixed-wrench"):
        mixed = kind.startswith("mixed")
        force = (0, 0, 0) if kind.endswith("moment") else (0, 0, 5)
        moment = (0, 0, 0) if kind.endswith("force") else (0, 0, 13)
        if mixed:
            force, moment = (2, -3, 5), (7, 11, 13)
        project = REFERENCE["spatial_pair"](mixed=mixed, force=force, moment=moment)
        model = ArticulatedModel(project)
        state = {
            "q": [0.3, -0.4],
            "velocity": [0.2, -0.1],
            "acceleration": [0.7, -0.3],
            "time_s": 0.5,
        }
        result = run_analysis(project, directory / kind, state)
        arm = np.array(project.bodies[1].position) + project.loads[0].point
        *_, effort, expected = REFERENCE["rotation_reference"](
            state["q"], arm, np.array(force), np.array(moment), mixed=mixed
        )
        actual = np.array(result["external_effort_derivatives"]["q"])
        np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=3e-14)
        np.testing.assert_allclose(
            result["external_effort"], effort, rtol=1e-13, atol=3e-14
        )
        matrix = np.array(result["loaded_inverse_derivatives"]["q"])
        sweep = []
        for step in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
            numerical = {
                "external_effort": np.zeros((2, 2)),
                "inverse_effort": np.zeros((2, 2)),
            }
            for j in range(2):
                delta = np.eye(2)[j] * step
                plus = model.analyse(
                    **(state | {"q": (np.array(state["q"]) + delta).tolist()})
                )
                minus = model.analyse(
                    **(state | {"q": (np.array(state["q"]) - delta).tolist()})
                )
                for channel in numerical:
                    numerical[channel][:, j] = (
                        np.array(plus[channel]) - minus[channel]
                    ) / (2 * step)
            sweep.append(
                {
                    "step_coordinate_units": step,
                    "external_max_abs_error": float(
                        np.max(np.abs(numerical["external_effort"] - expected))
                    ),
                    "loaded_max_abs_error": float(
                        np.max(np.abs(numerical["inverse_effort"] - matrix))
                    ),
                }
            )
        assert min(row["external_max_abs_error"] for row in sweep) < 1e-8
        assert min(row["loaded_max_abs_error"] for row in sweep) < 1e-8
        condition = float(np.linalg.cond(expected))
        rows.append(
            {
                "case": kind,
                "state": state,
                "coordinate_units": [
                    row["coordinate_unit"] for row in result["coordinate_order"]
                ],
                "effort_units": [
                    row["effort_unit"] for row in result["coordinate_order"]
                ],
                "analytic_external_effort": effort.tolist(),
                "analytic_external_derivative_q": expected.tolist(),
                "analytic_max_abs_error": float(np.max(np.abs(actual - expected))),
                "scaled_derivative_norm_2": float(np.linalg.norm(expected, 2)),
                "scaled_derivative_condition_2": (
                    condition if np.isfinite(condition) else None
                ),
                "derivative_rank": int(np.linalg.matrix_rank(expected)),
                "step_sweep": sweep,
            }
        )
    return rows


def hessian_binding_probe():
    """Record the observed 4.1.0/EigenPy layout; never use it in production."""
    project = REFERENCE["spatial_pair"]()
    model = ArticulatedModel(project)
    q = np.array([0.3, -0.4])
    pin, data = model.pin, model.model.createData()
    pin.computeJointKinematicHessians(model.model, data, q)
    pin.updateFramePlacements(model.model, data)
    body = project.bodies[1]
    raw = np.asarray(
        pin.getFrameKinematicHessian(
            model.model, data, model.frame_ids[body.id], pin.LOCAL_WORLD_ALIGNED
        )
    )
    *_, H, dW, _, _ = REFERENCE["rotation_reference"](
        q, np.array(body.position), np.zeros(3), np.zeros(3)
    )
    expected = np.concatenate((H, dW), axis=0)
    unpacked = raw.ravel(order="C").reshape((6, 2, 2), order="F")
    raw_error = float(np.max(np.abs(raw - expected)))
    unpacked_error = float(np.max(np.abs(unpacked - expected)))
    # A future binding may expose correct ndarray strides; record either result.
    assert min(raw_error, unpacked_error) < 1e-13
    return {
        "axes": "[world motion component vx,vy,vz,wx,wy,wz; Jacobian column i; differentiated q_j] at body COM",
        "q": q.tolist(),
        "body_centre_local_m": list(body.position),
        "raw_shape": list(raw.shape),
        "raw_strides_bytes": list(raw.strides),
        "raw": raw.tolist(),
        "expected": expected.tolist(),
        "raw_max_abs_error": raw_error,
        "unpacked_max_abs_error": unpacked_error,
        "diagnostic_unpacking": "raw.ravel(order='C').reshape((6,nv,nv), order='F')",
        "production_uses_binding_hessian": False,
    }


def chain(count):
    uid = REFERENCE["fixture_id"]
    bodies = tuple(
        Body(
            uid(f"chain-{count}-{i}"),
            f"Body {i}",
            mass=1,
            dimensions=(0.03, 0.04, 0.05),
            position=(0.06 * i, 0.02 * (i % 3), 0.01 * (i % 2)),
        )
        for i in range(count)
    )
    project = Project(uid(f"chain-{count}"), bodies=bodies)
    axes = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    joints = tuple(
        replace(
            joint_at(
                project,
                "glissiere" if i % 4 == 3 else "pivot",
                None if i == 0 else bodies[i - 1].id,
                body.id,
                point=(0, 0, 0) if i == 0 else bodies[i - 1].position,
                axis=axes[i % 3],
            ),
            id=uid(f"chain-{count}-joint-{i}"),
        )
        for i, body in enumerate(bodies)
    )
    loads = tuple(
        Load(
            uid(f"chain-{count}-load-{i}"),
            f"Wrench {i}",
            bodies[i].id,
            point=(0.02, -0.03, 0.04),
            force=tuple(Law(values=(v,)) for v in (2, -3, 5)),
            moment=tuple(Law(values=(v,)) for v in (7, 11, 13)),
        )
        for i in sorted({count // 2, count - 1})
    )
    return replace(project, joints=joints, loads=loads)


def timings(directory, count, repeats):
    project = chain(count)
    model = ArticulatedModel(project)
    state = {
        "q": (0.1 * np.sin(np.arange(count))).tolist(),
        "velocity": (0.02 * np.cos(np.arange(count))).tolist(),
        "acceleration": (0.03 * np.sin(np.arange(count))).tolist(),
    }
    result = model.analyse(**state)
    encoded = json.dumps(result, allow_nan=False).encode()
    kinematics = model.kinematics(np.array(state["q"]))
    inputs = []
    for load in project.loads:
        pose, J = kinematics[load.body]
        point = J[:3] + np.cross(J[3:].T, pose.rotation @ load.point).T
        inputs.append(
            (
                point,
                J[3:],
                np.array([law.value(0) for law in load.force]),
                np.array([law.value(0) for law in load.moment]),
                model.ancestor_coordinates[load.body],
            )
        )

    def loads_only():
        return sum(
            (world_load_derivative(*args) for args in inputs), np.zeros((count, count))
        )

    np.testing.assert_allclose(
        loads_only(), result["external_effort_derivatives"]["q"], atol=1e-13
    )
    folder = directory / f"tree-{count}"
    folder.mkdir()
    save_project(folder / "project.vinkulum.json", project)
    write(folder / "state.json", state)
    # All timed subprocesses use the same separate worker interpreter.
    child = os.environ.copy()
    for name in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV"):
        child.pop(name, None)

    def child_call(args):
        subprocess.run(
            [sys.executable, *args],
            check=True,
            env=child,
            cwd=folder,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )

    plain = samples(lambda: child_call(["-c", "pass"]), repeats)
    imports = samples(
        lambda: child_call(
            ["-c", "import pinocchio; import vinkulum_studio.pinocchio_backend"]
        ),
        repeats,
    )
    operators = samples(lambda: model.analyse(**state), repeats)
    loads = samples(loads_only, repeats)
    encoding = samples(lambda: json.dumps(result, allow_nan=False).encode(), repeats)
    file_write = samples(
        lambda: (folder / "payload.json").write_bytes(encoded), repeats
    )
    start = time.perf_counter()
    child_call(
        [
            "-m",
            "vinkulum_studio.pinocchio_backend",
            str(folder / "project.vinkulum.json"),
            "--state",
            str(folder / "state.json"),
            "--output",
            str(folder / "operators"),
        ]
    )
    end_to_end = time.perf_counter() - start
    reopened = load_operators(folder / "operators")
    for channel in ("external_effort_derivatives", "loaded_inverse_derivatives"):
        self_contained = reopened.report[channel]
        for variable in self_contained:
            np.testing.assert_array_equal(
                self_contained[variable], result[channel][variable]
            )
    reading = samples(lambda: load_operators(folder / "operators"), repeats)
    return {
        "coordinates": count,
        "loads": len(project.loads),
        "operator_payload_bytes": len(encoded),
        "captured_report_bytes": (folder / "operators/result.json").stat().st_size,
        "python_startup": plain,
        "python_and_imports": imports,
        "model_construction": samples(lambda: ArticulatedModel(project), repeats),
        "all_operators_including_derivatives_and_finite_json_check": operators,
        "load_derivatives_with_precomputed_jacobians": loads,
        "json_encoding": encoding,
        "payload_file_write_without_fsync": file_write,
        "reader_validation": reading,
        "full_cli_once_s": end_to_end,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=7)
    args = parser.parse_args()
    if args.samples < 3:
        parser.error("At least three timing samples are required.")
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    report = {
        "format": "vinkulum-applied-load-qualification",
        "schema": 1,
        "environment": {
            name: version(name) for name in ("pin", "libpinocchio", "eigenpy", "numpy")
        },
        "python": sys.version,
        "platform": platform.platform(),
        "available_cpus": os.process_cpu_count(),
        "affinity": (
            sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else None
        ),
        "threads": {
            key: os.environ[key]
            for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        },
        "source_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [
                Path(__file__),
                ROOT / "apps/studio/tests/fixtures/applied_loads/reference.py",
                *[
                    ROOT / "apps/studio/vinkulum_studio" / name
                    for name in (
                        "applied_loads.py",
                        "articulated.py",
                        "articulated_result.py",
                        "pinocchio_backend.py",
                    )
                ],
            ]
        },
    }
    report["references"] = differences(directory)
    report["hessian_binding_probe"] = hessian_binding_probe()
    report["timings"] = [
        timings(directory, count, args.samples) for count in (2, 16, 32)
    ]
    write(directory / "qualification.json", report)
    print(
        json.dumps(
            {
                "directory": str(directory),
                "cases": len(report["references"]),
                "max_analytic_error": max(
                    row["analytic_max_abs_error"] for row in report["references"]
                ),
                "hessian_raw_error": report["hessian_binding_probe"][
                    "raw_max_abs_error"
                ],
                "timings_ms": [
                    {
                        "n": row["coordinates"],
                        "operators": row[
                            "all_operators_including_derivatives_and_finite_json_check"
                        ]["median_s"]
                        * 1000,
                        "load_derivatives": row[
                            "load_derivatives_with_precomputed_jacobians"
                        ]["median_s"]
                        * 1000,
                    }
                    for row in report["timings"]
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
