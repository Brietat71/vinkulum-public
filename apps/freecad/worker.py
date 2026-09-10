"""Experimental FreeCAD bridge: run in Vinkulum's separate Python environment.

SPDX-License-Identifier: Apache-2.0
"""

import hashlib
import sys
import uuid
from dataclasses import replace
from pathlib import Path

import numpy as np


def import_body(request, step):
    """Reimport a captured solid and check all SI properties against FreeCAD."""
    from vinkulum_studio.cad import execute
    from vinkulum_studio.document import Body

    body = Body.from_dict(
        execute(
            {
                "operation": "import_step",
                "path": str(step),
                "density": request["density_kg_m3"],
                "name": request["label"],
            }
        )["body"]
    )
    body = replace(body, id=request["body_id"])
    expected = request["properties_si"]
    epsilon, length = 1e-8, float(np.linalg.norm(body.dimensions))
    for key in ("volume_m3", "mass_kg"):
        value = expected[key]
        if type(value) not in (float, int) or not np.isfinite(value) or value <= 0:
            raise ValueError("Invalid captured physical property")
    for actual, target, tolerance in (
        (body.cad.volume_m3, expected["volume_m3"], epsilon * expected["volume_m3"]),
        (body.mass, expected["mass_kg"], epsilon * expected["mass_kg"]),
        (body.position, expected["centre_m"], epsilon * length),
        (body.inertia(), expected["inertia_kg_m2"], epsilon * body.mass * length**2),
    ):
        np.testing.assert_allclose(actual, target, atol=tolerance, rtol=0)
    return body


def solve(directory):
    from vinkulum_studio.document import Project, joint_at, save_project
    from vinkulum_studio.engine_threads import configure_occt_threads
    from vinkulum_studio.execution import ExecutionPlan
    from vinkulum_studio.mechanism import MechanicalResult, simulate_project
    from vinkulum_studio.model import read_json, write_json

    path = directory / "request.json"
    request = read_json(path, 100_000)
    if request["format"] != "vinkulum-freecad-prototype-1":
        raise ValueError("Unsupported bridge request")
    if str(uuid.UUID(request["run_id"])) != request["run_id"]:
        raise ValueError("Invalid request identity")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    step = directory / "part.step"
    if step.stat().st_size > 8_000_000:
        raise ValueError("STEP input exceeds the existing Studio import budget")
    if hashlib.sha256(step.read_bytes()).hexdigest() != request["step_sha256"]:
        raise ValueError("Captured STEP changed")
    output = directory / "result"
    output.mkdir(exist_ok=False)
    threads = request.get("threads", 2)
    if type(threads) is not int or not 1 <= threads <= 64:
        raise ValueError("Invalid engine thread allocation")
    plan = ExecutionPlan(threads, threads)
    cad_runtime = configure_occt_threads(plan.threads)
    body = import_body(request, step)
    project = Project(
        request["project_id"],
        request["label"],
        bodies=(body,),
        duration=request["duration_s"],
        step=request["step_s"],
    )
    joint = joint_at(
        project,
        "pivot",
        None,
        body.id,
        point=request["pivot_m"],
        axis=request["axis_world"],
    )
    project = project.replace_object(joint)
    save_project(directory / "project.vinkulum.json", project)
    result = simulate_project(project, request["run_id"], output, execution=plan)
    write_json(output / "result.json", result)
    admitted = MechanicalResult.read(result, request["run_id"], project, output)
    # FreeCAD needs no Vinkulum, numpy or OCCT-8 imports to replay these samples.
    playback = {
        "format": "vinkulum-freecad-playback-1",
        "run_id": request["run_id"],
        "request_sha256": digest,
        "body_id": body.id,
        "properties_si": {
            "volume_m3": body.cad.volume_m3,
            "mass_kg": body.mass,
            "centre_m": body.position,
            "inertia_kg_m2": body.inertia(),
        },
        "time_s": admitted.time.tolist(),
        "position_m": admitted.position[:, 0].tolist(),
        "rotation": admitted.rotation[:, 0].tolist(),
        "manifest": result["manifest"],
        "cad_runtime": cad_runtime,
    }
    write_json(output / "playback.json", playback)


if __name__ == "__main__":
    try:
        solve(Path(sys.argv[1]).resolve())
    except (
        Exception
    ) as error:  # noqa: BLE001 - report native errors at the process/UI boundary
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1)
