"""Calculate a native FreeCAD Assembly capture in the external Vinkulum runtime.

SPDX-License-Identifier: Apache-2.0
"""

import hashlib
import sys
import uuid
from dataclasses import replace
from pathlib import Path

from worker import import_body


def solve(directory):
    from vinkulum_studio.document import Project, joint_at, save_project
    from vinkulum_studio.engine_threads import configure_occt_threads
    from vinkulum_studio.execution import ExecutionPlan
    from vinkulum_studio.mechanism import MechanicalResult, simulate_project
    from vinkulum_studio.model import read_json, write_json

    path = directory / "request.json"
    request = read_json(path, 2_000_000)
    if request["format"] != "vinkulum-freecad-assembly-1":
        raise ValueError("Unsupported native Assembly capture")
    for key in ("run_id", "project_id"):
        if str(uuid.UUID(request[key])) != request[key]:
            raise ValueError("Invalid capture identity")
    captures = request["bodies"]
    ids = [body["body_id"] for body in captures]
    if not 1 <= len(captures) <= 16 or len(set(ids)) != len(ids):
        raise ValueError("Invalid body identities or capture size")
    threads = request["threads"]
    if type(threads) is not int or not 1 <= threads <= 64:
        raise ValueError("Invalid thread allocation")
    plan = ExecutionPlan(threads, threads)
    cad_runtime = configure_occt_threads(threads)
    bodies = []
    for capture in captures:
        name = capture["step_file"]
        if Path(name).name != name:
            raise ValueError("STEP must belong to its capture directory")
        step = directory / name
        if step.stat().st_size > 8_000_000:
            raise ValueError("STEP exceeds the importer budget")
        if hashlib.sha256(step.read_bytes()).hexdigest() != capture["step_sha256"]:
            raise ValueError("Captured STEP changed")
        bodies.append(
            import_body(capture | {"density_kg_m3": request["density_kg_m3"]}, step)
        )
    project = Project(
        request["project_id"],
        request["label"],
        bodies=tuple(bodies),
        duration=request["duration_s"],
        step=request["step_s"],
    )
    for captured in request["joints"]:
        for name in ("body_a", "body_b"):
            if captured[name] is not None and captured[name] not in ids:
                raise ValueError("A joint references an uncaptured body")
        joint = joint_at(
            project,
            "pivot",
            captured["body_a"],
            captured["body_b"],
            point=captured["pivot_m"],
            axis=captured["axis_world"],
        )
        joint = replace(
            joint,
            id=str(uuid.uuid5(uuid.UUID(request["project_id"]), captured["name"])),
            name=captured["label"],
        )
        project = project.replace_object(joint)
    output = directory / "result"
    output.mkdir(exist_ok=False)
    save_project(directory / "project.vinkulum.json", project)
    result = simulate_project(project, request["run_id"], output, execution=plan)
    write_json(output / "result.json", result)
    admitted = MechanicalResult.read(result, request["run_id"], project, output)
    playback = {
        "format": "vinkulum-freecad-assembly-playback-1",
        "run_id": request["run_id"],
        "request_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "body_ids": ids,
        "time_s": admitted.time.tolist(),
        "position_m": admitted.position.tolist(),
        "rotation": admitted.rotation.tolist(),
        "properties_si": [
            {
                "volume_m3": body.cad.volume_m3,
                "mass_kg": body.mass,
                "centre_m": body.position,
                "inertia_kg_m2": body.inertia(),
            }
            for body in bodies
        ],
        "manifest": result["manifest"],
        "cad_runtime": cad_runtime,
    }
    write_json(output / "playback.json", playback)


if __name__ == "__main__":
    try:
        solve(Path(sys.argv[1]).resolve())
    except Exception as error:
        print(f"{type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1)
