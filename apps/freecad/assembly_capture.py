"""Capture native FreeCAD Assembly revolute mechanisms for an external engine.

SPDX-License-Identifier: Apache-2.0
Runs only in FreeCAD. Calls its installed Assembly placement APIs.
"""

import json
import math
import uuid
from pathlib import Path

import FreeCAD as App
import Part
import UtilsAssembly

from .bridge import matrix_values, request_bytes, sha


def world_shape(component):
    """Include parent Assembly placement; App::Link has no getGlobalPlacement."""
    world = UtilsAssembly.getGlobalPlacement((component, ["", ""]))
    shape = component.Shape.copy()
    shape.Placement = world.multiply(component.Placement.inverse()).multiply(
        shape.Placement
    )
    if not shape.isValid() or len(shape.Solids) != 1 or shape.Volume <= 0:
        raise ValueError(f"{component.Label}: one valid positive-volume solid required")
    return shape.Solids[0]


def _frame(placement):
    return {
        "position_m": [v * 0.001 for v in placement.Base],
        "rotation": matrix_values(placement.toMatrix()),
    }


def snapshot(assembly):
    """Read a flat native Assembly without changing its feature graph or poses."""
    if assembly.TypeId != "Assembly::AssemblyObject":
        raise ValueError("Select a native Assembly object")
    groups = [o for o in assembly.Group if o.TypeId == "Assembly::JointGroup"]
    if len(groups) != 1:
        raise ValueError("The assembly requires one native joint group")
    components = UtilsAssembly.getMovablePartsWithin(assembly, partsAsSolid=True)
    if not 2 <= len(components) <= 17:
        raise ValueError(
            "Capture supports one grounded component and up to 16 moving solids"
        )
    if any(
        not hasattr(o, "Shape") or o.TypeId == "Assembly::AssemblyObject"
        for o in components
    ):
        raise ValueError(
            "Nested assemblies and compound parts require their own conversion"
        )
    shapes = {o.Name: world_shape(o) for o in components}
    if len(shapes) != len(components):
        raise ValueError("Duplicate component identity")
    ground = [j.ObjectToGround for j in groups[0].Group if hasattr(j, "ObjectToGround")]
    if len(ground) != 1 or ground[0] not in components:
        raise ValueError("One component must be grounded in the native Assembly")
    joints = []
    for joint in groups[0].Group:
        if hasattr(joint, "ObjectToGround"):
            continue
        if getattr(joint, "Suppressed", False):
            continue
        if str(getattr(joint, "JointType", "")) != "Revolute":
            raise ValueError(
                f"{joint.Label}: only native revolute joints are currently mapped"
            )
        if any(
            getattr(joint, key, False)
            for key in (
                "EnableAngleMin",
                "EnableAngleMax",
                "EnableLengthMin",
                "EnableLengthMax",
            )
        ):
            raise ValueError(f"{joint.Label}: joint limits are not yet mapped")
        refs = (joint.Reference1, joint.Reference2)
        if any(
            not UtilsAssembly.isRefValid(ref, 1) or ref[0] not in components
            for ref in refs
        ):
            raise ValueError(f"{joint.Label}: missing or external component reference")
        frames = [
            _frame(UtilsAssembly.getJcsGlobalPlc(placement, ref))
            for placement, ref in zip((joint.Placement1, joint.Placement2), refs)
        ]
        points = [frame["position_m"] for frame in frames]
        axes = [[frame["rotation"][i] for i in (2, 5, 8)] for frame in frames]
        if (
            math.dist(*points) > 1e-8
            or abs(abs(sum(a * b for a, b in zip(*axes))) - 1) > 1e-8
        ):
            raise ValueError(
                f"{joint.Label}: revolute connector origins and axes must agree"
            )
        names = [ref[0].Name for ref in refs]
        if names[0] == names[1]:
            raise ValueError(f"{joint.Label}: a joint needs two distinct components")
        joints.append(
            {
                "name": joint.Name,
                "label": joint.Label,
                "components": names,
                "references": [list(ref[1]) for ref in refs],
                "frames": frames,
                "axis_world": axes[0],
                "pivot_m": points[0],
            }
        )
    if not joints:
        raise ValueError("The assembly has no active revolute joints")
    state = {
        "assembly_name": assembly.Name,
        "grounded_name": ground[0].Name,
        "components": [
            {
                "name": obj.Name,
                "label": obj.Label,
                "geometry_sha256": sha(
                    shapes[obj.Name].exportBrepToString().encode("ascii")
                ),
            }
            for obj in sorted(components, key=lambda o: o.Name)
        ],
        "joints": sorted(joints, key=lambda j: j["name"]),
    }
    return state, shapes


def capture(assembly, directory, *, density=2700, duration=0.5, step=0.005, threads=2):
    if not math.isfinite(density) or density <= 0:
        raise ValueError("A finite positive density is required")
    if not (0 < step <= duration <= 10 and math.ceil(duration / step) <= 2000):
        raise ValueError("Capture is limited to 10 seconds and 2001 native samples")
    if type(threads) is not int or not 1 <= threads <= 64:
        raise ValueError("Select between 1 and 64 engine threads")
    state, shapes = snapshot(assembly)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    project_id = str(uuid.uuid4())
    bodies = []
    ids = {state["grounded_name"]: None}
    for component in state["components"]:
        name = component["name"]
        if name == state["grounded_name"]:
            continue
        shape = shapes[name]
        body_id = str(uuid.uuid4())
        ids[name] = body_id
        step_file = body_id + ".step"
        shape.exportStep(str(directory / step_file))
        bodies.append(
            component
            | {
                "body_id": body_id,
                "step_file": step_file,
                "step_sha256": sha((directory / step_file).read_bytes()),
                "properties_si": {
                    "volume_m3": shape.Volume * 1e-9,
                    "mass_kg": density * shape.Volume * 1e-9,
                    "centre_m": [v * 0.001 for v in shape.CenterOfMass],
                    "inertia_kg_m2": [
                        v * density * 1e-15
                        for v in matrix_values(shape.MatrixOfInertia)
                    ],
                },
            }
        )
    request = {
        "format": "vinkulum-freecad-assembly-1",
        "run_id": str(uuid.uuid4()),
        "project_id": project_id,
        "label": assembly.Label,
        "source_state": state,
        "source_state_sha256": sha(request_bytes(state)),
        "bodies": bodies,
        "joints": [
            joint
            | {
                "body_a": ids[joint["components"][0]],
                "body_b": ids[joint["components"][1]],
            }
            for joint in state["joints"]
        ],
        "density_kg_m3": density,
        "duration_s": duration,
        "step_s": step,
        "threads": threads,
        "freecad_version": App.Version(),
        "freecad_occt": Part.OCC_VERSION,
    }
    (directory / "request.json").write_bytes(request_bytes(request))
    return request


def admit(assembly, directory, request):
    state, _ = snapshot(assembly)
    if sha(request_bytes(state)) != request["source_state_sha256"]:
        raise ValueError(
            "Assembly geometry, grounding or joints changed during calculation"
        )
    directory = Path(directory)
    captured = request_bytes(request)
    if (directory / "request.json").read_bytes() != captured:
        raise ValueError("Captured assembly request changed")
    path = directory / "result/playback.json"
    if path.stat().st_size > 32_000_000:
        raise ValueError("Assembly playback exceeds its input budget")
    result = json.loads(path.read_text())
    ids = [b["body_id"] for b in request["bodies"]]
    if (
        result.get("format") != "vinkulum-freecad-assembly-playback-1"
        or result.get("run_id") != request["run_id"]
        or result.get("request_sha256") != sha(captured)
        or result.get("body_ids") != ids
    ):
        raise ValueError("Playback does not identify the captured assembly")
    times, positions, rotations = (
        result[k] for k in ("time_s", "position_m", "rotation")
    )
    if not (2 <= len(times) <= 2001 and len(times) == len(positions) == len(rotations)):
        raise ValueError("Invalid native sample count")
    if times[0] != 0 or abs(times[-1] - request["duration_s"]) > 1e-9:
        raise ValueError("Incomplete native trajectory")
    for n, (t, ps, rs) in enumerate(zip(times, positions, rotations)):
        if (
            type(t) not in (int, float)
            or not math.isfinite(t)
            or (n and times[n - 1] >= t)
            or len(ps) != len(ids)
            or len(rs) != len(ids)
        ):
            raise ValueError("Invalid times or body count")
        for i, (p, r) in enumerate(zip(ps, rs)):
            if (
                len(p) != 3
                or len(r) != 9
                or not all(
                    type(v) in (int, float) and math.isfinite(v) for v in (*p, *r)
                )
            ):
                raise ValueError("Invalid native pose")
            for a in range(3):
                for b in range(3):
                    if (
                        abs(
                            sum(r[3 * k + a] * r[3 * k + b] for k in range(3))
                            - (a == b)
                        )
                        > 1e-8
                    ):
                        raise ValueError("Native rotation is not orthogonal")
            det = sum(
                r[k]
                * (
                    r[3 + (k + 1) % 3] * r[6 + (k + 2) % 3]
                    - r[3 + (k + 2) % 3] * r[6 + (k + 1) % 3]
                )
                for k in range(3)
            )
            if abs(det - 1) > 1e-8:
                raise ValueError("Native rotation is not proper")
            if not n and (
                math.dist(p, request["bodies"][i]["properties_si"]["centre_m"]) > 1e-9
                or any(
                    abs(a - b) > 1e-9 for a, b in zip(r, (1, 0, 0, 0, 1, 0, 0, 0, 1))
                )
            ):
                raise ValueError("Initial pose differs from the captured body")
    return result
