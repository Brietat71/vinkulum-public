"""External OCCT 8 / Gmsh / CalculiX experiment for captured FreeCAD faces.

SPDX-License-Identifier: Apache-2.0
"""

import argparse
import hashlib
import math
from dataclasses import replace
from pathlib import Path

import build123d as bd
import numpy as np
from vinkulum_studio.cad import execute
from vinkulum_studio.calculix import load_static_result, run_static
from vinkulum_studio.document import Body
from vinkulum_studio.engine_threads import configure_occt_threads
from vinkulum_studio.mesh_binding import MeshBinding, MeshCondition, surface_integrals
from vinkulum_studio.meshing import HxtMeshRequest, load_mesh, run_mesh
from vinkulum_studio.model import read_json, write_json


def checked_step(path, expected):
    if path.stat().st_size > 8_000_000:
        raise ValueError("Captured STEP exceeds its input budget.")
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError("Captured STEP changed.")


def import_body(root, request):
    checked_step(root / "part.step", request["step_sha256"])
    body = Body.from_dict(
        execute(
            {
                "operation": "import_step",
                "path": str(root / "part.step"),
                "density": request["density_kg_m3"],
                "name": request["label"],
            }
        )["body"]
    )
    body = replace(body, id=request["body_id"])
    expected = request["properties_si"]
    length = float(np.linalg.norm(body.dimensions))
    for actual, target, scale in (
        (body.cad.volume_m3, expected["volume_m3"], body.cad.volume_m3),
        (body.mass, expected["mass_kg"], body.mass),
        (body.position, expected["centre_m"], length),
        (body.inertia(), expected["inertia_kg_m2"], body.mass * length**2),
    ):
        if not np.isfinite(target).all():
            raise ValueError("Invalid captured SI property.")
        np.testing.assert_allclose(actual, target, atol=1e-8 * scale, rtol=0)
    return body


def import_faces(root, request, length_m):
    captures = request["faces"]
    if not isinstance(captures, list) or not 1 <= len(captures) <= 32:
        raise ValueError("Capture between 1 and 32 boundary faces.")
    faces = []
    for index, capture in enumerate(captures):
        path = root / f"face-{index:03}.step"
        checked_step(path, capture["step_sha256"])
        imported = bd.import_step(path).faces()
        if len(imported) != 1:
            raise ValueError("Each boundary STEP must contain exactly one face.")
        face = imported[0]
        if not math.isfinite(capture["area_m2"]) or capture["area_m2"] <= 0:
            raise ValueError("Invalid captured face area.")
        np.testing.assert_allclose(
            face.area * 1e-6, capture["area_m2"], rtol=1e-8, atol=0
        )
        np.testing.assert_allclose(
            np.array(tuple(face.center(bd.CenterOf.MASS))) * 1e-3,
            capture["centre_m"],
            atol=1e-8 * length_m,
            rtol=0,
        )
        faces.append(face)
    return faces


def match_faces(solid, faces, captures):
    """Use trimmed geometry membership and area coverage, never CAD/Gmsh indices."""
    points = np.array(solid.mesh.nodes)
    length = float(np.linalg.norm(solid.request.body.dimensions))
    tolerance_mm = max(1e-6, length * 1e-5)
    boundary_nodes = sorted({i for s in solid.surfaces for t in s.triangles for i in t})
    membership = [
        {
            i: face.distance_to(tuple(points[i - 1] * 1000)) <= tolerance_mm
            for i in boundary_nodes
        }
        for face in faces
    ]
    mapping = {index: [] for index in range(len(faces))}
    for surface in solid.surfaces:
        ids = {i for triangle in surface.triangles for i in triangle}
        eligible = [
            index
            for index, inside in enumerate(membership)
            if all(inside[i] for i in ids)
        ]
        if len(eligible) > 1:
            raise ValueError("Boundary face mapping is ambiguous.")
        if eligible:
            mapping[eligible[0]].append(surface.id)
    groups = {s.id: s for s in solid.surfaces}
    evidence = []
    for index, selected in mapping.items():
        if not selected:
            raise ValueError(
                "Captured boundary face has no complete matching mesh surface."
            )
        area = 0.0
        moment = np.zeros(3)
        for surface_id in selected:
            for triangle in groups[surface_id].triangles:
                p = points[np.array(triangle) - 1]
                weights, _ = surface_integrals(tuple(map(tuple, p)))
                area += sum(weights)
                moment += np.array(weights) @ p
        reference = captures[index]
        if abs(area / reference["area_m2"] - 1) > 1e-6:
            raise ValueError(
                "Mapped mesh area does not cover the captured face within 1e-6 relative error."
            )
        if (
            np.linalg.norm(moment / area - np.array(reference["centre_m"]))
            > 1e-6 * length
        ):
            raise ValueError("Mapped mesh centroid disagrees with the captured face.")
        evidence.append(
            {
                "source_reference": reference["source_reference"],
                "mesh_surfaces": selected,
                "mesh_area_m2": area,
                "captured_area_m2": reference["area_m2"],
                "point_distance_budget_mm": tolerance_mm,
            }
        )
    return mapping, evidence


def binding_for(solid, mapping, request):
    conditions = []
    for condition in request["conditions"]:
        references = condition["faces"]
        if not references or any(
            type(i) is not int or i not in mapping for i in references
        ):
            raise ValueError("A condition must reference captured faces.")
        surfaces = tuple(s for i in references for s in mapping[i])
        conditions.append(
            MeshCondition(
                condition["id"],
                condition["name"],
                condition["kind"],
                surfaces,
                tuple(condition.get("axes", [])),
                tuple(condition.get("values", [])),
            )
        )
    return MeshBinding(solid, tuple(conditions))


def solve(root, gmsh, ccx):
    raw = (root / "request.json").read_bytes()
    request = read_json(root / "request.json", 100_000)
    if request["format"] != "vinkulum-freecad-static-1":
        raise ValueError("Unsupported static capture.")
    threads = request["threads"]
    if type(threads) is not int or not 1 <= threads <= 64:
        raise ValueError("Choose between 1 and 64 engine threads.")
    configure_occt_threads(threads)
    body = import_body(root, request)
    faces = import_faces(root, request, float(np.linalg.norm(body.dimensions)))
    run_mesh(
        HxtMeshRequest(body, request["mesh_size_mm"], threads=threads),
        root / "mesh",
        executable=gmsh,
    )
    solid, _ = load_mesh(root / "mesh")
    mapping, evidence = match_faces(solid, faces, request["faces"])
    binding = binding_for(solid, mapping, request)
    study = binding.study(request["young_pa"], request["poisson"])
    run_static(study, root / "static", executable=ccx, threads=threads)
    admitted, result, _ = load_static_result(root / "static")
    if (root / "request.json").read_bytes() != raw:
        raise ValueError("Captured request changed during calculation.")
    preview = {
        "format": "vinkulum-freecad-static-preview-1",
        "run_id": request["run_id"],
        "request_sha256": hashlib.sha256(raw).hexdigest(),
        "nodes_m": admitted.solver_study.nodes,
        "elements": admitted.elements,
        "displacements_m": result["displacements"],
        "face_mapping": evidence,
        "strain_energy_J": result["strain_energy_J"],
        "reactions_N": result["reactions"],
        "gmsh_version": solid.gmsh_version,
        "gmsh_occt": solid.occ_version,
        "calculix_version": result["engine_version"],
        "execution": result.get("execution"),
    }
    write_json(root / "preview.json", preview)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--ccx", required=True)
    args = parser.parse_args()
    solve(args.directory, args.gmsh, args.ccx)
