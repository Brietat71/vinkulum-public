"""Experimental FreeCAD capture of explicit boundary faces for linear statics.

SPDX-License-Identifier: Apache-2.0
This module is not yet registered as an interactive extension command.
"""

import hashlib
import json
import math
import uuid
from pathlib import Path

import bridge
import FreeCAD as App


def capture(
    source,
    directory,
    conditions,
    *,
    young_pa,
    poisson,
    mesh_size_mm,
    density=7800.0,
    threads=2,
):
    shape = source.Shape
    if not shape.isValid() or len(shape.Solids) != 1 or shape.Volume <= 0:
        raise ValueError("Select one valid solid.")
    local, world = source.Placement.toMatrix(), source.getGlobalPlacement().toMatrix()
    if any(
        abs(getattr(local, f"A{i}{j}") - getattr(world, f"A{i}{j}")) > 1e-12
        for i in range(1, 5)
        for j in range(1, 5)
    ):
        raise ValueError("Select a top-level solid.")
    if not math.isfinite(density) or density <= 0:
        raise ValueError("A finite positive density is required.")
    references = sorted({ref for c in conditions for ref in c["faces"]})
    if not references or len(references) > 32:
        raise ValueError("Capture between 1 and 32 boundary faces.")
    for ref in references:
        if type(ref) is not int or not 1 <= ref <= len(shape.Faces):
            raise ValueError("Boundary references must identify existing source faces.")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    shape.exportStep(str(directory / "part.step"))
    faces = []
    for index, ref in enumerate(references):
        face = shape.Faces[ref - 1]
        path = directory / f"face-{index:03}.step"
        face.exportStep(str(path))
        faces.append(
            {
                "source_reference": f"{source.Name}.Face{ref}",
                "step_sha256": bridge.sha(path.read_bytes()),
                "area_m2": face.Area * 1e-6,
                "centre_m": [v * 1e-3 for v in face.CenterOfMass],
            }
        )
    request = {
        "format": "vinkulum-freecad-static-1",
        "run_id": str(uuid.uuid4()),
        "body_id": str(uuid.uuid4()),
        "label": source.Label,
        "source_geometry_sha256": bridge.signature(source),
        "step_sha256": bridge.sha((directory / "part.step").read_bytes()),
        "density_kg_m3": density,
        "properties_si": {
            "volume_m3": shape.Volume * 1e-9,
            "mass_kg": shape.Volume * density * 1e-9,
            "centre_m": [v * 1e-3 for v in shape.CenterOfMass],
            "inertia_kg_m2": [
                v * density * 1e-15 for v in bridge.matrix_values(shape.MatrixOfInertia)
            ],
        },
        "young_pa": young_pa,
        "poisson": poisson,
        "mesh_size_mm": mesh_size_mm,
        "threads": threads,
        "faces": faces,
        "conditions": [
            {
                **c,
                "id": str(uuid.uuid4()),
                "faces": [references.index(ref) for ref in c["faces"]],
            }
            for c in conditions
        ],
    }
    (directory / "request.json").write_bytes(bridge.request_bytes(request))
    return request


def import_result(source, directory, request):
    """Import checked displacement data into native FreeCAD FEM objects."""
    import Fem
    import ObjectsFem

    directory = Path(directory)
    if bridge.signature(source) != request["source_geometry_sha256"]:
        raise ValueError("Source geometry changed; recapture before importing results.")
    raw = (directory / "request.json").read_bytes()
    if raw != bridge.request_bytes(request):
        raise ValueError("Captured request changed.")
    path = directory / "preview.json"
    if path.stat().st_size > 4_000_000:
        raise ValueError("Static preview exceeds the input budget.")
    data = json.loads(path.read_text())
    if (
        data["format"] != "vinkulum-freecad-static-preview-1"
        or data["run_id"] != request["run_id"]
        or data["request_sha256"] != hashlib.sha256(raw).hexdigest()
    ):
        raise ValueError("Static result belongs to another capture.")
    nodes, displacements = data["nodes_m"], data["displacements_m"]
    if not 4 <= len(nodes) == len(displacements) <= 6000:
        raise ValueError("Invalid result node count.")
    if any(
        len(row) != 3 or not all(math.isfinite(v) for v in row)
        for rows in (nodes, displacements)
        for row in rows
    ):
        raise ValueError("Nonfinite result coordinates or displacements.")
    mesh = Fem.FemMesh()
    for index, point in enumerate(nodes, 1):
        mesh.addNode(*(v * 1000 for v in point), index)
    for index, element in enumerate(data["elements"], 1):
        if len(element) != 10 or any(
            type(i) is not int or not 1 <= i <= len(nodes) for i in element
        ):
            raise ValueError("Invalid quadratic tetrahedron in result.")
        mesh.addVolume(list(element), index)
    document = source.Document
    document.openTransaction("Import captured CalculiX result")
    try:
        mesh_object = ObjectsFem.makeMeshResult(document, "VinkulumStaticMesh")
        mesh_object.FemMesh = mesh
        result = ObjectsFem.makeResultMechanical(document, "VinkulumStaticResult")
        result.Mesh = mesh_object
        result.NodeNumbers = list(range(1, len(nodes) + 1))
        result.DisplacementVectors = [
            App.Vector(*(v * 1000 for v in row)) for row in displacements
        ]
        result.DisplacementLengths = [math.hypot(*row) * 1000 for row in displacements]
        result.addProperty("App::PropertyString", "CaptureDirectory", "Vinkulum")
        result.CaptureDirectory = str(directory)
        result.addProperty("App::PropertyLink", "Source", "Vinkulum")
        result.Source = source
        result.addProperty("App::PropertyString", "SourceGeometrySha256", "Vinkulum")
        result.SourceGeometrySha256 = request["source_geometry_sha256"]
        result.addProperty("App::PropertyString", "ScientificScope", "Vinkulum")
        result.ScientificScope = (
            "Linear statics; captured mesh; no general FE error bound"
        )
        document.commitTransaction()
    except Exception:
        document.abortTransaction()
        raise
    return result, data
