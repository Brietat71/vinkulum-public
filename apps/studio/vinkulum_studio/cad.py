"""OCCT 8 / build123d adapter. BREP is in mm, mechanics and display are SI.

Import this module in the CAD process only. The mechanical solver consumes the
captured mass tensor, never a tessellation or a mutable OCCT shape.
"""

import json
import math
import tempfile
from dataclasses import asdict
from io import BytesIO
from pathlib import Path

import build123d as bd
import numpy as np
import OCP
from OCP.BRep import BRep_Builder
from OCP.BRepTools import BRepTools
from OCP.gp import gp_Trsf
from OCP.TopoDS import TopoDS_Shape

from .cad_data import CadGeometry
from .document import IDENTITY, Body, new_id

if int(OCP.__version__.split(".")[0]) < 8:
    raise RuntimeError("OCCT 8 minimum requis ; aucun repli vers OCCT 7.")


def location(rotation=IDENTITY, position_m=(0, 0, 0)):
    matrix = np.column_stack(
        (np.array(rotation).reshape(3, 3), np.array(position_m) * 1000)
    )
    transform = gp_Trsf()
    transform.SetValues(*matrix.flat)
    return bd.Location(transform)


def read_brep(text):
    shape = TopoDS_Shape()
    BRepTools.Read_s(shape, BytesIO(text.encode("ascii")), BRep_Builder())
    if shape.IsNull():
        raise ValueError("BREP vide ou illisible.")
    return bd.Compound.cast(shape)


def body_shape(body):
    if body.cad is not None:
        return read_brep(body.cad.brep_mm)
    d = tuple(v * 1000 for v in body.dimensions)
    if body.shape == "box":
        return bd.Box(*d)
    if body.shape == "sphere":
        return bd.Sphere(d[0])
    return bd.Cylinder(*d)


def capture_shape(
    shape,
    *,
    name,
    density,
    position=(0, 0, 0),
    orientation=IDENTITY,
    identifier=None,
    operations=(),
):
    solids = shape.solids()
    if len(solids) != 1 or not shape.is_valid:
        raise ValueError(
            "Une pièce doit être un solide fermé valide et connexe. Importez les pièces séparément."
        )
    solid = solids[0]
    volume = solid.volume
    if not math.isfinite(volume) or volume <= 1e-9:
        raise ValueError("Solide vide ou de volume trop petit.")
    if not math.isfinite(density) or not 0 < density <= 1e6:
        raise ValueError("Masse volumique positive requise, au plus 10⁶ kg/m³.")
    centre = np.array(tuple(solid.center(bd.CenterOf.MASS)))
    tensor = np.array(solid.matrix_of_inertia) / volume * 1e-6
    centered = solid.moved(bd.Location(tuple(-centre)))
    buffer = BytesIO()
    if not bd.export_brep(centered, buffer):
        raise ValueError("Échec de capture BREP.")
    size = tuple(centered.bounding_box().size)
    if min(size) <= 1e-6 or max(size) > 1e6:
        raise ValueError("Dimensions CAD hors plage (1 nm à 1 km).")
    # Deflection is a display tolerance in mm, never the source of mass properties.
    vertices, triangles = centered.tessellate(max(0.02, max(size) * 0.001), 0.15)
    data = CadGeometry(
        buffer.getvalue().decode("ascii"),
        volume * 1e-9,
        tuple(float(v) for v in tensor.flat),
        tuple(tuple(float(v) * 0.001 for v in point) for point in vertices),
        tuple(tuple(int(i) for i in row) for row in triangles),
        tuple(operations),
        OCP.__version__,
        bd.__version__,
    )
    world_centre = np.array(position) + np.array(orientation).reshape(3, 3) @ (
        centre * 0.001
    )
    return Body(
        identifier or new_id(),
        name,
        "cad",
        tuple(v * 0.001 for v in size),
        volume * 1e-9 * density,
        tuple(world_centre),
        orientation,
        cad=data,
    )


def execute(request):
    operation = request["operation"]
    if operation == "self_check":
        shape = bd.Box(100, 60, 20) - bd.Cylinder(10, 40)
        expected = 100 * 60 * 20 - math.pi * 10**2 * 20
        if not shape.is_valid or not math.isclose(
            shape.volume, expected, rel_tol=1e-10
        ):
            raise ValueError("Contrôle CAD : volume de la pièce percée incorrect.")
        with tempfile.TemporaryDirectory(prefix="vinkulum-cad-check-") as directory:
            path = Path(directory) / "check.step"
            bd.export_step(shape, path, unit=bd.Unit.MM)
            restored = bd.import_step(path)
            if not restored.is_valid or not math.isclose(
                restored.volume, expected, rel_tol=1e-10
            ):
                raise ValueError("Contrôle CAD : aller-retour STEP incorrect.")
        body = capture_shape(shape, name="Contrôle CAD", density=7800)
        return {
            "cad_check": {
                "occt_version": OCP.__version__,
                "build123d_version": bd.__version__,
                "volume_m3": body.cad.volume_m3,
                "mass_kg": body.mass,
                "triangles": len(body.cad.triangles),
                "step_round_trip": True,
            }
        }
    density = float(request.get("density", 7800))
    name = request.get("name", "Pièce CAD")
    a = Body.from_dict(request["a"]) if request.get("a") is not None else None
    if operation == "export_step":
        if a is None:
            raise ValueError("Sélectionnez un corps à exporter.")
        shape = body_shape(a).moved(location(a.orientation, a.position))
        shape.label = a.name
        with tempfile.TemporaryDirectory(prefix="vinkulum-step-") as directory:
            path = Path(directory) / "part.step"
            if not bd.export_step(shape, path, unit=bd.Unit.MM):
                raise ValueError("Échec de l’export STEP.")
            return {"step": path.read_text(), "occt_version": OCP.__version__}
    position, orientation, identifier, operations = (0, 0, 0), IDENTITY, None, ()
    if operation == "import_step":
        path = Path(request["path"])
        if path.stat().st_size > 8 * 1024 * 1024:
            raise ValueError("Import STEP limité à 8 Mo pour cette version.")
        shape = bd.import_step(path)
    elif operation in (
        "box",
        "cylinder",
        "sphere",
        "extrude_rectangle",
        "extrude_circle",
    ):
        values = tuple(float(v) for v in request["dimensions_mm"])
        expected = {
            "box": 3,
            "cylinder": 2,
            "sphere": 1,
            "extrude_rectangle": 3,
            "extrude_circle": 2,
        }[operation]
        if len(values) != expected or not all(
            math.isfinite(v) and 0.001 <= v <= 1e6 for v in values
        ):
            raise ValueError("Dimensions positives de 0,001 à 10⁶ mm requises.")
        if operation == "box":
            shape = bd.Box(*values)
        elif operation == "cylinder":
            shape = bd.Cylinder(*values)
        elif operation == "sphere":
            shape = bd.Sphere(*values)
        elif operation == "extrude_rectangle":
            shape = bd.extrude(bd.Rectangle(*values[:2]), values[2])
        else:
            shape = bd.extrude(bd.Circle(values[0]), values[1])
        position = tuple(
            float(v) * 0.001 for v in request.get("position_mm", (0, 0, 0))
        )
    elif operation in ("cut", "fuse", "common", "fillet"):
        if a is None:
            raise ValueError("Sélectionnez le corps à modifier.")
        shape = body_shape(a)
        position, orientation, identifier, name = (
            a.position,
            a.orientation,
            a.id,
            a.name,
        )
        operations = a.cad.operations if a.cad else ()
        if operation == "fillet":
            radius = float(request["radius_mm"])
            if not math.isfinite(radius) or radius <= 0:
                raise ValueError("Rayon de congé strictement positif requis.")
            shape = bd.fillet(shape.edges(), radius)
        else:
            b = Body.from_dict(request["b"])
            if b.id == a.id:
                raise ValueError("Deux corps distincts sont requis.")
            R = np.array(a.orientation).reshape(3, 3)
            relative_r = R.T @ np.array(b.orientation).reshape(3, 3)
            relative_p = R.T @ (np.array(b.position) - a.position)
            tool = body_shape(b).moved(
                location(tuple(relative_r.flat), tuple(relative_p))
            )
            shape = {
                "cut": lambda: shape - tool,
                "fuse": lambda: shape + tool,
                "common": lambda: shape & tool,
            }[operation]()
    else:
        raise ValueError("Opération CAD inconnue.")
    record = {k: v for k, v in request.items() if k not in ("a", "b", "path")}
    if request.get("b"):
        record["tool_id"] = request["b"]["id"]
    if operation == "import_step":
        record["file"] = Path(request["path"]).name
    body = capture_shape(
        shape,
        name=name,
        density=density,
        position=position,
        orientation=orientation,
        identifier=identifier,
        operations=(*operations, json.dumps(record, ensure_ascii=False)),
    )
    return {"body": asdict(body)}
