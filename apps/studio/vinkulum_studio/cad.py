"""OCCT 8 / build123d adapter. BREP is in mm, mechanics and display are SI.

Import this module in the CAD process only. The mechanical solver consumes the
captured mass tensor, never a tessellation or a mutable OCCT shape.
"""

import hashlib
import json
import math
import tempfile
from dataclasses import asdict, replace
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
from .cad_history import CadFeature, CadRecipe, recipe_for_body
from .document import IDENTITY, Body, new_id

if int(OCP.__version__.split(".")[0]) < 8:
    raise RuntimeError("OCCT 8 or later is required; OCCT 7 is not supported.")


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
        raise ValueError("Empty or unreadable BREP.")
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


def brep_text(shape):
    buffer = BytesIO()
    if not bd.export_brep(shape, buffer):
        raise ValueError("BREP capture failed.")
    return buffer.getvalue().decode("ascii")


def regenerate_shape(recipe):
    """Evaluate a solid graph by stable input identity in its design frame [mm]."""
    shapes = {}
    for feature in recipe.ordered():
        try:
            kind, values = feature.kind, feature.dimensions_mm
            if kind == "box":
                shape = bd.Box(*values)
            elif kind == "cylinder":
                shape = bd.Cylinder(*values)
            elif kind == "sphere":
                shape = bd.Sphere(*values)
            elif kind == "extrude_rectangle":
                shape = bd.extrude(bd.Rectangle(*values[:2]), values[2])
            elif kind == "extrude_circle":
                shape = bd.extrude(bd.Circle(values[0]), values[1])
            elif kind == "snapshot":
                shape = read_brep(feature.brep_mm)
            elif kind == "fillet":
                base = shapes[feature.inputs[0]]
                shape = bd.fillet(base.edges(), values[0])
            else:
                a, b = (shapes[key] for key in feature.inputs)
                shape = a - b if kind == "cut" else a + b if kind == "fuse" else a & b
            if not feature.inputs:
                shape = shape.moved(
                    location(
                        feature.orientation,
                        tuple(v * 0.001 for v in feature.position_mm),
                    )
                )
            if len(shape.solids()) != 1 or not shape.is_valid or shape.volume <= 1e-9:
                raise ValueError("The feature must produce one valid nonempty solid.")
            shapes[feature.id] = shape
        except Exception as error:
            raise ValueError(
                f"Feature '{feature.name}' ({feature.id}) failed: {error}"
            ) from error
    return shapes[recipe.root]


def capture_shape(
    shape,
    *,
    name,
    density,
    position=(0, 0, 0),
    orientation=IDENTITY,
    identifier=None,
    operations=(),
    recipe=None,
):
    solids = shape.solids()
    if len(solids) != 1 or not shape.is_valid:
        raise ValueError(
            "A part must be one valid, closed, connected solid. Import separate parts individually."
        )
    solid = solids[0]
    volume = solid.volume
    if not math.isfinite(volume) or volume <= 1e-9:
        raise ValueError("Empty solid or volume too small.")
    if not math.isfinite(density) or not 0 < density <= 1e6:
        raise ValueError("Density must be positive and at most 10⁶ kg/m³.")
    centre = np.array(tuple(solid.center(bd.CenterOf.MASS)))
    tensor = np.array(solid.matrix_of_inertia) / volume * 1e-6
    centered = solid.moved(bd.Location(tuple(-centre)))
    size = tuple(centered.bounding_box().size)
    if min(size) <= 1e-6 or max(size) > 1e6:
        raise ValueError("CAD dimensions out of range (1 nm to 1 km).")
    # Deflection is a display tolerance in mm, never the source of mass properties.
    vertices, triangles = centered.tessellate(max(0.02, max(size) * 0.001), 0.15)
    data = CadGeometry(
        brep_text(centered),
        volume * 1e-9,
        tuple(float(v) for v in tensor.flat),
        tuple(tuple(float(v) * 0.001 for v in point) for point in vertices),
        tuple(tuple(int(i) for i in row) for row in triangles),
        tuple(operations),
        OCP.__version__,
        bd.__version__,
        None
        if recipe is None
        else replace(recipe, origin_in_body_m=tuple(float(v) * -0.001 for v in centre)),
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
            raise ValueError("CAD check: incorrect perforated-part volume.")
        with tempfile.TemporaryDirectory(prefix="vinkulum-cad-check-") as directory:
            path = Path(directory) / "check.step"
            bd.export_step(shape, path, unit=bd.Unit.MM)
            restored = bd.import_step(path)
            if not restored.is_valid or not math.isclose(
                restored.volume, expected, rel_tol=1e-10
            ):
                raise ValueError("CAD check: STEP round-trip failed.")
        body = capture_shape(shape, name="CAD check", density=7800)
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
    name = request.get("name", "CAD part")
    a = Body.from_dict(request["a"]) if request.get("a") is not None else None
    if operation == "export_step":
        if a is None:
            raise ValueError("Select a body to export.")
        shape = body_shape(a).moved(location(a.orientation, a.position))
        shape.label = a.name
        with tempfile.TemporaryDirectory(prefix="vinkulum-step-") as directory:
            path = Path(directory) / "part.step"
            if not bd.export_step(shape, path, unit=bd.Unit.MM):
                raise ValueError("STEP export failed.")
            return {"step": path.read_text(), "occt_version": OCP.__version__}
    position, orientation, identifier, operations = (0, 0, 0), IDENTITY, None, ()
    recipe = None
    if operation == "import_step":
        path = Path(request["path"])
        if path.stat().st_size > 8 * 1024 * 1024:
            raise ValueError("STEP import is limited to 8 MB in this version.")
        shape = bd.import_step(path)
        # Imported file coordinates may be far from the origin. Capture the
        # solid in a local design frame while retaining its original world pose.
        centre = np.array(tuple(shape.center(bd.CenterOf.MASS)))
        shape = shape.moved(bd.Location(tuple(-centre)))
        position = tuple(float(v) * 0.001 for v in centre)
        feature = CadFeature(
            new_id(), Path(path).name[:128], "snapshot", brep_mm=brep_text(shape)
        )
        recipe = CadRecipe((feature,), feature.id)
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
            raise ValueError("Positive dimensions from 0.001 to 10⁶ mm are required.")
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
        feature = CadFeature(new_id(), name, operation, dimensions_mm=values)
        recipe = CadRecipe((feature,), feature.id)
    elif operation in ("cut", "fuse", "common", "fillet", "regenerate"):
        if a is None:
            raise ValueError("Select the body to modify.")
        recipe = recipe_for_body(a)
        position, orientation, identifier, name = (
            a.position,
            a.orientation,
            a.id,
            a.name,
        )
        operations = a.cad.operations if a.cad else ()
        origin = np.array(recipe.origin_in_body_m)
        R = np.array(a.orientation).reshape(3, 3)
        position = tuple(np.array(a.position) + R @ origin)
        if operation == "regenerate":
            updated = CadRecipe.from_dict(request["recipe"])
            if updated.origin_in_body_m != recipe.origin_in_body_m:
                raise ValueError(
                    "Regeneration must preserve the captured design origin."
                )
            recipe = updated
        if operation == "fillet":
            radius = float(request["radius_mm"])
            if not math.isfinite(radius) or radius <= 0:
                raise ValueError("Fillet radius must be positive.")
            feature = CadFeature(
                new_id(),
                "Fillet · all edges",
                "fillet",
                inputs=(recipe.root,),
                dimensions_mm=(radius,),
            )
            recipe = replace(
                recipe, features=(*recipe.features, feature), root=feature.id
            )
        elif operation in ("cut", "fuse", "common"):
            b = Body.from_dict(request["b"])
            if b.id == a.id:
                raise ValueError("Two distinct bodies are required.")
            relative_r = R.T @ np.array(b.orientation).reshape(3, 3)
            relative_p = R.T @ (np.array(b.position) - a.position) - origin
            tool = CadFeature(
                new_id(),
                ("Captured tool · " + b.name)[:128],
                "snapshot",
                position_mm=tuple(float(v) * 1000 for v in relative_p),
                orientation=tuple(float(v) for v in relative_r.flat),
                brep_mm=brep_text(body_shape(b)),
                source_body_id=b.id,
                source_body_sha256=hashlib.sha256(
                    json.dumps(asdict(b), sort_keys=True, allow_nan=False).encode()
                ).hexdigest(),
            )
            feature = CadFeature(
                new_id(),
                {"cut": "Subtract", "fuse": "Union", "common": "Intersection"}[
                    operation
                ],
                operation,
                inputs=(recipe.root, tool.id),
            )
            recipe = replace(
                recipe, features=(*recipe.features, tool, feature), root=feature.id
            )
        shape = regenerate_shape(recipe)
    else:
        raise ValueError("Unknown CAD operation.")
    record = {k: v for k, v in request.items() if k not in ("a", "b", "path", "recipe")}
    if operation == "regenerate":
        record["recipe_sha256"] = recipe.fingerprint()
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
        recipe=recipe,
    )
    return {"body": asdict(body)}
