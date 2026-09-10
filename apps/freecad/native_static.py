"""Read native FreeCAD FEM inputs for the experimental linear-static bridge.

SPDX-License-Identifier: Apache-2.0
Uses native FEM document objects; defines no custom input widgets or solver UI.
"""

import math
import re
from pathlib import Path

import FreeCAD as App
import bridge
import static_bridge


def read_inputs(analysis, source):
    """Return deterministic SI inputs, refusing unsupported analysis members."""
    if analysis.TypeId != "Fem::FemAnalysis" or analysis.Document != source.Document:
        raise ValueError("Select a native FEM analysis in the source document.")
    materials, conditions = [], []
    for obj in analysis.Group:
        if getattr(obj, "Suppressed", False):
            raise ValueError("Suppressed FEM members are not supported.")
        if (
            obj.TypeId == "App::MaterialObjectPython"
            and getattr(obj.Proxy, "Type", None) == "Fem::MaterialCommon"
        ):
            materials.append(obj)
            continue
        if obj.TypeId not in ("Fem::ConstraintFixed", "Fem::ConstraintPressure"):
            raise ValueError("Unsupported FEM analysis member: " + obj.Name)
        if getattr(obj, "EnableAmplitude", False):
            raise ValueError("Time-dependent FEM amplitudes are not supported.")
        faces = []
        for owner, names in obj.References:
            if owner != source or not names:
                raise ValueError(
                    "Every boundary must reference faces of the source solid."
                )
            for name in names:
                match = re.fullmatch(r"Face([1-9][0-9]*)", name)
                if match is None or int(match[1]) > len(source.Shape.Faces):
                    raise ValueError("Only existing source faces are supported.")
                faces.append(int(match[1]))
        if not faces or len(faces) != len(set(faces)):
            raise ValueError("Boundary faces must be nonempty and unique.")
        condition = {"name": obj.Name, "faces": sorted(faces)}
        if obj.TypeId == "Fem::ConstraintFixed":
            condition.update(kind="support", axes=[1, 2, 3])
        else:
            pressure = obj.Pressure.getValueAs("Pa").Value
            if not math.isfinite(pressure) or pressure <= 0:
                raise ValueError("Use positive pressure and Reversed for tension.")
            condition.update(
                kind="pressure", values=[-pressure if obj.Reversed else pressure]
            )
        conditions.append(condition)
    if len(materials) != 1:
        raise ValueError("Exactly one native solid material is required.")
    material = materials[0]
    if material.Category != "Solid" or material.References:
        raise ValueError("Use one solid material applied to the entire analysis.")
    try:
        young = (
            App.Units.Quantity(material.Material["YoungsModulus"])
            .getValueAs("Pa")
            .Value
        )
        density = (
            App.Units.Quantity(material.Material["Density"]).getValueAs("kg/m^3").Value
        )
        poisson = float(material.Material["PoissonRatio"])
    except (KeyError, ValueError, TypeError) as error:
        raise ValueError(
            "Material needs YoungsModulus, PoissonRatio and Density with valid units."
        ) from error
    if not all(math.isfinite(v) for v in (young, density, poisson)) or not (
        young > 0 and density > 0 and -1 < poisson < 0.5
    ):
        raise ValueError("Invalid linear isotropic material constants.")
    if {c["kind"] for c in conditions} != {"support", "pressure"}:
        raise ValueError("At least one fixed boundary and one pressure are required.")
    return {
        "analysis": analysis.Name,
        "source": source.Name,
        "material": material.Name,
        "young_pa": young,
        "poisson": poisson,
        "density": density,
        "conditions": sorted(conditions, key=lambda c: c["name"]),
    }


def capture(analysis, source, directory, *, mesh_size_mm, threads=2):
    inputs = read_inputs(analysis, source)
    request = static_bridge.capture(
        source,
        directory,
        inputs["conditions"],
        young_pa=inputs["young_pa"],
        poisson=inputs["poisson"],
        density=inputs["density"],
        mesh_size_mm=mesh_size_mm,
        threads=threads,
    )
    request["native_fem_inputs"] = inputs
    (Path(directory) / "request.json").write_bytes(bridge.request_bytes(request))
    return request


def import_result(analysis, source, directory, request):
    if read_inputs(analysis, source) != request.get("native_fem_inputs"):
        raise ValueError(
            "Native FEM inputs changed; recapture before importing results."
        )
    return static_bridge.import_result(source, directory, request)
