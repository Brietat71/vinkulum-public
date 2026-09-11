"""Numerical CAD coverage of mesher entities, separate from mesh area error.

SPDX-License-Identifier: Apache-2.0
"""

import hashlib

import build123d as bd
import numpy as np
from OCP.BRepAlgoAPI import BRepAlgoAPI_Common, BRepAlgoAPI_Cut
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.gp import gp_Trsf


def coverage(solid, directory, surface_ids, reference):
    if not np.isfinite(reference.area) or reference.area <= 0 or not reference.is_valid:
        raise ValueError("Invalid captured reference face.")
    surface_ids = tuple(surface_ids)
    if len(surface_ids) != len(set(surface_ids)):
        raise ValueError("Repeated CAD witness surface identifier.")
    remaining = reference.wrapped
    rotation = np.array(solid.request.body.orientation).reshape(3, 3)
    transform = np.eye(4)
    transform[:3, :3] = rotation
    transform[:3, 3] = np.array(solid.request.body.position) * 1000
    rigid = gp_Trsf()
    rigid.SetValues(*(float(v) for v in transform[:3, :].flat))
    hashes = dict(solid.cad_faces)
    area = 0.0
    common_area = 0.0
    moment = np.zeros(3)
    for identifier in surface_ids:
        path = directory / f"cad-face-{identifier}.brep"
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != hashes[identifier]:
            raise ValueError("CAD face witness changed before coverage validation.")
        imported = bd.import_brep(path)
        if len(imported.faces()) != 1 or imported.solids():
            raise ValueError("A CAD witness must contain exactly one face.")
        if path.read_bytes() != raw:
            raise ValueError("CAD face witness changed during import.")
        face = bd.Compound.cast(
            BRepBuilderAPI_Transform(imported.faces()[0].wrapped, rigid, True).Shape()
        )
        if not face.is_valid or not np.isfinite(face.area) or face.area <= 0:
            raise ValueError("Invalid CAD witness face.")
        intersection = BRepAlgoAPI_Common(reference.wrapped, face.wrapped)
        intersection.Build()
        if not intersection.IsDone():
            raise ValueError("CAD coverage operation failed.")
        common = bd.Compound.cast(intersection.Shape())
        covered = common.area
        if not np.isfinite(covered) or abs(covered / face.area - 1) > 1e-6:
            raise ValueError("Mesher CAD face is not covered by the captured face.")
        difference = BRepAlgoAPI_Cut(remaining, face.wrapped)
        difference.Build()
        if not difference.IsDone():
            raise ValueError("CAD uncovered-area operation failed.")
        remaining = difference.Shape()
        common_area += covered
        area += face.area
        moment += face.area * np.array(tuple(face.center(bd.CenterOf.MASS)))
    if area <= 0 or abs(common_area / reference.area - 1) > 1e-6:
        raise ValueError("Mesher CAD faces do not cover the complete captured face.")
    uncovered = bd.Compound.cast(remaining).area
    if not np.isfinite(uncovered) or uncovered / reference.area > 1e-6:
        raise ValueError("Mesher CAD faces leave part of the captured face uncovered.")
    return area * 1e-6, moment / area * 1e-3
