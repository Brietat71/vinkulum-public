"""Persistent native linear-static inputs and conservative result invalidation.

SPDX-License-Identifier: Apache-2.0
"""

import json
import math

from . import bridge
from .analysis import transaction


class StaticAnalysis:
    def execute(self, obj):
        refresh_result(obj)

    def dumps(self):
        return None

    def loads(self, state):
        pass


class Boundary:
    def execute(self, obj):
        pass

    def dumps(self):
        return None

    def loads(self, state):
        pass


class StaticView:
    def getIcon(self):
        return ":/icons/FEM_Analysis.svg"

    def doubleClicked(self, view):
        from .static_host import open_analysis

        open_analysis(view.Object)
        return True

    def setupContextMenu(self, view, menu):
        from .static_host import open_analysis

        menu.addAction("Edit static analysis…", lambda: open_analysis(view.Object))

    def claimChildren(self):
        return list(self.Object.Boundaries)

    def attach(self, view):
        self.Object = view.Object

    def dumps(self):
        return None

    def loads(self, state):
        pass


def is_analysis(obj):
    return isinstance(getattr(obj, "Proxy", None), StaticAnalysis)


def create(source):
    if not hasattr(source, "Shape") or len(source.Shape.Solids) != 1:
        raise ValueError("Select one solid or PartDesign Body.")
    with transaction(source.Document, "Create static analysis"):
        obj = source.Document.addObject("App::FeaturePython", "VinkulumStatic")
        obj.Label = "Static · " + source.Label
        for kind, name, value in (
            ("App::PropertyInteger", "SchemaVersion", 1),
            ("App::PropertyLink", "Source", source),
            ("App::PropertyLinkList", "Boundaries", []),
            ("App::PropertyFloat", "YoungPa", 210e9),
            ("App::PropertyFloat", "Poisson", 0.3),
            ("App::PropertyFloat", "MeshSizeMm", 10.0),
            ("App::PropertyInteger", "Threads", 2),
            ("App::PropertyString", "LastCapture", ""),
            ("App::PropertyString", "CapturedInputsSha256", ""),
            ("App::PropertyLink", "Result", None),
            ("App::PropertyString", "ResultState", "No result"),
            ("App::PropertyBool", "SourceHidden", False),
            ("App::PropertyBool", "OriginalSourceVisibility", True),
        ):
            obj.addProperty(kind, name, "Vinkulum")
            setattr(obj, name, value)
        for name in (
            "SchemaVersion",
            "LastCapture",
            "CapturedInputsSha256",
            "Result",
            "ResultState",
            "SourceHidden",
            "OriginalSourceVisibility",
        ):
            obj.setEditorMode(name, 1)
        obj.Proxy = StaticAnalysis()
        obj.ViewObject.Proxy = StaticView()
    return obj


def for_selection(selected):
    if is_analysis(selected):
        return selected
    if selected is None or not hasattr(selected, "Shape"):
        raise ValueError("Select one solid or static analysis first.")
    matches = [
        o for o in selected.Document.Objects if is_analysis(o) and o.Source is selected
    ]
    if len(matches) > 1:
        raise ValueError("Select the desired static analysis in the tree.")
    return matches[0] if matches else create(selected)


def selected_faces(obj, selection):
    refs = []
    for row in selection:
        if row.Object is not obj.Source:
            raise ValueError("Select faces belonging only to this analysis solid.")
        refs.extend(row.SubElementNames)
    refs = sorted(set(refs))
    if not refs or any(
        not ref.startswith("Face")
        or not ref[4:].isdigit()
        or not 1 <= int(ref[4:]) <= len(obj.Source.Shape.Faces)
        for ref in refs
    ):
        raise ValueError("Select one or more faces in the 3D view (Ctrl for several).")
    return refs


def add_boundary(obj, refs, kind, pressure_pa=0):
    if kind not in ("Fixed", "Pressure") or not math.isfinite(pressure_pa):
        raise ValueError("Choose a fixed support or finite pressure.")
    if not refs or any(
        r not in {f"Face{i + 1}" for i in range(len(obj.Source.Shape.Faces))}
        for r in refs
    ):
        raise ValueError("Choose existing source faces.")
    with transaction(obj.Document, "Add static boundary condition"):
        row = obj.Document.addObject("App::FeaturePython", "VinkulumBoundary")
        row.Label = kind + " · " + ", ".join(refs)
        row.addProperty("App::PropertyLinkSub", "Faces", "Boundary")
        row.Faces = (obj.Source, refs)
        row.addProperty("App::PropertyString", "Kind", "Boundary")
        row.Kind = kind
        row.addProperty("App::PropertyFloat", "PressurePa", "Boundary")
        row.PressurePa = pressure_pa
        row.addProperty("App::PropertyString", "GeometrySha256", "Boundary")
        row.GeometrySha256 = bridge.signature(obj.Source)
        row.Proxy = Boundary()
        for name in ("Faces", "Kind", "GeometrySha256"):
            row.setEditorMode(name, 1)
        obj.Boundaries = [*obj.Boundaries, row]
    refresh_result(obj)
    return row


def inputs(obj):
    if not is_analysis(obj) or obj.SchemaVersion != 1:
        raise ValueError("Unsupported static analysis document schema.")
    source = obj.Source
    if source is None or source.Document is not obj.Document:
        raise ValueError("Assign an existing solid in this document.")
    digest = bridge.signature(source)
    conditions = []
    for row in obj.Boundaries:
        linked, refs = row.Faces
        if linked is not source or row.GeometrySha256 != digest:
            raise ValueError(
                "Geometry changed: remove and reselect the boundary faces."
            )
        if not refs or any(
            not r.startswith("Face")
            or not r[4:].isdigit()
            or not 1 <= int(r[4:]) <= len(source.Shape.Faces)
            for r in refs
        ):
            raise ValueError("A boundary references a missing face.")
        condition = {"name": row.Label, "faces": [int(r[4:]) for r in refs]}
        if row.Kind == "Fixed":
            condition.update(kind="support", axes=[1, 2, 3])
        elif row.Kind == "Pressure" and math.isfinite(row.PressurePa):
            condition.update(kind="pressure", values=[row.PressurePa])
        else:
            raise ValueError("Invalid boundary condition.")
        conditions.append(condition)
    if not any(c["kind"] == "support" for c in conditions) or not any(
        c["kind"] == "pressure" for c in conditions
    ):
        raise ValueError("Add at least one fixed support and one pressure.")
    if not math.isfinite(obj.YoungPa) or obj.YoungPa <= 0 or not -1 < obj.Poisson < 0.5:
        raise ValueError("Use positive finite Young modulus and −1 < Poisson < 0.5.")
    if (
        not math.isfinite(obj.MeshSizeMm)
        or obj.MeshSizeMm <= 0
        or not 1 <= obj.Threads <= 64
    ):
        raise ValueError("Use positive finite mesh size and 1–64 engine threads.")
    return {
        "conditions": conditions,
        "young_pa": obj.YoungPa,
        "poisson": obj.Poisson,
        "mesh_size_mm": obj.MeshSizeMm,
        "threads": obj.Threads,
    }


def signature(obj):
    return bridge.sha(
        json.dumps(
            {"geometry": bridge.signature(obj.Source), "inputs": inputs(obj)},
            sort_keys=True,
            allow_nan=False,
        ).encode()
    )


def refresh_result(obj):
    if not obj.Result:
        if obj.ResultState != "No result":
            obj.ResultState = "No result"
        restore_source(obj)
        return
    try:
        current = signature(obj) == obj.CapturedInputsSha256
    except Exception:  # noqa: BLE001 - invalid or deleted native dependencies make results stale
        current = False
    value = "Current capture" if current else "Out of date — recalculate"
    if obj.ResultState != value:
        obj.ResultState = value
    if not current and obj.Result.Mesh:
        obj.Result.Mesh.Visibility = False
        restore_source(obj)


def restore_source(obj):
    if obj.SourceHidden and obj.Source:
        obj.Source.Visibility = obj.OriginalSourceVisibility
        obj.SourceHidden = False
