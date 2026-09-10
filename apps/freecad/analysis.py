"""Persistent native FreeCAD analysis data; restoring a document runs no solver.

SPDX-License-Identifier: Apache-2.0
"""

import os
import uuid
from contextlib import contextmanager

import FreeCAD as App


@contextmanager
def transaction(document, label):
    document.openTransaction(label)
    try:
        yield
    except Exception:
        document.abortTransaction()
        raise
    else:
        document.commitTransaction()


class MotionAnalysis:
    """All persisted state is in typed FreeCAD properties, not the Python proxy."""

    def execute(self, obj):
        pass

    def dumps(self):
        return None

    def loads(self, state):
        pass


class MotionView:
    def getIcon(self):
        return ":/icons/preferences-motion.svg"

    def doubleClicked(self, view):
        from .host import open_analysis

        open_analysis(view.Object)
        return True

    def setupContextMenu(self, view, menu):
        from .host import open_analysis

        menu.addAction("Edit motion analysis…", lambda: open_analysis(view.Object))

    def dumps(self):
        return None

    def loads(self, state):
        pass


def is_analysis(obj):
    return isinstance(getattr(obj, "Proxy", None), MotionAnalysis)


def validate(obj):
    if not is_analysis(obj) or obj.SchemaVersion != 1:
        raise ValueError("This motion analysis requires a supported document schema.")
    source = obj.Source
    if (
        source is None
        or source.Document is not obj.Document
        or not hasattr(source, "Shape")
    ):
        raise ValueError(
            "Assign a solid in this document to the analysis Source property."
        )
    return source


def create(source):
    document = source.Document
    with transaction(document, "Create Vinkulum motion analysis"):
        obj = document.addObject("App::FeaturePython", "VinkulumMotion")
        obj.Label = "Motion · " + source.Label
        properties = (
            ("App::PropertyInteger", "SchemaVersion", "Identity", "Document schema", 1),
            (
                "App::PropertyString",
                "AnalysisId",
                "Identity",
                "Stable analysis UUID",
                str(uuid.uuid4()),
            ),
            (
                "App::PropertyString",
                "BodyId",
                "Identity",
                "Stable captured body UUID",
                str(uuid.uuid4()),
            ),
            ("App::PropertyLink", "Source", "Motion", "Top-level rigid solid", source),
            (
                "App::PropertyFloat",
                "Density",
                "Motion",
                "Uniform density in kg/m³",
                7800.0,
            ),
            (
                "App::PropertyVector",
                "Pivot",
                "Motion",
                "World-frame pivot in mm",
                App.Vector(),
            ),
            (
                "App::PropertyVector",
                "Axis",
                "Motion",
                "World-frame axis direction (normalised at capture)",
                App.Vector(0, 1, 0),
            ),
            (
                "App::PropertyTime",
                "Duration",
                "Calculation",
                "Duration in seconds; at most 10 s",
                2.0,
            ),
            (
                "App::PropertyTime",
                "Step",
                "Calculation",
                "Native time step in seconds; at most 2001 samples",
                0.005,
            ),
            (
                "App::PropertyInteger",
                "Threads",
                "Calculation",
                "Explicit engine thread budget, 1 to 64",
                min(2, os.cpu_count() or 1),
            ),
            (
                "App::PropertyString",
                "LastCapture",
                "Captured result",
                "External calculation folder; use Reopen last calculation",
                "",
            ),
        )
        for kind, name, group, description, value in properties:
            obj.addProperty(kind, name, group, description)
            setattr(obj, name, value)
        for name in ("SchemaVersion", "AnalysisId", "BodyId", "LastCapture"):
            obj.setEditorMode(name, 1)
        obj.Proxy = MotionAnalysis()
        obj.ViewObject.Proxy = MotionView()
    return obj


def for_selection(selected):
    if is_analysis(selected):
        validate(selected)
        return selected
    if not hasattr(selected, "Shape"):
        raise ValueError(
            "Select one solid, PartDesign Body or Vinkulum motion analysis."
        )
    matches = [
        obj
        for obj in selected.Document.Objects
        if is_analysis(obj) and obj.Source is selected
    ]
    if len(matches) > 1:
        raise ValueError("Select the desired motion analysis in the document tree.")
    obj = matches[0] if matches else create(selected)
    validate(obj)
    return obj


def inputs(obj):
    """Read full document precision, independent of the task panel's formatting."""
    validate(obj)
    return {
        "density": obj.Density,
        "body_id": obj.BodyId,
        "project_id": obj.AnalysisId,
        "duration": obj.Duration.Value,
        "step": obj.Step.Value,
        "pivot_mm": list(obj.Pivot),
        "axis_world": list(obj.Axis),
        "threads": obj.Threads,
    }


def update(obj, name, value, component=None):
    if component is not None:
        vector = list(getattr(obj, name))
        vector[component] = value
        value = App.Vector(*vector)
    with transaction(obj.Document, "Edit Vinkulum " + name):
        setattr(obj, name, value)


def remember_capture(obj, directory):
    value = str(directory)
    if obj.LastCapture != value:
        with transaction(obj.Document, "Record Vinkulum calculation"):
            obj.LastCapture = value
