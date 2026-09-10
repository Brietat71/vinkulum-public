"""Editable native static reference, ready for the existing task panel.

SPDX-License-Identifier: Apache-2.0
"""

import FreeCAD as App
import FreeCADGui as Gui
from PySide6 import QtWidgets

from . import static_analysis as model


def create_tension_reference():
    """Create the affine-tension reference without starting any external process."""
    document = App.newDocument("TensionReference")
    document.Label = "Tension reference · 2 MPa"
    document.UndoMode = 1
    try:
        source = document.addObject("Part::Box", "Bar")
        source.Label = "Bar · 120 × 20 × 15 mm"
        source.Length, source.Width, source.Height = 120, 20, 15
        document.recompute()
        analysis = model.create(source)
        analysis.Poisson = 0.0  # The reference has an exact affine displacement field.
        for target, kind, pressure in (
            (App.Vector(0, 10, 7.5), "Fixed", 0),
            (App.Vector(120, 10, 7.5), "Pressure", -2e6),
        ):
            face = min(
                range(len(source.Shape.Faces)),
                key=lambda i: (source.Shape.Faces[i].CenterOfMass - target).Length,
            )
            if (source.Shape.Faces[face].CenterOfMass - target).Length > 1e-7:
                raise ValueError("The reference end face could not be identified.")
            model.add_boundary(analysis, [f"Face{face + 1}"], kind, pressure)
        return analysis
    except Exception:
        App.closeDocument(document.Name)
        raise


def open_tension_reference():
    try:
        if Gui.Control.activeDialog():
            raise ValueError("Finish or close the current FreeCAD task first.")
        analysis = create_tension_reference()
        view = Gui.activeDocument().activeView()
        animated = view.isAnimationEnabled()
        view.setAnimationEnabled(False)
        try:
            view.viewAxonometric()
            view.fitAll(1.15)
        finally:
            view.setAnimationEnabled(animated)
        Gui.Selection.clearSelection()
        Gui.Selection.addSelection(analysis)
        return analysis
    except Exception as error:  # noqa: BLE001 - native command error boundary
        QtWidgets.QMessageBox.information(Gui.getMainWindow(), "Vinkulum", str(error))
        return None
