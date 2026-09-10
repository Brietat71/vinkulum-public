"""Reproduce the real CAD-face study workflow with Gmsh/OCCT 8 and CalculiX.

Run from an installed Studio environment under Xvfb on Linux. Every execution
uses a new output directory and preserves raw mesh/calculation files on failure.
The loaded plate is clamped at its minimum-X end and loaded by 100 kPa inward
pressure on its top face. This is a workflow example, not an analytic reference.
"""

import argparse
import json
import time
from pathlib import Path
from unittest.mock import patch

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QDialogButtonBox, QFileDialog
from vinkulum_studio import __version__
from vinkulum_studio.calculix import load_static_result, load_study
from vinkulum_studio.document import History, load_project
from vinkulum_studio.editor import EditorWindow
from vinkulum_studio.mesh_binding import mesh_measurements
from vinkulum_studio.mesh_window import ConditionDialog


def main(output, gmsh, ccx):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    app = QApplication.instance() or QApplication([])
    project = load_project(
        Path(__file__).resolve().parents[1]
        / "examples/studio/platine-percee.vinkulum.json"
    )
    editor = EditorWindow()
    editor.history = History(project)
    editor._saved = project
    editor._refresh(fit=True)
    editor.select_object(project.bodies[0].id)
    editor.show()
    editor.commands["cad_study"].trigger()
    window = editor._mesh_window
    if window is None:
        raise RuntimeError(editor.status.text())

    def wait(predicate, diagnostic, seconds=180):
        deadline = time.monotonic() + seconds
        while not predicate() and time.monotonic() < deadline:
            QTest.qWait(20)
        if not predicate():
            raise RuntimeError(diagnostic())

    window.executable.setText(gmsh)
    window.output.setText(str(output / "meshes"))
    window.size.setText("40")
    QTest.mouseClick(window.generate, Qt.MouseButton.LeftButton)
    wait(lambda: not window.busy, window.status.text)
    solid = window.solid
    if solid is None:
        raise RuntimeError(window.status.text())
    points = np.array(solid.mesh.nodes)

    def plane(axis, maximum):
        coordinate = points[:, axis].max() if maximum else points[:, axis].min()
        found = [
            f.id
            for f in solid.surfaces
            if np.max(abs(points[np.array(f.triangles) - 1, axis] - coordinate)) < 1e-10
        ]
        if len(found) != 1:
            raise RuntimeError("Expected exactly one planar boundary group.")
        return found[0]

    clamp, top = plane(0, False), plane(2, True)
    for kind, face, name in (
        ("support", clamp, "Clamped end"),
        ("pressure", top, "Top pressure"),
    ):
        window.set_selected_faces((face,))

        def fill(kind=kind, name=name):
            dialog = app.activeModalWidget()
            if not isinstance(dialog, ConditionDialog):
                raise TypeError("Boundary editor did not open.")
            dialog.name.setText(name)
            if kind == "pressure":
                dialog.pressure.setText("100000")
            QTest.mouseClick(
                dialog.findChild(QDialogButtonBox).button(
                    QDialogButtonBox.StandardButton.Ok
                ),
                Qt.MouseButton.LeftButton,
            )

        QTimer.singleShot(30, fill)
        QTest.mouseClick(window.add_buttons[kind], Qt.MouseButton.LeftButton)
    window.set_selected_faces(())
    camera = window.viewport.renderer.GetActiveCamera()
    camera.SetFocalPoint(0, 0, 0)
    camera.SetPosition(-1.0, -1.0, 0.8)
    window.viewport.fit_scene()
    QTest.qWait(200)
    if (
        not window.screen()
        .grabWindow(window.winId())
        .save(str(output / "mesh-workspace.png"))
    ):
        raise RuntimeError("Could not capture the mesh workspace.")
    with patch.object(
        QFileDialog,
        "getSaveFileName",
        return_value=(str(output / "plate.ccx.json"), "JSON"),
    ):
        QTest.mouseClick(window.save, Qt.MouseButton.LeftButton)
        wait(lambda: not window.busy, window.status.text)
    saved = load_study(output / "plate.ccx.json")
    QTest.mouseClick(window.open_calculix, Qt.MouseButton.LeftButton)
    wait(lambda: not window.busy, window.status.text)
    static = window._static_window
    if static is None or static.study != saved:
        raise RuntimeError("The solver workspace did not retain the saved study.")
    static.executable.setText(ccx)
    static.output.setText(str(output / "calculations"))
    static.timeout.setText("120")
    QTest.mouseClick(static.run_button, Qt.MouseButton.LeftButton)
    wait(lambda: static.controller.process is None, static.status.text)
    if static.result is None:
        raise RuntimeError(static.status.text())
    study, report, calculation = static.result
    if load_static_result(calculation)[1] != report or study != saved:
        raise RuntimeError("The captured calculation did not reopen identically.")
    static._auto_scale()
    static.viewport.setFocus()
    QTest.qWait(200)
    if (
        not static.screen()
        .grabWindow(static.winId())
        .save(str(output / "static-workspace.png"))
    ):
        raise RuntimeError("Could not capture the calculated displacement field.")
    data = {
        "studio_version": __version__,
        "mesh": str(window.mesh_root.relative_to(output)),
        "calculation": str(calculation.relative_to(output)),
        "nodes": len(solid.mesh.nodes),
        "elements": len(solid.mesh.elements),
        "face_groups": len(solid.surfaces),
        "clamp_face": clamp,
        "pressure_face": top,
        "energy_J": report["strain_energy_J"],
        "input_transport": report["input_transport"],
        **mesh_measurements(solid),
        "qualification": "Workflow, local geometry and balance checks; no mesh-convergence claim for this plate.",
    }
    (output / "recipe.json").write_text(json.dumps(data, indent=2) + "\n")
    print(json.dumps(data, indent=2), flush=True)
    if not editor.close():
        raise RuntimeError("The completed workspace did not close cleanly.")
    QTest.qWait(50)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--gmsh", required=True)
    parser.add_argument("--ccx", default="ccx")
    args = parser.parse_args()
    main(args.output, args.gmsh, args.ccx)
