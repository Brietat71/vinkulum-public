"""Create and edit a dimensioned bracket through real Qt dialogs and the OCCT worker.

Use an installed Studio CAD wheel under a desktop or Xvfb. Output includes
schema-4 projects, unmodified window captures and a machine-readable report.
"""

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QDialog
from vinkulum_studio import __version__
from vinkulum_studio.cad_dialog import CadDialog
from vinkulum_studio.cad_history_dialog import CadHistoryDialog
from vinkulum_studio.document import (
    History,
    Project,
    joint_at,
    load_project,
    new_id,
    save_project,
)
from vinkulum_studio.editor import EditorWindow
from vinkulum_studio.pinocchio_controller import PinocchioController
from vinkulum_studio.sketch import solve
from vinkulum_studio.sketch_dialog import SketchDialog


def require(value, message):
    if not value:
        raise RuntimeError(message)


def wait(predicate, diagnostic, seconds=30):
    deadline = time.monotonic() + seconds
    while not predicate() and time.monotonic() < deadline:
        QTest.qWait(20)
    require(predicate(), diagnostic())


def main(output, pinocchio_python=None):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    app = QApplication.instance() or QApplication([])
    editor = EditorWindow()
    editor.history = History(
        Project(new_id(), name="Dimensioned bracket", gravity=(0.0, 0.0, -9.81))
    )
    editor._saved = editor.project
    editor._refresh(fit=True)
    editor.show()
    editor.activateWindow()
    QTest.qWait(50)

    def capture(window, name):
        QTest.qWait(150)
        require(
            window.screen().grabWindow(window.winId()).save(str(output / name)),
            "Window capture failed.",
        )

    def modal(action, callback):
        errors = []
        guard = QTimer(editor)
        guard.setSingleShot(True)

        def reject():
            active = app.activeModalWidget()
            if isinstance(active, QDialog):
                active.reject()

        def expired():
            errors.append("The modal interaction timed out.")
            reject()

        def fill():
            try:
                callback(app.activeModalWidget())
            except (
                RuntimeError,
                ValueError,
                TypeError,
                OSError,
                AssertionError,
                AttributeError,
                StopIteration,
            ) as error:
                errors.append(f"{type(error).__name__}: {error}")
                reject()

        guard.timeout.connect(expired)
        guard.start(40_000)
        QTimer.singleShot(40, fill)
        action()
        guard.stop()
        require(not errors, "; ".join(errors))

    def inspect_sketch(dialog):
        require(isinstance(dialog, SketchDialog), "The sketch editor did not open.")
        width = next(c for c in dialog.profile.constraints if c.name == "Width")
        dialog.select((("constraint", width.id),))
        capture(dialog, "sketch-80.png")
        QTest.mouseClick(dialog.accept_button, Qt.MouseButton.LeftButton)

    def create(dialog):
        require(isinstance(dialog, CadDialog), "The CAD dialog did not open.")
        require(
            dialog.operation.currentData() == "extrude_sketch",
            "The sketch command selected the wrong operation.",
        )
        dialog.name.setText("Dimensioned bracket")
        dialog.sketch_height.setText("10")
        modal(
            lambda: QTest.mouseClick(dialog.profile_button, Qt.MouseButton.LeftButton),
            inspect_sketch,
        )
        QTest.mouseClick(dialog.apply_button, Qt.MouseButton.LeftButton)
        wait(lambda: dialog.process is None, dialog.status.text)
        require(dialog.result_data is not None, dialog.status.text())

    modal(
        lambda: QTest.keyClick(
            editor,
            Qt.Key.Key_G,
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier,
        ),
        create,
    )
    require(
        len(editor.project.bodies) == 1,
        "The CAD solid was not applied to the document.",
    )
    initial = editor.project.bodies[0]
    require(
        abs(initial.mass - 0.2184) < 1e-12,
        "Initial bracket mass differs from the two-cuboid reference.",
    )
    save_project(output / "bracket-80.vinkulum.json", editor.project)
    initial_feature = initial.cad.recipe.features[0]
    point_ids = tuple(p.id for p in initial_feature.profile.points)
    origin = np.array(initial.position) + initial.cad.recipe.origin_in_body_m
    history_length = len(editor.history.past)

    def edit_sketch(dialog):
        require(
            isinstance(dialog, SketchDialog), "The upstream sketch editor did not open."
        )
        width = next(c for c in dialog.profile.constraints if c.name == "Width")
        dialog.select((("constraint", width.id),))
        QTest.qWait(20)
        field = dialog.value_fields[0]
        field.setFocus()
        field.selectAll()
        QTest.keyClicks(field, "100")
        QTest.mouseClick(dialog.update_button, Qt.MouseButton.LeftButton)
        require(
            solve(dialog.profile).degrees_of_freedom == 0,
            "The edited profile lost its constraints.",
        )
        dialog.canvas.setFocus()
        QTest.keyClick(dialog.canvas, Qt.Key.Key_F)
        capture(dialog, "sketch-100.png")
        QTest.mouseClick(dialog.accept_button, Qt.MouseButton.LeftButton)

    def regenerate(dialog):
        require(
            isinstance(dialog, CadHistoryDialog), "The CAD feature editor did not open."
        )
        dialog.tree.setCurrentItem(dialog.items[initial_feature.id])
        modal(
            lambda: QTest.mouseClick(dialog.profile_button, Qt.MouseButton.LeftButton),
            edit_sketch,
        )
        QTest.mouseClick(dialog.preview_button, Qt.MouseButton.LeftButton)
        wait(lambda: dialog.controller.process is None, dialog.status.text)
        require(
            dialog.preview_body is not None and dialog.apply_button.isEnabled(),
            dialog.status.text(),
        )
        require(
            editor.project.bodies[0] == initial,
            "Preview changed the document before Apply.",
        )
        dialog.viewport.camera("iso")
        dialog.viewport.fit_scene()
        capture(dialog, "extruded-bracket.png")
        QTest.mouseClick(dialog.apply_button, Qt.MouseButton.LeftButton)

    modal(editor.commands["cad_history"].trigger, regenerate)
    updated = editor.project.bodies[0]
    require(
        abs(updated.mass - 0.2496) < 1e-12,
        "Updated bracket mass differs from the two-cuboid reference.",
    )
    require(
        updated.id == initial.id
        and updated.cad.recipe.features[0].id == initial_feature.id,
        "CAD identities changed during regeneration.",
    )
    require(
        tuple(p.id for p in updated.cad.recipe.features[0].profile.points) == point_ids,
        "Sketch point identities changed during the dimension edit.",
    )
    np.testing.assert_allclose(
        np.array(updated.position) + updated.cad.recipe.origin_in_body_m,
        origin,
        rtol=0,
        atol=1e-13,
    )
    require(
        len(editor.history.past) == history_length + 1,
        "Apply must create one document history entry.",
    )
    editor.activateWindow()
    editor.viewport.setFocus()
    QTest.qWait(20)
    QTest.keyClick(editor, Qt.Key.Key_Z, Qt.KeyboardModifier.ControlModifier)
    require(
        editor.project.bodies[0] == initial, "Undo did not restore the original solid."
    )
    QTest.keyClick(
        editor,
        Qt.Key.Key_Z,
        Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier,
    )
    require(
        editor.project.bodies[0] == updated, "Redo did not restore the edited solid."
    )
    save_project(output / "bracket-100.vinkulum.json", editor.project)
    require(
        load_project(output / "bracket-100.vinkulum.json") == editor.project,
        "The schema-4 document did not reopen identically.",
    )
    report = {
        "studio_version": __version__,
        "occt_version": updated.cad.occt_version,
        "build123d_version": updated.cad.build123d_version,
        "source_width_mm": 80,
        "edited_width_mm": 100,
        "extrusion_height_mm": 10,
        "mass_before_kg": initial.mass,
        "mass_after_kg": updated.mass,
        "volume_before_m3": initial.cad.volume_m3,
        "volume_after_m3": updated.cad.volume_m3,
        "point_identities_preserved": True,
        "feature_identity_preserved": True,
        "design_origin_preserved": True,
        "preview_does_not_apply": True,
        "one_document_undo_transaction": True,
        "schema_4_reopened": True,
        "degrees_of_freedom": solve(
            updated.cad.recipe.features[0].profile
        ).degrees_of_freedom,
        "scope": "Closed line profiles with affine X/Y dimensions; not general curved sketches or persistent face naming.",
    }
    if pinocchio_python is not None:
        project = editor.project.replace_object(
            joint_at(editor.project, "pivot", None, updated.id, axis=(0, 1, 0))
        )
        controller = PinocchioController()
        failures = []
        controller.problem.connect(failures.append)
        angle = 0.2
        try:
            controller.start(
                project,
                {"q": [angle]},
                output / "pinocchio",
                interpreter=pinocchio_python,
            )
            wait(lambda: controller.process is None, lambda: "; ".join(failures))
            require(not failures and controller.last_result is not None, str(failures))
            result = controller.last_result
            require(result.project == project, "Pinocchio changed the CAD capture.")
            # Independent two-cuboid integrals about a Y pivot at the design origin.
            pieces = (
                (0.10, 0.02, 0.01, (0.05, 0.01, 0.005)),
                (0.03, 0.04, 0.01, (0.015, 0.04, 0.005)),
            )
            masses = [7800 * x * y * z for x, y, z, _ in pieces]
            centre = sum(m * np.array(p) for m, (*_, p) in zip(masses, pieces)) / sum(
                masses
            )
            expected_mass_matrix = sum(
                m * ((x * x + z * z) / 12 + p[0] ** 2 + p[2] ** 2)
                for m, (x, y, z, p) in zip(masses, pieces)
            )
            expected_bias = (
                -sum(masses)
                * 9.81
                * (centre[0] * math.cos(angle) + centre[2] * math.sin(angle))
            )
            x, y, z = centre
            expected_position = (
                x * math.cos(angle) + z * math.sin(angle),
                y,
                -x * math.sin(angle) + z * math.cos(angle),
            )
            np.testing.assert_allclose(
                result.report["mass_matrix"],
                [[expected_mass_matrix]],
                rtol=1e-11,
                atol=1e-15,
            )
            np.testing.assert_allclose(
                result.report["intrinsic_bias"], [expected_bias], rtol=1e-11, atol=1e-14
            )
            np.testing.assert_allclose(
                result.report["bodies"][0]["position_m"],
                expected_position,
                rtol=0,
                atol=1e-13,
            )
            report["pinocchio"] = {
                "engine_version": result.report["engine_version"],
                "schema_4_capture_preserved": True,
                "mass_matrix_reference_kg_m2": expected_mass_matrix,
                "gravity_effort_reference_N_m": expected_bias,
                "com_position_reference_m": list(expected_position),
                "reference": "Independent two-cuboid integrals, fixed-base Y revolute joint at design origin, q=0.2 rad.",
            }
        finally:
            controller.shutdown()
    (output / "recipe.json").write_text(json.dumps(report, indent=2) + "\n")
    editor._saved = editor.project
    require(editor.close(), "The completed application did not close cleanly.")
    QTest.qWait(30)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--pinocchio-python",
        type=Path,
        help="Optional separate Pinocchio worker with the current Studio wheel",
    )
    args = parser.parse_args()
    main(args.output, args.pinocchio_python)
