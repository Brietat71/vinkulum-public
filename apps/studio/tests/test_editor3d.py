import os
import tempfile
import time
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1",
    "Editor integration requires a graphics-enabled test run",
)
class EditorRecipe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def wait(self, predicate, timeout=10):
        deadline = time.monotonic() + timeout
        while not predicate() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(predicate(), "Missing event before timeout")

    def dialog(self, invoke, fill):
        failures = []

        def edit():
            dialog = self.app.activeModalWidget()
            try:
                self.assertIsNotNone(dialog)
                fill(dialog)
                dialog._accept()
                self.assertEqual(dialog.result(), 1)
            except Exception as exc:
                failures.append(exc)
                if dialog:
                    dialog.reject()

        QTimer.singleShot(0, edit)
        invoke()
        if failures:
            raise failures[0]

    def type_field(self, field, text):
        from vinkulum_studio.controls import VectorField

        if isinstance(field, VectorField):
            for component, value in zip(field.components, text.split(","), strict=True):
                self.type_field(component, value.strip())
            return
        field.setFocus()
        field.selectAll()
        QTest.keyClicks(field, text)

    def test_author_save_run_edit_cancel_and_captured_scene(self):
        from vinkulum_studio.editor import EditorWindow

        window = EditorWindow()
        window._discard_allowed = lambda: True
        try:
            window.show()
            QTest.qWait(50)
            window.new_project()
            window.add_body("box")
            self.assertEqual(len(window.project.bodies), 1)
            field = window.fields["mass"]
            field.setFocus()
            field.selectAll()
            QTest.keyClicks(field, "2")
            self.assertTrue(window.dirty_fields)
            self.assertTrue(window.apply_properties())
            body = window.project.bodies[0]
            self.assertEqual(body.mass, 2.0)
            window.duration.setText("0.2")
            window.run()
            self.assertIsNotNone(window.controller.process)
            self.assertFalse(window.run_button.isEnabled())
            window._commit(window.project.replace_object(replace(body, mass=3.0)))
            self.wait(lambda: window.controller.process is None)
            self.assertIsNotNone(window.result, window.status.text())
            self.assertEqual(window.result.project.bodies[0].mass, 2.0)
            self.assertEqual(window.project.bodies[0].mass, 3.0)
            self.assertEqual(window.display_project.bodies[0].mass, 2.0)
            self.assertFalse(window.properties.isEnabled())
            self.assertEqual(window.result.project.duration, 0.2)
            self.assertIn("Previous", window.result_label.text())
            window.toggle_play()
            QTest.qWait(80)
            self.assertGreater(window.slider.value(), 0)
            window.pause()
            window.mode.setCurrentIndex(0)
            self.assertEqual(window.display_project.bodies[0].mass, 3.0)
            self.assertTrue(window.properties.isEnabled())
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "mechanism.json"
                self.assertTrue(window.save(path))
                before = window.project
                window.add_body("sphere")
                window.load(path)
                self.assertEqual(window.project, before)
                window.result.export_csv(Path(directory) / "samples.csv")
            previous = window.result
            window.run()
            window.stop()
            self.wait(lambda: window.controller.process is None)
            self.assertIs(window.result, previous)
            self.assertIn("previous", window.status.text())
            screenshot = os.environ.get("VINKULUM_EDITOR_SCREENSHOT")
            if screenshot:
                window.mode.setCurrentIndex(1)
                QTest.qWait(50)
                window.viewport.screenshot(
                    str(Path(screenshot).with_suffix(".viewport.png"))
                )
                self.assertTrue(
                    self.app.primaryScreen().grabWindow(window.winId()).save(screenshot)
                )
        finally:
            window.close()

    def test_double_pendulum_authored_through_properties_and_joint_dialogs(self):
        from vinkulum_studio.dialogs import numeric_text
        from vinkulum_studio.editor import EditorWindow, matrix_euler
        from vinkulum_studio.examples3d import double_pendulum

        reference = double_pendulum()
        window = EditorWindow()
        window._discard_allowed = lambda: True
        try:
            window.show()
            window.new_project()
            ids = []
            for body in reference.bodies:
                window.add_body("box")
                ids.append(window.selection)
                self.type_field(
                    window.fields["dimensions"], numeric_text(body.dimensions)
                )
                self.type_field(window.fields["position"], numeric_text(body.position))
                self.type_field(
                    window.fields["orientation"],
                    numeric_text(matrix_euler(body.orientation)),
                )
                self.assertTrue(window.apply_properties())
            for a, b, point in (
                (None, ids[0], (0.0, 0.0, 0.0)),
                (ids[0], ids[1], tuple(v * 2 for v in reference.bodies[0].position)),
            ):

                def fill(dialog):
                    dialog.kind.setCurrentIndex(dialog.kind.findData("pivot"))
                    dialog.a.setCurrentIndex(dialog.a.findData(a))
                    dialog.b.setCurrentIndex(dialog.b.findData(b))
                    dialog.point.setText(numeric_text(point))
                    dialog.axis.setText("0,1,0")

                self.dialog(window.add_joint_dialog, fill)
            self.assertEqual(len(window.project.joints), 2)
            self.assertFalse(window.project.diagnostics())
            window.duration.setText("0.2")
            window.run()
            self.wait(lambda: window.controller.process is None)
            self.assertIsNotNone(window.result, window.status.text())
            self.assertEqual(window.result.position.shape[1], 2)
            window.mode.setCurrentIndex(0)
            window.select_object(ids[0])
            window.add_load()

            def force_law(dialog):
                dialog.kind.setCurrentIndex(dialog.kind.findData("lineaire"))
                dialog.input.setPlainText("0, 1")

            self.dialog(lambda: window.edit_load_law("force", 0, "N"), force_law)
            self.assertEqual(window.project.loads[0].force[0].kind, "lineaire")
            window.run()
            self.wait(lambda: window.controller.process is None)
            self.assertEqual(window.result.project.loads, window.project.loads)
        finally:
            window.close()

    def test_pending_fields_selection_and_drag_survive_completion(self):
        from vinkulum_studio.document import IDENTITY
        from vinkulum_studio.editor import EditorWindow

        window = EditorWindow()
        window._discard_allowed = lambda: True
        try:
            window.show()
            window.new_project()
            window.add_body("box")
            first = window.selection
            window.add_body("sphere")
            second = window.selection
            self.type_field(window.fields["mass"], "invalid")
            window.select_object(first)
            self.assertEqual(window.selection, second)
            self.assertEqual(window.fields["mass"].text(), "invalid")
            self.type_field(window.fields["mass"], "3")
            window.move_body(second, (1, 0, 0), IDENTITY)
            self.assertEqual(window.object(second).mass, 3)
            self.assertEqual(window.object(second).position, (1, 0, 0))
            window.duration.setText("0.1")
            window.run()
            self.type_field(window.fields["mass"], "4")
            self.wait(lambda: window.controller.process is None)
            self.assertIsNotNone(window.result, window.status.text())
            self.assertEqual(window.mode.currentIndex(), 0)
            self.assertEqual(window.fields["mass"].text(), "4")
            self.assertTrue(window.dirty_fields)
            self.assertEqual(window.result.project.bodies[1].mass, 3)
            window.apply_properties()
            window.run()
            window.viewport.dragging = True
            self.wait(lambda: window.controller.process is None)
            self.assertEqual(window.mode.currentIndex(), 0)
            self.assertTrue(window.viewport.dragging)
            window.viewport.dragging = False
            self.assertEqual(window.result.project.bodies[1].mass, 4)
        finally:
            window.close()

    def test_replacing_visible_result_updates_scene_before_slider_signals(self):
        from vinkulum_studio.editor import EditorWindow
        from vinkulum_studio.examples3d import double_pendulum, slider_crank
        from vinkulum_studio.mechanism import MechanicalResult, simulate_project

        window = EditorWindow()
        window._discard_allowed = lambda: True
        try:
            window.show()
            with patch("sys.excepthook") as unhandled:
                for factory in (double_pendulum, slider_crank):
                    project = replace(factory(), duration=0.1)
                    with tempfile.TemporaryDirectory() as directory:
                        metadata = simulate_project(project, "replacement", directory)
                        result = MechanicalResult.read(
                            metadata, "replacement", project, directory
                        )
                    window._commit(project, project.bodies[0].id)
                    # Keep the previous result visible while the next arrives.
                    if window.result:
                        window.mode.setCurrentIndex(1)
                    window._completed(result)
                    window.slider.setValue(window.slider.maximum())
                    QTest.qWait(10)
                unhandled.assert_not_called()
                self.assertEqual(
                    window.viewport._body_ids, {b.id for b in project.bodies}
                )
                group = window.tree.topLevelItem(0)
                self.assertEqual(
                    [group.child(i).text(0) for i in range(group.childCount())],
                    [b.name for b in project.bodies],
                )
                self.assertEqual(window.fields["name"].text(), project.bodies[0].name)
        finally:
            window.close()


if __name__ == "__main__":
    unittest.main()
