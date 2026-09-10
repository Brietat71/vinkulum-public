"""Desktop acceptance checks for the 2026 workspace, not just widget construction."""

import os
import tempfile
import time
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QKeySequence
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QMessageBox


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires a desktop OpenGL context"
)
class WorkspaceRecipe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.exception_patch = patch("sys.excepthook")
        self.unhandled = self.exception_patch.start()
        self.addCleanup(self.exception_patch.stop)
        self.addCleanup(self.assert_no_qt_exceptions)

    def assert_no_qt_exceptions(self):
        import gc

        gc.collect()
        QTest.qWait(10)
        self.unhandled.assert_not_called()

    def make_window(self, settings=None):
        from vinkulum_studio.editor import EditorWindow

        window = EditorWindow(settings)
        window._discard_allowed = lambda: True
        window.show()
        QTest.qWait(30)
        self.addCleanup(window.close)
        return window

    def test_palette_keyboard_search_and_disabled_commands(self):
        from vinkulum_studio.commands import CommandPalette

        window = self.make_window()
        window.new_project()
        palette = CommandPalette(window.commands.values(), window)
        palette.show()
        QTest.keyClicks(palette.query, "sphere")
        self.assertEqual(palette.results.count(), 1)
        QTest.keyClick(palette.query, Qt.Key.Key_Return)
        self.assertEqual(window.project.bodies[0].shape, "sphere")
        self.assertFalse(palette.isVisible())
        self.assertEqual(
            window.commands["save"].shortcut(),
            QKeySequence(QKeySequence.StandardKey.Save),
        )
        # Keyboard editing a number must not invoke the single-key camera shortcuts.
        window.fields["mass"].setFocus()
        window.fields["mass"].selectAll()
        QTest.keyClicks(window.fields["mass"], "137")
        self.assertEqual(window.fields["mass"].text(), "137")
        self.assertTrue(window.apply_properties())
        palette = CommandPalette(window.commands.values(), window)
        palette.show()
        palette.query.setText("stop")
        self.assertIn("unavailable", palette.results.item(0).text())
        with patch.object(window.controller, "cancel") as cancel:
            QTest.keyClick(palette.query, Qt.Key.Key_Return)
            cancel.assert_not_called()
        palette.close()

    def test_statics_shortcut_opens_window_and_commands_leave_menu_mnemonics_free(self):
        window = self.make_window()
        window.raise_()
        window.activateWindow()
        window.viewport.setFocus()
        QTest.qWait(50)
        QTest.keyClick(
            window,
            Qt.Key.Key_E,
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier,
        )
        QTest.qWait(100)
        self.assertIsNotNone(window._static_window)
        static = window._static_window
        static._discard_allowed = lambda: True
        self.assertTrue(static.isVisible())
        static.close()
        self.assertIsNone(window._static_window)
        mnemonics = {
            QKeySequence.mnemonic(action.text()).toString()
            for action in window.menuBar().actions()
        } - {""}
        for name, action in window.commands.items():
            for shortcut in action.shortcuts():
                self.assertNotIn(shortcut.toString(), mnemonics, name)

    def test_vector_display_never_rounds_untouched_model_components(self):
        window = self.make_window()
        window.new_project()
        window.add_body("box")
        exact = (1.234567890123456, -0.123456789876543, 3.141592653589793)
        body = replace(window.project.bodies[0], position=exact)
        window._commit(window.project.replace_object(body))
        vector = window.fields["position"]
        self.assertNotEqual(float(vector.components[0].text()), exact[0])
        component = vector.components[1]
        component.setFocus()
        component.selectAll()
        QTest.keyClicks(component, "-2")
        self.assertTrue(window.apply_properties())
        self.assertEqual(window.project.bodies[0].position, (exact[0], -2.0, exact[2]))
        window.undo()
        self.assertEqual(window.project.bodies[0].position, exact)
        vector = window.fields["position"]
        vector.components[0].setFocus()
        vector.components[0].selectAll()
        QTest.keyClicks(vector.components[0], "nan")
        before = window.project
        self.assertFalse(window.apply_properties())
        self.assertEqual(window.project, before)
        self.assertEqual(vector.components[0].text(), "nan")

    def test_inspector_resize_and_focus_preserve_exact_numbers(self):
        from vinkulum_studio.editor import euler_matrix

        window = self.make_window()
        window.resize(1280, 844)
        window.activateWindow()
        self.assertTrue(QTest.qWaitForWindowActive(window))
        window.new_project()
        window.add_body("box")
        exact = (-0.00134987654321, -2.01565123456e-19, 1234567.891234)
        body = replace(
            window.project.bodies[0],
            position=exact,
            mass=1.3267053771417465,
            orientation=euler_matrix((13.123456789, -21.2, 35.7)),
        )
        window._commit(window.project.replace_object(body))
        for theme in ("light", "dark"):
            window.set_theme(theme)
            for width in (280, 320, 420):
                window.resizeDocks(
                    [window.docks["inspector"]], [width], Qt.Orientation.Horizontal
                )
                window.fields["name"].setFocus()
                QTest.qWait(20)
                for key in ("position", "orientation"):
                    for field in window.fields[key].components:
                        self.assertNotEqual(field.text(), "…")
                        self.assertLessEqual(
                            field.fontMetrics().horizontalAdvance(field.text()),
                            field.available_width(),
                        )
                        field.setFocus()
                        QTest.qWait(1)
                        self.assertTrue(field.hasFocus())
                        self.assertEqual(float(field.text()), field.value())
                        window.fields["name"].setFocus()
                self.assertFalse(window.dirty_fields)
                self.assertEqual(window.fields["position"].values(), exact)
        # Applying an unrelated edit must not round the mass or reconstruct SO(3).
        window.fields["name"].selectAll()
        QTest.keyClicks(window.fields["name"], "Precision")
        self.assertTrue(window.apply_properties())
        self.assertEqual(window.project.bodies[0], replace(body, name="Precision"))
        # A new high-precision value also survives compact preview and refocus.
        field = window.fields["position"].components[1]
        field.setFocus()
        field.selectAll()
        QTest.keyClicks(field, "-2.123456789123456e-19")
        window.fields["name"].setFocus()
        self.assertEqual(field.value(), -2.123456789123456e-19)
        self.assertTrue(window.apply_properties())
        self.assertEqual(window.project.bodies[0].position[1], -2.123456789123456e-19)

    def test_english_labels_preserve_saved_joint_and_law_types(self):
        from vinkulum_studio.dialogs import JointDialog, LawDialog
        from vinkulum_studio.document import Law
        from vinkulum_studio.labels import JOINT_LABELS, LAW_LABELS

        window = self.make_window()
        window.new_project()
        window.add_body("box")
        for value, label in JOINT_LABELS.items():
            dialog = JointDialog(window.project, window.selection, window)
            dialog.kind.setCurrentIndex(dialog.kind.findData(value))
            self.assertEqual(dialog.kind.currentText(), label)
            dialog._accept()
            self.assertEqual(dialog.joint.kind, value)
            self.assertEqual(dialog.joint.name, label)
            dialog.deleteLater()
        for law in (Law(), Law("lineaire", (0, 1)), Law("table", (0, 0, 1, 2))):
            dialog = LawDialog(law, parent=window)
            self.assertEqual(dialog.kind.currentText(), LAW_LABELS[law.kind])
            self.assertEqual(dialog.read(), law)
            dialog.deleteLater()

    def test_search_visibility_isolation_and_transform_are_view_state(self):
        window = self.make_window()
        window.new_project()
        window.add_body("box")
        box = window.selection
        window.add_body("sphere")
        sphere = window.selection
        original = window.project
        window.tree_search.setText("sphere")
        group = window.tree.topLevelItem(0)
        self.assertTrue(group.child(0).isHidden())
        self.assertFalse(group.child(1).isHidden())
        window.isolate_selection()
        self.assertFalse(window.viewport._actors[box][0].GetVisibility())
        self.assertTrue(window.viewport._actors[sphere][0].GetVisibility())
        window.toggle_visibility()
        self.assertFalse(window.viewport._actors[sphere][0].GetVisibility())
        window.show_all_objects()
        self.assertTrue(
            all(
                actor.GetVisibility()
                for actors in window.viewport._actors.values()
                for actor in actors
            )
        )
        window.transform_mode.setCurrentIndex(1)
        self.assertTrue(window.viewport.box.GetTranslationEnabled())
        self.assertFalse(window.viewport.box.GetRotationEnabled())
        window.transform_mode.setCurrentIndex(2)
        self.assertFalse(window.viewport.box.GetTranslationEnabled())
        self.assertTrue(window.viewport.box.GetRotationEnabled())
        self.assertEqual(window.project, original)

    def test_native_time_sync_and_readonly_result_after_pending_edit(self):
        window = self.make_window()
        window.duration.setText("0.1")
        window.run()
        deadline = time.monotonic() + 15
        while window.controller.process is not None and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertIsNotNone(window.result, window.status.text())
        self.assertEqual(window.workspace_tabs.currentIndex(), 2)
        times = window.result.time
        index = 7
        window.time_input.setValue(float(times[index]) + 0.0001)
        self.assertEqual(window.slider.value(), index)
        self.assertEqual(window.time_input.value(), float(times[index]))
        self.assertEqual(window.sample_table.currentIndex().row(), index)
        self.assertEqual(window.curve.index, index)
        window.curve.setFocus()
        QTest.keyClick(window.curve, Qt.Key.Key_Right)
        self.assertEqual(window.slider.value(), index + 1)
        self.assertEqual(window.sample_model.rowCount(), len(times))
        self.assertIn(
            "[m]", window.sample_model.headerData(2, Qt.Orientation.Horizontal)
        )
        window.workspace_tabs.setCurrentIndex(0)
        window.select_object(window.project.bodies[0].id)
        window.fields["mass"].setFocus()
        window.fields["mass"].selectAll()
        QTest.keyClicks(window.fields["mass"], "invalid")
        window.workspace_tabs.setCurrentIndex(2)
        self.assertEqual(window.mode.currentIndex(), 0)
        self.assertEqual(window.fields["mass"].text(), "invalid")
        self.assertFalse(window.timer.isActive())
        window.toggle_play()
        self.assertFalse(window.timer.isActive())

    def test_comparison_keeps_run_identities_samples_and_model_differences(self):
        window = self.make_window()
        window.duration.setText("0.1")

        def calculate():
            window.run()
            deadline = time.monotonic() + 15
            while window.controller.process is not None and time.monotonic() < deadline:
                QTest.qWait(10)
            self.assertIsNotNone(window.result, window.status.text())
            return window.result

        first = calculate()
        window.mode.setCurrentIndex(0)
        body = window.project.bodies[0]
        window._commit(
            window.project.replace_object(
                replace(body, mass=body.mass * 2, name="Bras renommé")
            )
        )
        window.step.setText("0.0025")
        second = calculate()
        self.assertNotEqual(first.run_id, second.run_id)
        self.assertEqual(window.run_combo.count(), 2)
        window.compare_combo.setCurrentIndex(
            window.compare_combo.findData(first.run_id)
        )
        self.assertIsNotNone(window.curve.reference)
        self.assertIs(window.curve.reference[0], first.time)
        self.assertIs(window.curve.times, second.time)
        self.assertIn("mass", window.comparison_label.text())
        self.assertIn("time step", window.comparison_label.text())
        self.assertIn(first.run_id[:8], window.comparison_label.text())
        self.assertIn(second.run_id[:8], window.comparison_label.text())
        window.run_combo.setCurrentIndex(window.run_combo.findData(first.run_id))
        self.assertIs(window.result, first)
        self.assertEqual(window.display_project.bodies[0].mass, body.mass)
        self.assertEqual(window.project.bodies[0].mass, body.mass * 2)
        window.mode.setCurrentIndex(0)
        self.assertEqual(window.display_project.bodies[0].mass, body.mass * 2)

    def test_workspace_persistence_and_small_window(self):
        with tempfile.TemporaryDirectory() as directory:
            settings = QSettings(
                str(Path(directory) / "ui.ini"), QSettings.Format.IniFormat
            )
            window = self.make_window(settings)
            window.resize(1280, 800)
            window.set_theme("light")
            window.docks["inspector"].hide()
            window.save_workspace()
            window.close()
            restored = self.make_window(settings)
            self.assertEqual(restored.theme_name, "light")
            self.assertFalse(restored.docks["inspector"].isVisible())
            self.assertLessEqual(restored.width(), 1280)
            restored.reset_workspace()
            self.assertTrue(restored.docks["inspector"].isVisible())
            self.assertGreaterEqual(restored.viewport.width(), 400)
            self.assertGreaterEqual(restored.viewport.height(), 300)
            restored.close()

    def test_save_discard_cancel_protect_document(self):
        from vinkulum_studio.editor import EditorWindow

        window = self.make_window()
        window.new_project()
        window.add_body("box")
        with patch.object(
            QMessageBox, "question", return_value=QMessageBox.StandardButton.Save
        ):
            with patch.object(window, "save_dialog", return_value=False):
                self.assertFalse(EditorWindow._discard_allowed(window))
            with patch.object(window, "save_dialog", return_value=True):
                self.assertTrue(EditorWindow._discard_allowed(window))
        with patch.object(
            QMessageBox, "question", return_value=QMessageBox.StandardButton.Cancel
        ):
            self.assertFalse(EditorWindow._discard_allowed(window))

    def test_curve_zoom_and_native_samples_remain_exact(self):
        from vinkulum_studio.series_view import SeriesView

        curve = SeriesView()
        curve.resize(900, 260)
        times = np.linspace(0, 1, 20001)
        values = np.zeros_like(times)
        values[10001] = (
            1e6  # Narrow feature must not disappear from the displayed range.
        )
        curve.set_series(times, values, "Impulsion", "N")
        curve.show()
        self.addCleanup(curve.close)
        QTest.qWait(20)
        self.assertGreater(curve._geometry()[2], 1e6)
        curve._set_domain(0.4999, 0.5002)
        self.assertGreater(curve._geometry()[2], 1e6)
        selected = []
        curve.sample_selected.connect(selected.append)
        curve.index = 10000
        curve.setFocus()
        QTest.keyClick(curve, Qt.Key.Key_Right)
        self.assertEqual(selected, [10001])
        self.assertIs(curve.values, values)
        self.assertEqual(values[10001], 1e6)
        curve.reset_view()
        self.assertEqual(curve.domain, (0.0, 1.0))


class RunArchiveContracts(unittest.TestCase):
    def test_memory_budget_evicts_oldest_without_copying_arrays(self):
        from types import SimpleNamespace

        from vinkulum_studio.run_archive import ARRAYS, RunArchive, result_bytes

        array = np.zeros(10)

        def result(identifier):
            return SimpleNamespace(
                run_id=identifier, **{name: array for name in ARRAYS}
            )

        first, second, third = map(result, ("first", "second", "third"))
        archive = RunArchive(max_runs=8, max_bytes=2 * result_bytes(first))
        archive.add(first)
        archive.add(second)
        self.assertEqual(archive.add(third), ["first"])
        self.assertIsNone(archive.find("first"))
        self.assertIs(archive.find("second"), second)
        self.assertIs(archive.find("third").time, array)
        with self.assertRaises(ValueError):
            RunArchive(max_bytes=1).add(first)
