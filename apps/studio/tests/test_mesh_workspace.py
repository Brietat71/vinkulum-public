"""Real Qt/VTK face picking, condition editing and captured-study handoff."""

import json
import os
import tempfile
import threading
import time
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QPoint, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog
from shiboken6 import isValid
from test_meshing import condition, solid
from vinkulum_studio.calculix import load_study, study_document
from vinkulum_studio.mesh_binding import MeshBinding, mesh_measurements
from vinkulum_studio.mesh_view import SURFACE_COLORS
from vinkulum_studio.mesh_window import (
    ConditionDialog,
    ConditionEdit,
    MeshWindow,
    surface_areas,
)


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires a desktop OpenGL context"
)
class MeshWorkspace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.windows = []
        self.catch = patch("sys.excepthook")
        self.exceptions = self.catch.start()
        self.addCleanup(self.catch.stop)
        self.addCleanup(self.exceptions.assert_not_called)
        self.addCleanup(self.close_windows)

    def close_windows(self):
        for window in reversed(self.windows):
            if isValid(window):
                window.cancel_work()
                self.wait_until(lambda window=window: not window.busy)
                with patch.object(window, "_discard_allowed", return_value=True):
                    window.close()
        QTest.qWait(20)

    def wait_until(self, condition, seconds=10):
        deadline = time.monotonic() + seconds
        while not condition() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(condition(), "Workspace did not reach the expected state.")

    def window(self, *, capture=None):
        s = solid(curve=0.05)
        window = MeshWindow(s.request.body, capture=capture)
        self.windows.append(window)
        window.show()
        window._present_mesh((s, self.root, surface_areas(s), mesh_measurements(s)))
        QTest.qWait(60)
        return window

    @staticmethod
    def screen_point(view, point):
        view.renderer.SetWorldPoint(*point, 1.0)
        view.renderer.WorldToDisplay()
        x, y, _ = view.renderer.GetDisplayPoint()
        width, height = view.view.GetRenderWindow().GetSize()
        return QPoint(
            round(x * view.view.width() / width),
            round((height - 1 - y) * view.view.height() / height),
        )

    def test_real_face_pick_toggle_orbit_and_true_quadratic_edge_display(self):
        w = self.window()
        v = w.viewport
        v.camera("top")
        QTest.qWait(30)
        # Sloping face opposite node 1 is exposed from above.
        point = self.screen_point(v, (1 / 3, 1 / 3, 1 / 3))
        self.assertEqual(v.pick_surface(point.x(), point.y()), 3)
        QTest.mouseClick(v.view, Qt.MouseButton.LeftButton, pos=point)
        self.assertEqual(w.selected_faces, {3})
        colors = v.mesh.GetCellData().GetScalars()
        for i, face in enumerate(v._cell_faces):
            if face == 3:
                self.assertEqual(tuple(colors.GetTuple3(i)), SURFACE_COLORS["selected"])
        QTest.mouseClick(
            v.view,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.ControlModifier,
            point,
        )
        self.assertEqual(w.selected_faces, set())
        w.set_selected_faces((3,))
        QTest.mousePress(v.view, Qt.MouseButton.LeftButton, pos=point)
        QTest.mouseMove(v.view, point + QPoint(60, 30), 30)
        QTest.mouseRelease(
            v.view, Qt.MouseButton.LeftButton, pos=point + QPoint(60, 30)
        )
        self.assertEqual(w.selected_faces, {3})
        # One tetrahedron has six edges, irrespective of display tessellation.
        self.assertEqual(v.edge_mesh.GetNumberOfCells(), 6)
        self.assertEqual(v.edge_mesh.GetNumberOfPoints(), 30)
        self.assertFalse(v.actor.GetProperty().GetEdgeVisibility())
        w.edges.setChecked(False)
        self.assertFalse(v.edge_actor.GetVisibility())
        w.filter.setText("Face 4")
        self.assertEqual(w.selected_faces, {3})
        w._face_items[4].setSelected(True)
        self.assertEqual(w.selected_faces, {3, 4})

    def test_condition_dialog_validates_and_keeps_identity_and_world_components(self):
        w = self.window()
        old = condition("total_force", (1, 3), values=(1.2345678901234567, -2.0, 3.0))
        dialog = ConditionDialog("total_force", (2,), w, old)
        dialog.show()
        QTest.qWait(10)
        dialog.accept()
        self.assertEqual(dialog.condition.id, old.id)
        self.assertEqual(dialog.condition.surfaces, (2,))
        self.assertEqual(dialog.condition.values, old.values)
        support = ConditionDialog("support", (1,), w)
        for field in support.axes:
            field.setChecked(False)
        support.accept()
        self.assertIsNone(support.condition)
        self.assertIn("world axes", support.error.text())
        support.axes[2].setChecked(True)
        support.accept()
        self.assertEqual(support.condition.axes, (3,))

    def test_undo_save_reopen_multiplier_and_handoff_preserve_captured_intent(self):
        w = self.window()
        support = condition("support", (1,), axes=(1, 2, 3))
        pressure = condition("pressure", (3,), values=(120.0,))
        w.undo.push(ConditionEdit(w, (support, pressure), "Add conditions"))
        w.undo.undo()
        self.assertEqual(w.conditions, ())
        w.undo.redo()
        self.assertEqual(w.conditions, (support, pressure))
        w.load_factor.setText("-1.2345678901234567")
        expected = MeshBinding(w.solid, w.conditions, w.load_factor.value()).study(
            210e9, 0.3
        )
        path = self.root / "saved.ccx.json"
        with patch.object(
            QFileDialog, "getSaveFileName", return_value=(str(path), "JSON")
        ):
            QTest.mouseClick(w.save, Qt.MouseButton.LeftButton)
            self.wait_until(lambda: not w.busy)
        self.assertEqual(load_study(path), expected)
        w.undo.push(ConditionEdit(w, (support,), "Remove pressure"))
        with (
            patch.object(w, "_discard_allowed", return_value=True),
            patch.object(
                QFileDialog, "getOpenFileName", return_value=(str(path), "JSON")
            ),
        ):
            w.open_study()
            self.wait_until(lambda: not w.busy)
        self.assertEqual(w.conditions, expected.mesh_binding.conditions)
        self.assertEqual(w.load_factor.value(), expected.mesh_binding.load_factor)
        QTest.mouseClick(w.open_calculix, Qt.MouseButton.LeftButton)
        self.wait_until(lambda: not w.busy)
        self.assertIsNotNone(w._static_window, w.status.text())
        self.assertEqual(w._static_window.study, expected)
        self.assertIn("face conditions included", w._static_window.limits.text())

    def test_invalid_archive_is_transactional_and_new_source_clears_conditions(self):
        changed = replace(
            solid().request.body, name="Changed source", position=(3.0, 4.0, 5.0)
        )
        w = self.window(capture=lambda: changed)
        support = condition("support", (1,), axes=(1, 2, 3))
        w.undo.push(ConditionEdit(w, (support,), "Support"))
        before = w._state()
        path = self.root / "bad.ccx.json"
        data = study_document(MeshBinding(w.solid, w.conditions).study(210e9, 0.3))
        data["study"]["mesh_binding"]["solid"]["input_sha256"] = "f" * 64
        path.write_text(json.dumps(data))
        with (
            patch.object(w, "_discard_allowed", return_value=True),
            patch.object(
                QFileDialog, "getOpenFileName", return_value=(str(path), "JSON")
            ),
        ):
            w.open_study()
            self.wait_until(lambda: not w.busy)
        self.assertEqual(w._state(), before)
        self.assertIn("fingerprint", w.status.text())
        with patch.object(w, "_discard_allowed", return_value=True):
            w.capture_selected()
        self.assertEqual(w.source, changed)
        self.assertIsNone(w.solid)
        self.assertFalse(w.conditions)
        self.assertFalse(w.undo.canUndo())
        self.assertFalse(w.open_calculix.isEnabled())

    def test_close_waits_for_background_checks_and_does_not_publish_after_cancel(self):
        w = self.window()
        started, finish = threading.Event(), threading.Event()
        self.addCleanup(finish.set)
        published, closed = [], []
        w.closed.connect(lambda: closed.append(True))

        def slow():
            started.set()
            finish.wait(5)
            return "must not be presented"

        w._background(slow, published.append, "Checking fixture")
        self.wait_until(started.is_set)
        viewport = w.viewport
        self.assertFalse(w.close())
        self.assertFalse(viewport._closed)
        self.assertTrue(w._close_pending)
        finish.set()
        self.wait_until(lambda: bool(closed))
        self.assertFalse(published)
        self.assertTrue(viewport._closed)


if __name__ == "__main__":
    unittest.main()
