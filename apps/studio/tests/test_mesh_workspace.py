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
from PySide6.QtGui import QImage
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog
from shiboken6 import isValid
from test_meshing import condition, solid
from vinkulum_studio.calculix import load_study, study_document
from vinkulum_studio.mesh_binding import MeshBinding
from vinkulum_studio.mesh_view import SURFACE_COLORS
from vinkulum_studio.mesh_window import (
    ConditionDialog,
    ConditionEdit,
    MeshWindow,
    mesh_presentation,
)


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires a desktop OpenGL context"
)
class MeshWorkspace(unittest.TestCase):
    def test_rendered_symbols_remain_distinct_from_colored_faces(self):
        w = self.window()
        v = w.viewport
        w.edges.setChecked(False)
        camera = v.renderer.GetActiveCamera()
        camera.SetParallelProjection(True)
        camera.SetParallelScale(0.8)
        camera.SetPosition(0.4, -3, 2)
        camera.SetFocalPoint(0.3, 0, 0.3)
        camera.SetViewUp(0, 0, 1)

        def luminance(color):
            channels = (color.redF(), color.greenF(), color.blueF())
            linear = [
                c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
                for c in channels
            ]
            return sum(
                c * weight for c, weight in zip(linear, (0.2126, 0.7152, 0.0722))
            )

        for row in (
            condition("pressure", (2,), values=(10.0,)),
            condition("total_force", (2,), values=(10.0, 0.0, 0.0)),
            condition("support", (2,), axes=(1, 3)),
        ):
            with self.subTest(kind=row.kind):
                w._set_conditions((row,))
                pictures = []
                for enabled in (False, True):
                    w.symbols.setChecked(enabled)
                    path = self.root / f"contrast-{row.kind}-{enabled}.png"
                    v.screenshot(path)
                    pictures.append(QImage(str(path)))
                point = v._frames[2][0].point
                v.renderer.SetWorldPoint(*point, 1)
                v.renderer.WorldToDisplay()
                x, y, _ = v.renderer.GetDisplayPoint()
                x, y = round(x), pictures[0].height() - 1 - round(y)
                radius = round(45 * pictures[0].width() / v.view.width())
                distinct = 0
                for a in range(x - radius, x + radius):
                    for b in range(y - radius, y + radius):
                        first, second = (
                            luminance(im.pixelColor(a, b)) for im in pictures
                        )
                        if (max(first, second) + 0.05) / (
                            min(first, second) + 0.05
                        ) >= 1.8:
                            distinct += 1
                self.assertGreater(
                    distinct, 12, "Symbols are lost against their face colour"
                )

    def test_rendered_force_arrow_changes_world_axis_without_pickable_geometry(self):
        w = self.window()
        v = w.viewport
        v.actor.VisibilityOff()
        v.edge_actor.VisibilityOff()
        camera = v.renderer.GetActiveCamera()
        camera.SetParallelProjection(True)
        camera.SetParallelScale(0.8)
        camera.SetPosition(0.4, 0.4, 3)
        camera.SetFocalPoint(0.4, 0.4, 0)
        camera.SetViewUp(0, 1, 0)
        camera.SetClippingRange(0.01, 10)
        for axis, values in ((0, (5.0, 0.0, 0.0)), (1, (0.0, 5.0, 0.0))):
            w._set_conditions((condition("total_force", (2,), values=values),))
            v.render()
            group = v.glyph_layer.groups["total_force"]
            v.renderer.SetWorldPoint(*group["positions"][0], 1)
            v.renderer.WorldToDisplay()
            x, y, _ = v.renderer.GetDisplayPoint()
            path = self.root / f"force-axis-{axis}.png"
            v.screenshot(path)
            image = QImage(str(path))
            self.assertFalse(image.isNull())
            ratio = image.width() / v.view.width()
            x, y = round(x), image.height() - 1 - round(y)
            radius = round(45 * ratio)
            colored = []
            for a in range(max(0, x - radius), min(image.width(), x + radius)):
                for b in range(max(0, y - radius), min(image.height(), y + radius)):
                    r, g, blue, _ = image.pixelColor(a, b).getRgb()
                    if r > 70 and r > 1.15 * g and r > 1.1 * blue:
                        colored.append((a, b))
            self.assertGreater(len(colored), 10)
            widths = [
                (max(p[i] for p in colored) - min(p[i] for p in colored) + 1) / ratio
                for i in (0, 1)
            ]
            self.assertAlmostEqual(widths[axis], 30, delta=3)
            self.assertLess(widths[1 - axis], 12)
            self.assertFalse(group["actor"].GetPickable())

    def test_boundary_symbols_follow_multiplier_and_preserve_physical_conditions(self):
        w = self.window()
        conditions = (
            condition("support", (1,), axes=(1, 3)),
            condition("pressure", (2,), values=(20.0,)),
            condition("total_force", (3,), values=(3.0, 4.0, 0.0)),
        )
        w._set_conditions(conditions)
        layer = w.viewport.glyph_layer
        before = {
            kind: [g.direction for g in group["rows"]]
            for kind, group in layer.groups.items()
        }
        w.load_factor.setText("-2")
        for kind, group in layer.groups.items():
            expected = (
                before[kind]
                if kind == "support"
                else [tuple(-v for v in d) for d in before[kind]]
            )
            self.assertEqual([g.direction for g in group["rows"]], expected)
            self.assertFalse(group["actor"].GetPickable())
        self.assertIn("Pressure -40 Pa", w.condition_list.item(1).text())
        self.assertEqual(w.conditions, conditions)
        w.load_factor.setText("1.000000000000001")
        displayed, tooltip = w.load_factor.text(), w.condition_list.item(1).toolTip()
        w.load_factor.setText("1.000000000000002")
        self.assertEqual(w.load_factor.text(), displayed)
        self.assertNotEqual(w.condition_list.item(1).toolTip(), tooltip)
        self.assertIn(
            repr(20.0 * w.load_factor.value()), w.condition_list.item(1).toolTip()
        )
        for value in ("invalid", "nan", "inf"):
            w.load_factor.setText(value)
            self.assertFalse(layer.groups["pressure"]["rows"])
            self.assertFalse(layer.groups["total_force"]["rows"])
            self.assertTrue(layer.groups["support"]["rows"])
            self.assertFalse(w.save.isEnabled() or w.open_calculix.isEnabled())
            self.assertIn("Correct", w.symbol_hint.text())
        w.load_factor.setText("0")
        self.assertFalse(
            layer.groups["pressure"]["rows"] or layer.groups["total_force"]["rows"]
        )
        self.assertTrue(w.save.isEnabled())
        w.load_factor.setText("1")
        w.symbols.setFocus()
        QTest.keyClick(w.symbols, Qt.Key.Key_Space)
        self.assertFalse(any(g["actor"].GetVisibility() for g in layer.groups.values()))
        QTest.keyClick(w.symbols, Qt.Key.Key_Space)
        self.assertTrue(all(g["actor"].GetVisibility() for g in layer.groups.values()))
        visible = []
        for identifier, frames in w.viewport._frames.items():
            for frame in frames:
                point = self.screen_point(w.viewport, frame.point)
                if w.viewport.pick_surface(point.x(), point.y()) == identifier:
                    visible.append((point, identifier))
        self.assertTrue(visible)
        point, identifier = visible[0]
        QTest.mouseClick(w.viewport.view, Qt.MouseButton.LeftButton, pos=point)
        self.assertEqual(w.selected_faces, {identifier})
        self.assertEqual(w.conditions, conditions)

    def test_direction_symbol_projection_keeps_logical_size_during_zoom_and_resize(
        self,
    ):
        import numpy as np

        w = self.window()
        w._set_conditions((condition("total_force", (2,), values=(3.0, 4.0, 0.0)),))
        v = w.viewport
        camera = v.renderer.GetActiveCamera()
        camera.SetPosition(0.4, 0.4, 3)
        camera.SetFocalPoint(0.4, 0.4, 0)
        camera.SetViewUp(0, 1, 0)

        def projected(point):
            v.renderer.SetWorldPoint(*point, 1)
            v.renderer.WorldToDisplay()
            x, y, _ = v.renderer.GetDisplayPoint()
            width, height = v.view.GetRenderWindow().GetSize()
            return np.array((x * v.view.width() / width, y * v.view.height() / height))

        for parallel in (True, False):
            camera.SetParallelProjection(parallel)
            camera.SetParallelScale(0.8)
            camera.SetViewAngle(35)
            for zoom, size in ((1, (1280, 850)), (2, (1100, 740)), (0.5, (1400, 950))):
                camera.Zoom(zoom)
                w.resize(*size)
                QTest.qWait(30)
                v.render()
                group = v.glyph_layer.groups["total_force"]
                data = group["mapper"].GetInput()
                start = np.array(data.GetPoint(0))
                length = data.GetPointData().GetArray("symbol_scale").GetValue(0)
                direction = np.array(
                    data.GetPointData().GetArray("direction").GetTuple3(0)
                )
                self.assertAlmostEqual(
                    np.linalg.norm(
                        projected(start + length * direction) - projected(start)
                    ),
                    30,
                    delta=0.1,
                )
        # The fitting bounds belong to the captured geometry, independently of
        # symbol size or the magnitude of the total-force components.
        v.camera("iso")
        before = camera.GetPosition(), camera.GetFocalPoint(), camera.GetParallelScale()
        w._set_conditions((condition("total_force", (2,), values=(3e300, 4e300, 0.0)),))
        self.assertEqual(
            before,
            (camera.GetPosition(), camera.GetFocalPoint(), camera.GetParallelScale()),
        )
        # Compare fits from the same canonical direction. Repeated ResetCamera
        # normalisation alone can move a component by one ULP on ARM64.
        v.camera("iso")
        self.assertEqual(
            before,
            (camera.GetPosition(), camera.GetFocalPoint(), camera.GetParallelScale()),
        )

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
        window._present_mesh(mesh_presentation(s, self.root))
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
