"""User-visible selection, context actions and continuous camera navigation."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from PySide6.QtCore import QPoint, QPointF, QSettings, Qt
from PySide6.QtGui import QNativeGestureEvent, QPointingDevice, QWheelEvent
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QMenu


@unittest.skipUnless(os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires OpenGL")
class PointerNavigation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        from vinkulum_studio.editor import EditorWindow

        self.errors = []
        self.hook = patch(
            "sys.excepthook", side_effect=lambda *args: self.errors.append(args)
        )
        self.hook.start()
        self.directory = tempfile.TemporaryDirectory()
        self.window = EditorWindow(
            QSettings(
                str(Path(self.directory.name) / "ui.ini"), QSettings.Format.IniFormat
            )
        )
        self.window._discard_allowed = lambda: True
        self.window.new_project()
        self.window.add_body("box")
        self.window.show()
        self.window.activateWindow()
        self.assertTrue(QTest.qWaitForWindowActive(self.window))
        QTest.qWait(60)
        self.view = self.window.viewport
        self.body = self.window.project.bodies[0]

    def tearDown(self):
        menu = self.app.activePopupWidget()
        if menu:
            menu.close()
        self.window.close()
        QTest.qWait(30)
        self.hook.stop()
        self.directory.cleanup()
        self.assertFalse(self.errors, self.errors)

    def body_point(self):
        self.view.renderer.SetWorldPoint(*self.body.position, 1.0)
        self.view.renderer.WorldToDisplay()
        x, y, _ = self.view.renderer.GetDisplayPoint()
        ratio = self.view.view._getPixelRatio()
        return QPoint(round(x / ratio), round(self.view.view.height() - 1 - y / ratio))

    def test_right_click_selects_target_and_exposes_actions_without_moving_camera(self):
        self.window.select_object(None)
        camera = self.view.renderer.GetActiveCamera()
        before = camera.GetPosition(), camera.GetFocalPoint()
        project = self.window.project
        QTest.mouseClick(
            self.view.view, Qt.MouseButton.RightButton, pos=self.body_point()
        )
        QTest.qWait(30)
        menu = self.app.activePopupWidget()
        self.assertIsInstance(menu, QMenu)
        self.assertEqual(self.window.selection, self.body.id)
        self.assertEqual((camera.GetPosition(), camera.GetFocalPoint()), before)
        self.assertIs(self.window.project, project)
        labels = {action.text() for action in menu.actions()}
        self.assertIn("Edit properties…", labels)
        self.assertIn("Duplicate selection", labels)
        self.assertIn("Fit all visible objects", labels)
        menu.close()
        QTest.qWait(20)
        self.assertIsNone(self.app.activePopupWidget())

    def test_drag_orbits_without_selecting_and_click_selects_only_on_release(self):
        self.window.select_object(None)
        point = self.body_point()
        before = self.view.renderer.GetActiveCamera().GetPosition()
        QTest.mousePress(self.view.view, Qt.MouseButton.LeftButton, pos=point)
        self.assertIsNone(self.window.selection)
        QTest.mouseMove(self.view.view, point + QPoint(80, 30), 20)
        QTest.mouseRelease(
            self.view.view, Qt.MouseButton.LeftButton, pos=point + QPoint(80, 30)
        )
        self.assertIsNone(self.window.selection)
        self.assertNotEqual(self.view.renderer.GetActiveCamera().GetPosition(), before)
        point = self.body_point()
        QTest.mousePress(self.view.view, Qt.MouseButton.LeftButton, pos=point)
        self.assertIsNone(self.window.selection)
        QTest.mouseRelease(self.view.view, Qt.MouseButton.LeftButton, pos=point)
        self.assertEqual(self.window.selection, self.body.id)

    def test_selection_preserves_browser_state_and_tabs_preserve_camera(self):
        joints = self.window.tree.topLevelItem(1)
        joints.setExpanded(False)
        self.window.select_object(None)
        self.window.select_object(self.body.id)
        self.assertFalse(joints.isExpanded())
        camera = self.view.renderer.GetActiveCamera()
        camera.Azimuth(22)
        camera.Elevation(13)
        before = camera.GetPosition(), camera.GetFocalPoint(), camera.GetParallelScale()
        self.window.workspace_tabs.setCurrentIndex(1)
        QTest.qWait(40)
        self.window.workspace_tabs.setCurrentIndex(0)
        QTest.qWait(40)
        self.assertEqual(
            (camera.GetPosition(), camera.GetFocalPoint(), camera.GetParallelScale()),
            before,
        )
        direction = camera.GetDirectionOfProjection()
        self.window.commands["fit"].trigger()
        np.testing.assert_allclose(
            camera.GetDirectionOfProjection(), direction, atol=1e-14
        )

    def test_fractional_wheel_and_native_pinch_zoom_without_changing_document(self):
        camera = self.view.renderer.GetActiveCamera()
        before = camera.GetDistance()
        project = self.window.project
        point = QPointF(self.body_point())
        wheel = QWheelEvent(
            point,
            point,
            QPoint(),
            QPoint(0, 15),
            Qt.MouseButton.NoButton,
            Qt.KeyboardModifier.NoModifier,
            Qt.ScrollPhase.ScrollUpdate,
            False,
        )
        self.app.sendEvent(self.view.view, wheel)
        self.assertLess(camera.GetDistance(), before)
        before = camera.GetDistance()
        pinch = QNativeGestureEvent(
            Qt.NativeGestureType.ZoomNativeGesture,
            QPointingDevice.primaryPointingDevice(),
            2,
            point,
            point,
            point,
            0.1,
            QPointF(),
        )
        self.app.sendEvent(self.view.view, pinch)
        self.assertLess(camera.GetDistance(), before)
        self.assertIs(self.window.project, project)

    def test_mac_control_click_and_pixel_scroll_use_context_and_pan(self):
        point = self.body_point()
        camera = self.view.renderer.GetActiveCamera()
        before = camera.GetPosition(), camera.GetFocalPoint()
        with patch("vinkulum_studio.render_interactor.sys.platform", "darwin"):
            QTest.mouseClick(
                self.view.view,
                Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.MetaModifier,
                point,
            )
            QTest.qWait(20)
            menu = self.app.activePopupWidget()
            self.assertIsInstance(menu, QMenu)
            self.assertEqual((camera.GetPosition(), camera.GetFocalPoint()), before)
            menu.close()
            QTest.qWait(20)
            distance = camera.GetDistance()
            scroll = QWheelEvent(
                QPointF(point),
                QPointF(point),
                QPoint(12, 6),
                QPoint(),
                Qt.MouseButton.NoButton,
                Qt.KeyboardModifier.NoModifier,
                Qt.ScrollPhase.ScrollUpdate,
                False,
            )
            self.app.sendEvent(self.view.view, scroll)
            self.assertNotEqual(camera.GetFocalPoint(), before[1])
            self.assertAlmostEqual(camera.GetDistance(), distance, places=12)
