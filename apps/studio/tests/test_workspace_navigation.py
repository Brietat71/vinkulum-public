"""Returning between real workspaces preserves the user's pending study."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QToolButton


@unittest.skipUnless(os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires OpenGL")
class WorkspaceNavigation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_visible_workspace_menu_and_return_keep_the_same_edited_study(self):
        from vinkulum_studio.editor import EditorWindow

        with (
            tempfile.TemporaryDirectory() as directory,
            patch("sys.excepthook") as errors,
        ):
            window = EditorWindow(
                QSettings(str(Path(directory) / "ui.ini"), QSettings.Format.IniFormat)
            )
            window._discard_allowed = lambda: True
            window.show()
            window.activateWindow()
            self.assertTrue(QTest.qWaitForWindowActive(window))
            try:
                buttons = [
                    b
                    for b in window.findChildren(QToolButton)
                    if b.accessibleName() == "Open an analysis workspace"
                ]
                self.assertEqual(len(buttons), 1)
                menu = buttons[0].menu()
                self.assertIn(window.commands["static_study"], menu.actions())
                window.commands["static_study"].trigger()
                study = window._static_window
                study._discard_allowed = lambda: True
                self.assertTrue(QTest.qWaitForWindowActive(study))
                fields = study.findChildren(QToolButton)
                back = study.findChild(QAction, "back_to_workspace")
                self.assertIsNotNone(back)
                button = next(b for b in fields if b.defaultAction() is back)
                camera = study.viewport.renderer.GetActiveCamera()
                camera.Azimuth(25)
                original_camera = camera.GetPosition(), camera.GetFocalPoint()
                captured = study.study
                study.young.setFocus()
                study.young.selectAll()
                QTest.keyClicks(study.young, "123e9")
                self.assertEqual(study.young.text(), "123e9")
                QTest.mouseClick(button, Qt.MouseButton.LeftButton)
                QTest.qWait(30)
                self.assertFalse(study.isVisible())
                self.assertIs(window._static_window, study)
                window.commands["static_study"].trigger()
                self.assertTrue(QTest.qWaitForWindowActive(study))
                self.assertIs(window._static_window, study)
                self.assertIs(study.study, captured)
                self.assertEqual(study.young.text(), "123e9")
                self.assertEqual(
                    (camera.GetPosition(), camera.GetFocalPoint()), original_camera
                )
            finally:
                window.close()
                QTest.qWait(40)
            errors.assert_not_called()
