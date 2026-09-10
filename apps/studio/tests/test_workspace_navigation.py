"""A single project window retains edited analyses across real browser clicks."""

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog, QMainWindow, QMenu, QToolButton


@unittest.skipUnless(os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires OpenGL")
class WorkspaceNavigation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        from vinkulum_studio.editor import EditorWindow

        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.errors = patch("sys.excepthook").start()
        self.addCleanup(patch.stopall)
        self.window = EditorWindow(
            QSettings(
                str(Path(self.directory.name) / "ui.ini"), QSettings.Format.IniFormat
            )
        )
        self.window._discard_allowed = lambda: True
        self.window.resize(1280, 800)
        self.window.show()
        self.window.activateWindow()
        self.assertTrue(QTest.qWaitForWindowActive(self.window))
        QTest.qWait(40)

    def tearDown(self):
        menu = self.app.activePopupWidget()
        if menu:
            menu.close()
        for page, _ in self.window.analysis_pages.values():
            page._discard_allowed = lambda: True
        self.window.close()
        QTest.qWait(40)
        self.errors.assert_not_called()

    def click_browser(self, identifier):
        tree = self.window.tree
        for g in range(tree.topLevelItemCount()):
            group = tree.topLevelItem(g)
            for i in range(group.childCount()):
                item = group.child(i)
                if item.data(0, Qt.ItemDataRole.UserRole) == identifier:
                    group.setExpanded(True)
                    tree.scrollToItem(item)
                    QTest.mouseClick(
                        tree.viewport(),
                        Qt.MouseButton.LeftButton,
                        pos=tree.visualItemRect(item).center(),
                    )
                    QTest.qWait(30)
                    return
        self.fail(f"Missing browser entry: {identifier}")

    def assert_integrated(self, page):
        w = self.window
        self.assertFalse(page.isWindow())
        self.assertIs(page.window(), w)
        self.assertIs(w.project_pages.currentWidget(), page)
        self.assertTrue(w.tree.isVisible())
        self.assertTrue(w.rect().contains(w.tree.mapTo(w, w.tree.rect().center())))
        self.assertEqual(
            [
                p
                for p in self.app.topLevelWidgets()
                if isinstance(p, QMainWindow) and p.isVisible()
            ],
            [w],
        )
        self.assertIsNone(page.findChild(QAction, "back_to_workspace"))
        self.assertEqual(w.width(), 1280)

    def test_browser_keeps_same_edited_study_and_camera_in_one_window(self):
        w = self.window
        buttons = [
            b
            for b in w.findChildren(QToolButton)
            if b.accessibleName() == "Add or select an analysis"
        ]
        self.assertEqual(len(buttons), 1)
        self.assertIn(w.commands["static_study"], buttons[0].menu().actions())
        w.commands["static_study"].trigger()
        study = w._static_window
        self.assert_integrated(study)
        camera = study.viewport.renderer.GetActiveCamera()
        camera.Azimuth(25)
        before = camera.GetPosition(), camera.GetFocalPoint()
        captured = study.study
        study.young.setFocus()
        study.young.selectAll()
        QTest.keyClicks(study.young, "123e9")
        self.click_browser(w.project.bodies[0].id)
        self.assertIsNone(w.active_analysis)
        self.assertFalse(study.isVisible())
        self.click_browser("analysis:static")
        self.assert_integrated(study)
        self.assertIs(study.study, captured)
        self.assertEqual(study.young.value(), 123e9)
        self.assertEqual((camera.GetPosition(), camera.GetFocalPoint()), before)
        w.commands["fit_all"].trigger()
        self.assertEqual(w.viewport._project, w.project)

    def test_save_targets_displayed_analysis_and_background_refresh_keeps_context(self):
        from vinkulum_studio.calculix import load_study

        w = self.window
        w.commands["static_study"].trigger()
        page = w._static_window
        captured = w.project
        target = str(Path(self.directory.name) / "study.ccx.json")
        with patch.object(QFileDialog, "getSaveFileName", return_value=(target, "")):
            w.commands["save"].trigger()
        self.assertEqual(load_study(target), page.study)
        self.assertIs(w.project, captured)
        self.assertIsNone(w._path)
        w._refresh()  # Native background refresh must not steal this view.
        self.assert_integrated(page)
        self.assertEqual(w.active_analysis, "static")
        self.assertFalse(w.commands["delete"].isEnabled())
        self.click_browser("analysis:native")
        self.assertTrue(w.docks["analysis"].isVisible())
        self.assertIsNone(w.active_analysis)
        self.assertEqual(w.commands["save"].text(), "Save project…")

    def test_cad_to_statics_stays_embedded_and_background_completion_does_not_steal_view(
        self,
    ):
        from vinkulum_studio.static_window import tension_example

        w = self.window
        root = Path(__file__).resolve().parents[3]
        w.load(root / "examples/studio/platine-parametrique.vinkulum.json")
        w.select_object(w.project.bodies[0].id)
        w.commands["cad_study"].trigger()
        mesh = w._mesh_window
        self.assert_integrated(mesh)
        self.click_browser(w.project.bodies[0].id)
        mesh._study_ready(tension_example(), "open")
        self.assertIsNone(w.active_analysis)
        self.assertIs(w.project_pages.currentWidget(), w.model_page)
        self.click_browser("analysis:cad_static")
        self.assert_integrated(mesh._static_window)
        self.assertEqual(mesh._static_window.study, tension_example())
        self.click_browser("analysis:cad")
        mesh._study_ready(tension_example(), "open")
        self.assert_integrated(mesh._static_window)

    def test_hidden_model_panel_and_contextual_undo_survive_analysis_navigation(self):
        w = self.window
        w.docks["inspector"].hide()
        w.add_body("box")
        self.assertTrue(w.commands["undo"].isEnabled())
        project = w.project
        w.commands["static_study"].trigger()
        self.assertFalse(w.commands["undo"].isEnabled())
        w.set_theme("light")
        self.assertEqual(
            w._static_window.palette().window().color(), w.palette().window().color()
        )
        self.click_browser(project.bodies[-1].id)
        self.assertFalse(w.docks["inspector"].isVisible())
        self.assertTrue(w.commands["undo"].isEnabled())
        self.assertTrue(w.commands["box"].isEnabled())
        w.commands["undo"].trigger()
        self.assertEqual(len(w.project.bodies), len(project.bodies) - 1)

    def test_selecting_a_mesh_face_exposes_boundary_controls(self):
        from test_mesh_workspace import MeshWorkspace
        from test_meshing import solid
        from vinkulum_studio.mesh_window import mesh_presentation

        w = self.window
        capture = solid(curve=0.05)
        w._commit(w.project.replace_object(capture.request.body))
        w.select_object(capture.request.body.id)
        w.open_cad_study()
        mesh = w._mesh_window
        mesh._present_mesh(mesh_presentation(capture, Path(self.directory.name)))
        mesh.viewport.camera("top")
        QTest.qWait(40)
        point = MeshWorkspace.screen_point(mesh.viewport, (1 / 3, 1 / 3, 1 / 3))
        QTest.mouseClick(mesh.viewport.view, Qt.MouseButton.LeftButton, pos=point)
        QTest.qWait(40)
        self.assertEqual(mesh.selected_faces, {3})
        self.assert_integrated(mesh)
        button = mesh.add_buttons["support"]
        self.assertTrue(button.isVisible())
        self.assertTrue(button.isEnabled())
        self.assertTrue(
            mesh.faces_dock.visibleRegion().contains(
                button.mapTo(mesh.faces_dock, button.rect().center())
            )
        )

    def test_context_close_honours_pending_analysis_guard(self):
        from PySide6.QtWidgets import QMessageBox

        w = self.window
        w.commands["static_study"].trigger()
        page = w._static_window
        page.young.setText("123e9")
        w.analysis_menu("static", w.tree.mapToGlobal(w.tree.rect().center()))
        QTest.qWait(20)
        menu = self.app.activePopupWidget()
        self.assertIsInstance(menu, QMenu)
        with patch.object(
            QMessageBox, "question", return_value=QMessageBox.StandardButton.Cancel
        ):
            next(a for a in menu.actions() if a.text() == "Close analysis…").trigger()
        self.assertIs(w._static_window, page)
        self.assert_integrated(page)
        menu.close()
        page._discard_allowed = lambda: True
        self.assertTrue(page.close())
        QTest.qWait(30)
        self.assertIsNone(w._static_window)
        self.assertNotIn("static", w.analysis_pages)
        self.assertIs(w.project_pages.currentWidget(), w.model_page)
