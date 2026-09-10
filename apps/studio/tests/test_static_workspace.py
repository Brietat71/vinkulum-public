"""Real Qt/ccx lifecycle and captured finite-element results in the workbench."""

import csv
import os
import shutil
import sys
import tempfile
import time
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QImage
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog
from vinkulum_studio.calculix import load_study, run_static
from vinkulum_studio.static_controller import StaticController
from vinkulum_studio.static_window import tension_example


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires a desktop OpenGL context"
)
class StaticWorkspace(unittest.TestCase):
    @unittest.skipUnless(shutil.which("ccx"), "Requires CalculiX")
    def test_quadratic_curved_mesh_display_and_integration_point_table(self):
        from test_tetrahedra import tetra_specimen
        from vinkulum_studio.static_window import StaticWindow
        from vtkmodules.vtkCommonDataModel import VTK_QUADRATIC_TETRA

        study = tetra_specimen(curved=True, bending=True)
        root = self.root / "curved-bending"
        report = run_static(study, root)
        window = StaticWindow()
        self.addCleanup(window.close)
        window.show()
        window.set_study(study, "Curved quadratic tetrahedra")
        self.assertIn("C3D10", window.mesh_info.text())
        self.assertEqual(window.viewport.mesh.GetCellType(0), VTK_QUADRATIC_TETRA)
        self.assertEqual(window.viewport.mesh.GetCell(0).GetNumberOfPoints(), 10)
        window._present_result((study, report, root))
        window.channels.setCurrentIndex(1)
        self.assertEqual(window.table_model.rowCount(), 4 * len(study.elements))
        np.testing.assert_array_equal(
            window.table_model.rows[:, 1], np.tile(np.arange(1, 5), len(study.elements))
        )
        np.testing.assert_array_equal(window.table_model.rows[:, 2:8], report["stress"])
        window.scale.setText("10000")
        window._display(fit=True)
        np.testing.assert_allclose(
            window.viewport.display_positions,
            np.array(study.nodes) + 10000 * np.array(report["displacements"]),
            rtol=1e-14,
        )
        path = self.root / "tetra.ccx.json"
        with patch.object(
            QFileDialog, "getSaveFileName", return_value=(str(path), "JSON")
        ):
            window.save_study()
        self.assertEqual(load_study(path), study)
        window.mode.setCurrentIndex(1)
        self.assertEqual(window.result[0], study)

    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.catch = patch("sys.excepthook")
        self.exceptions = self.catch.start()
        self.addCleanup(self.catch.stop)
        self.addCleanup(self.check_exceptions)
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def check_exceptions(self):
        QTest.qWait(10)
        self.exceptions.assert_not_called()

    def wait_until(self, condition, seconds=10):
        deadline = time.monotonic() + seconds
        while not condition() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(condition(), "Qt process did not reach the expected state")

    def executable(self, name, body):
        path = self.root / name
        path.write_text(
            f'#!{sys.executable}\nimport sys, os, time\nif "-v" in sys.argv:\n print("Version 2.21", flush=True)\n sys.exit(201)\n{body}\n'
        )
        path.chmod(0o755)
        return path

    def test_direct_process_cancel_timeout_failure_and_old_result(self):
        controller = StaticController()
        self.addCleanup(controller.shutdown)
        sentinel = object()
        controller.last_result = sentinel
        problems, states, completed = [], [], []
        controller.problem.connect(problems.append)
        controller.busy_changed.connect(states.append)
        controller.completed.connect(completed.append)
        hanging = self.executable(
            "hanging",
            'from pathlib import Path\nPath("pid.txt").write_text(str(os.getpid()))\ntime.sleep(30)',
        )
        directory = self.root / "cancelled"
        controller.start(tension_example(), directory, executable=hanging)
        self.wait_until(lambda: (directory / "pid.txt").exists())
        pid = int((directory / "pid.txt").read_text())
        with self.assertRaises(RuntimeError):
            controller.start(
                tension_example(), self.root / "overlap", executable=hanging
            )
        controller.cancel()
        self.wait_until(lambda: controller.process is None)
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)
        self.assertIn("cancelled", problems[-1])
        self.assertIs(controller.last_result, sentinel)
        self.assertFalse((directory / "result.json").exists())
        controller.start(
            tension_example(), self.root / "timeout", executable=hanging, timeout=0.05
        )
        self.wait_until(lambda: controller.process is None)
        self.assertIn("timed out", problems[-1])
        failing = self.executable("failing", 'print("*ERROR: deliberate failure")')
        controller.start(tension_example(), self.root / "failure", executable=failing)
        self.wait_until(lambda: controller.process is None)
        self.assertIn("did not complete", problems[-1])
        self.assertTrue((self.root / "failure" / "solver.log").is_file())
        self.assertEqual(states, [True, False] * 3)
        self.assertFalse(completed)
        self.assertIs(controller.last_result, sentinel)

    def test_shutdown_stops_the_actual_solver(self):
        controller = StaticController()
        self.addCleanup(controller.shutdown)
        hanging = self.executable(
            "slow",
            'from pathlib import Path\nPath("pid.txt").write_text(str(os.getpid()))\ntime.sleep(30)',
        )
        output = self.root / "close"
        controller.start(tension_example(), output, executable=hanging)
        self.wait_until(lambda: (output / "pid.txt").exists())
        pid = int((output / "pid.txt").read_text())
        controller.shutdown()
        self.assertIsNone(controller.process)
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)
        self.assertFalse((output / "result.json").exists())

    @unittest.skipUnless(
        shutil.which("ccx"), "Install the separate CalculiX executable"
    )
    def test_open_archive_keeps_edited_study_and_rejects_bad_archive_transactionally(
        self,
    ):
        from vinkulum_studio.static_window import StaticWindow

        source = self.root / "previous calculation"
        captured = tension_example()
        expected = run_static(captured, source)
        window = StaticWindow()
        window._discard_allowed = lambda: True
        self.addCleanup(window.close)
        window.show()
        QTest.qWait(50)
        current = replace(
            captured, nodes=tuple((x + 4, y, z) for x, y, z in captured.nodes)
        )
        window.set_study(current, "Current edited study")
        window.young.setText("1e9")
        pending = window.edited_study()
        with (
            patch.object(
                QFileDialog,
                "getOpenFileName",
                return_value=(str(source / "result.json"), ""),
            ),
            patch(
                "vinkulum_studio.static_controller.identify_engine",
                side_effect=AssertionError(
                    "Archive viewing must not identify an executable"
                ),
            ),
        ):
            window.open_result_action.trigger()
            self.wait_until(lambda: not window.archive.busy)
        self.assertEqual(window.edited_study(), pending)
        self.assertEqual(window.result[0], captured)
        self.assertEqual(window.result[1], expected)
        self.assertIs(window.controller.last_result, window.result)
        self.assertEqual(window.mode.currentIndex(), 1)
        self.assertIn("Archived result loaded", window.status.text())
        np.testing.assert_allclose(
            window.viewport.display_positions,
            np.array(captured.nodes) + np.array(expected["displacements"]),
            rtol=1e-14,
            atol=1e-16,
        )
        previous = window.result
        position = window.viewport.display_positions.copy()
        damaged = self.root / "damaged"
        shutil.copytree(source, damaged)
        (damaged / "study.dat").write_text("invalid")
        with patch.object(
            QFileDialog,
            "getOpenFileName",
            return_value=(str(damaged / "result.json"), ""),
        ):
            window.open_result_action.trigger()
            self.wait_until(lambda: not window.archive.busy)
        self.assertIs(window.result, previous)
        self.assertEqual(window.edited_study(), pending)
        np.testing.assert_array_equal(window.viewport.display_positions, position)
        self.assertIn("Archive rejected", window.status.text())
        window.mode.setCurrentIndex(0)
        np.testing.assert_array_equal(window.viewport.display_positions, current.nodes)

    @unittest.skipUnless(
        shutil.which("ccx"), "Install the separate CalculiX executable"
    )
    def test_studio_entry_real_solve_scaled_view_capture_and_export(self):
        from vinkulum_studio.editor import EditorWindow

        editor = EditorWindow()
        editor._discard_allowed = lambda: True
        self.addCleanup(editor.close)
        editor.show()
        original = editor.project
        editor.commands["static_study"].trigger()
        window = editor._static_window
        window._discard_allowed = lambda: True
        QTest.qWait(100)
        self.addCleanup(
            lambda: window.close() if editor._static_window is window else None
        )
        window.output.setText(str(self.root))
        window.load_factor.setText("2")
        self.assertEqual(sum(r[1] for r in window.edited_study().forces), 2000)
        QTest.mouseClick(window.run_button, Qt.MouseButton.LeftButton)
        self.assertFalse(window.run_button.isEnabled())
        self.wait_until(lambda: window.result is not None)
        self.assertIsNone(window.controller.process)
        self.assertTrue(window.run_button.isEnabled())
        self.assertEqual(editor.project, original)
        study, report, directory = window.result
        self.assertEqual(load_study(directory / "study.json"), study)
        np.testing.assert_allclose(
            report["strain_energy_J"], 2000**2 / (2 * 210e9 * 0.01), rtol=5e-7
        )
        np.testing.assert_allclose(np.array(report["stress"])[:, 0], 200000, rtol=5e-7)
        self.assertEqual(window.table.model().rowCount(), 8)
        self.assertFalse(
            window.table.model().flags(window.table.model().index(0, 1))
            & Qt.ItemFlag.ItemIsEditable
        )
        window.scale.setText("10000")
        window._display(fit=True)
        np.testing.assert_allclose(
            window.viewport.display_positions,
            np.array(study.nodes) + 10000 * np.array(report["displacements"]),
            rtol=1e-14,
            atol=1e-16,
        )
        np.testing.assert_allclose(
            window.table_model.rows[:, 1:4], report["displacements"], rtol=0, atol=0
        )
        window.channels.setCurrentIndex(1)
        self.assertEqual(window.table_model.rowCount(), 8)
        np.testing.assert_array_equal(window.table_model.rows[:, 1], np.arange(1, 9))
        csv_path = self.root / "stress.csv"
        with patch.object(
            QFileDialog, "getSaveFileName", return_value=(str(csv_path), "CSV")
        ):
            window.export_table()
        with csv_path.open() as stream:
            rows = list(csv.reader(stream))
        self.assertIn("sxx [Pa]", rows[0])
        self.assertEqual(float(rows[1][2]), report["stress"][0][0])
        window.young.setText("1e9")
        self.assertIs(window.result[1], report)
        self.assertEqual(window.result[0].young_pa, 210e9)
        self.assertEqual(window.edited_study().young_pa, 1e9)
        moved = replace(
            tension_example(),
            nodes=tuple((x + 3, y, z) for x, y, z in tension_example().nodes),
        )
        window.set_study(moved, "A different loaded mesh")
        np.testing.assert_allclose(
            window.viewport.display_positions, moved.nodes, rtol=0, atol=0
        )
        window.mode.setCurrentIndex(1)
        np.testing.assert_allclose(
            window.viewport.display_positions,
            np.array(study.nodes) + 10000 * np.array(report["displacements"]),
            rtol=1e-14,
            atol=1e-16,
        )
        self.assertIs(window.result[1], report)
        self.assertIn("captured", window.result_summary.text())
        self.assertEqual(window.edited_study().young_pa, moved.young_pa)
        failed = self.executable("gui-failure", 'print("*ERROR: deliberate failure")')
        window.executable.setText(str(failed))
        window.start()
        self.wait_until(lambda: window.controller.process is None)
        self.assertIs(window.result[1], report)
        self.assertTrue(window.run_button.isEnabled())
        self.assertIn("did not complete", window.status.text())
        self.assertEqual(window.result[0].young_pa, 210e9)
        QTest.qWait(100)
        window.viewport.render()
        # QScreen captures the native OpenGL child; QWidget.grab() can instead
        # read undefined backing-store bytes beneath that child on X11.
        image = self.app.primaryScreen().grabWindow(window.winId()).toImage()
        self.assertFalse(image.isNull())
        field_path = self.root / "field.png"
        window.viewport.screenshot(field_path)
        field = QImage(str(field_path))
        background = field.pixelColor(10, 10)
        centre = field.pixelColor(field.width() // 2, field.height() // 2)
        self.assertLess(
            max(background.red(), background.green(), background.blue()), 100
        )
        self.assertGreater(max(centre.red(), centre.green(), centre.blue()), 120)
        path = os.environ.get("VINKULUM_STATIC_SCREENSHOT")
        if path:
            image.save(path)
            window.viewport.screenshot(str(Path(path).with_stem("static-field")))
        window.close()
        self.assertIsNone(editor._static_window)
