"""Real CAD preview transactions, editing precision and distinct worker failures."""

import importlib.util
import os
import sys
import tempfile
import time
import unittest
from dataclasses import asdict
from pathlib import Path
from unittest.mock import patch

import numpy as np
from PySide6.QtCore import QSettings, Qt, QTimer
from PySide6.QtGui import QPalette
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QDialog
from vinkulum_studio.cad_controller import CadController
from vinkulum_studio.cad_history_dialog import CadHistoryDialog
from vinkulum_studio.document import Body, Project, joint_at, new_id
from vinkulum_studio.editor import EditorWindow


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires a real desktop context"
)
class CadFeatureWorkspace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.errors = patch("sys.excepthook")
        self.caught = self.errors.start()
        self.addCleanup(self.errors.stop)
        self.addCleanup(self.check_exceptions)

    def check_exceptions(self):
        QTest.qWait(10)
        self.caught.assert_not_called()

    def wait(self, condition, seconds=30):
        deadline = time.monotonic() + seconds
        while not condition() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(condition(), "CAD transaction did not reach the expected state")

    def edit(self, field, value):
        field.setFocus()
        field.selectAll()
        QTest.keyClicks(field, str(value))

    def controller(self):
        result = CadController()
        self.addCleanup(result.shutdown)
        return result

    def script(self, name, body):
        path = self.root / name
        path.write_text(f"#!{sys.executable}\n" + body)
        path.chmod(0o755)
        return path

    def test_failure_causes_and_cancellation_preserve_previous_result(self):
        c = self.controller()
        sentinel = object()
        c.last_result = sentinel
        failures = []
        c.failed.connect(lambda kind, message: failures.append((kind, message)))
        slow = self.script("slow", "import time\ntime.sleep(30)\n")
        c.start({"operation": "box"}, executable=slow, arguments=[])
        self.wait(lambda: c.process is not None and c.process.processId() > 0)
        pid = c.process.processId()
        c.cancel()
        self.wait(lambda: c.process is None)
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)
        self.assertEqual(failures[-1][0], "cancelled")
        self.assertIsNone(c.temporary)
        c.start({}, executable=slow, arguments=[], timeout=0.05)
        self.wait(lambda: c.process is None)
        self.assertEqual(failures[-1][0], "timed_out")
        crashed = self.script(
            "crash", "import os,signal\nos.kill(os.getpid(), signal.SIGTERM)\n"
        )
        c.start({}, executable=crashed, arguments=[])
        self.wait(lambda: c.process is None)
        self.assertEqual(failures[-1][0], "crashed")
        missing = self.root / "does-not-exist"
        c.start({}, executable=missing, arguments=[])
        self.wait(lambda: c.process is None)
        self.assertEqual(failures[-1][0], "start_failed")
        for payload, expected in (
            ({"status": "failed", "error": "Impossible fillet"}, "operation_failed"),
            (
                {"status": "completed", "request_sha256": "0" * 64, "body": {}},
                "invalid_result",
            ),
        ):
            fake = self.script(
                expected,
                f"import sys,json\nfrom pathlib import Path\nPath(sys.argv[-1]).write_text(json.dumps({payload!r}))\n",
            )
            c.start({}, executable=fake, arguments=[])
            self.wait(lambda: c.process is None)
            self.assertEqual(failures[-1][0], expected)
        self.assertIs(c.last_result, sentinel)
        c.start({}, executable=slow, arguments=[])
        self.wait(lambda: c.process is not None and c.process.processId() > 0)
        pid = c.process.processId()
        c.shutdown()
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)

    def test_invalid_inputs_keep_selection_and_disable_preview(self):
        body = Body(new_id(), "Stock", dimensions=(0.12, 0.075, 0.02), mass=1.404)
        project = Project(new_id(), "Unchanged", bodies=(body,))
        dialog = CadHistoryDialog(project, body.id)
        self.addCleanup(dialog.reject)
        dialog.show()
        QTest.qWait(50)
        self.edit(dialog.dimensions[0], "invalid")
        self.assertFalse(dialog.preview_button.isEnabled())
        self.assertFalse(dialog.apply_button.isEnabled())
        self.edit(dialog.dimensions[0], "150.000000000001")
        self.assertEqual(
            dialog.edited_recipe().features[0].dimensions_mm[0], 150.000000000001
        )
        self.assertTrue(dialog.preview_button.isEnabled())
        self.edit(dialog.density, "nan")
        self.assertFalse(dialog.preview_button.isEnabled())
        dialog.accept()
        self.assertTrue(dialog.isVisible())
        self.assertIsNone(dialog.result_body)
        self.assertEqual(dialog.project, project)
        # The dialog has its own dark palette; opening a window must not leave
        # white form backgrounds behind white labels.
        palette = dialog.form.parentWidget().palette()
        self.assertGreater(
            palette.color(QPalette.ColorRole.WindowText).lightness(),
            palette.color(QPalette.ColorRole.Window).lightness(),
        )

    @unittest.skipUnless(importlib.util.find_spec("build123d"), "Install OCCT 8")
    def test_input_hash_does_not_admit_wrong_mass_frame_or_recipe(self):
        c = self.controller()
        c.start({"operation": "box", "dimensions_mm": [120, 75, 20], "density": 7800})
        self.wait(lambda: c.process is None)
        self.assertIsNone(c.failure_kind, c.last_log)
        previous = c.last_result
        source = Body.from_dict(previous["body"])
        failures = []
        c.failed.connect(lambda kind, message: failures.append(message))
        for modification, diagnostic in (
            ("body['mass'] *= 2", "mass properties"),
            ("body['position'][0] += 0.01", "world origin"),
            (
                "body['cad']['recipe']['features'][0]['name'] = 'Wrong graph'",
                "feature graph",
            ),
            ("body['cad']['volume_m3'] = True", "volume"),
        ):
            fake = self.script(
                "wrong-result",
                f"import sys,json,hashlib\nfrom pathlib import Path\np=Path(sys.argv[-2])\nrequest=json.loads(p.read_text())\nbody=request['a']\n{modification}\nPath(sys.argv[-1]).write_text(json.dumps({{'status':'completed','body':body,'request_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}}))\n",
            )
            c.start(
                {
                    "operation": "regenerate",
                    "a": asdict(source),
                    "recipe": asdict(source.cad.recipe),
                    "density": 7800,
                },
                executable=fake,
                arguments=[],
            )
            self.wait(lambda: c.process is None)
            self.assertEqual(c.failure_kind, "invalid_result")
            self.assertIn(diagnostic, failures[-1])
            self.assertIs(c.last_result, previous)

    @unittest.skipUnless(
        importlib.util.find_spec("build123d")
        and os.environ.get("VINKULUM_PINOCCHIO_PYTHON"),
        "Requires OCCT 8 and the separate current Pinocchio worker",
    )
    def test_parametric_cad_to_pinocchio_retains_mass_tensor_and_design_history(self):
        import math

        from vinkulum_studio.pinocchio_controller import PinocchioController

        c = self.controller()
        c.start(
            {
                "operation": "box",
                "name": "CAD rod",
                "dimensions_mm": [20, 30, 1000],
                "position_mm": [0, 0, -500],
                "density": 7800,
            }
        )
        self.wait(lambda: c.process is None)
        self.assertIsNone(c.failure_kind, c.last_log)
        body = Body.from_dict(c.last_result["body"])
        project = Project(new_id(), "CAD rod operators", bodies=(body,), gravity=(0.0, 0.0, -9.81))
        project = project.replace_object(
            joint_at(project, "pivot", None, body.id, axis=(0, 1, 0))
        )
        controller = PinocchioController()
        self.addCleanup(controller.shutdown)
        failures = []
        controller.problem.connect(failures.append)
        angle = 0.2
        controller.start(
            project,
            {"q": [angle]},
            self.root / "operators",
            interpreter=os.environ["VINKULUM_PINOCCHIO_PYTHON"],
        )
        self.wait(lambda: controller.process is None)
        self.assertFalse(failures, failures)
        result = controller.last_result
        self.assertIsNotNone(result)
        self.assertEqual(result.project, project)
        self.assertEqual(result.project.bodies[0].cad.recipe, body.cad.recipe)
        mass = 7800 * 0.02 * 0.03 * 1.0
        inertia_about_com_y = mass * (0.02**2 + 1.0**2) / 12
        self.assertAlmostEqual(
            result.report["mass_matrix"][0][0],
            inertia_about_com_y + mass * 0.5**2,
            places=10,
        )
        self.assertAlmostEqual(
            result.report["intrinsic_bias"][0],
            mass * 9.81 * 0.5 * math.sin(angle),
            places=10,
        )
        np.testing.assert_allclose(
            result.report["bodies"][0]["position_m"],
            [-0.5 * math.sin(angle), 0, -0.5 * math.cos(angle)],
            atol=1e-12,
        )

    @unittest.skipUnless(importlib.util.find_spec("build123d"), "Install OCCT 8")
    def test_real_preview_is_not_a_commit_and_keeps_exact_edits(self):
        c = self.controller()
        c.start(
            {
                "operation": "box",
                "name": "Stock",
                "dimensions_mm": [120, 75, 20],
                "density": 7800,
            }
        )
        self.wait(lambda: c.process is None)
        self.assertIsNone(c.failure_kind, c.last_log)
        body = Body.from_dict(c.last_result["body"])
        # Add a dependent feature so invalid input cannot silently disappear
        # when selecting another row.
        c.start(
            {"operation": "fillet", "a": asdict(body), "radius_mm": 2, "density": 7800}
        )
        self.wait(lambda: c.process is None)
        body = Body.from_dict(c.last_result["body"])
        project = Project(new_id(), "Model", bodies=(body,))
        dialog = CadHistoryDialog(project, body.id)
        self.addCleanup(dialog.reject)
        dialog.show()
        dialog.activateWindow()
        QTest.qWait(80)
        self.edit(dialog.dimensions[0], "invalid")
        dialog.tree.setCurrentItem(dialog.tree.topLevelItem(1))
        self.assertIs(dialog.tree.currentItem(), dialog.tree.topLevelItem(0))
        self.assertEqual(dialog.dimensions[0].text(), "invalid")
        self.edit(dialog.dimensions[0], "150.000000000001")
        self.edit(dialog.rotation.components[2], 30)
        QTest.keyClick(dialog, Qt.Key.Key_Return, Qt.KeyboardModifier.ControlModifier)
        self.wait(
            lambda: (
                dialog.controller.process is None and dialog.preview_body is not None
            )
        )
        self.assertEqual(dialog.project, project)
        self.assertIsNone(dialog.result_body)
        self.assertTrue(dialog.apply_button.isEnabled())
        self.assertEqual(
            dialog.preview_body.cad.recipe.features[0].dimensions_mm[0],
            150.000000000001,
        )
        self.assertGreater(dialog.preview_body.mass, body.mass)
        # Reverting a rotation after a preview must restore the orientation that
        # was originally presented in the fields, not the latest preview matrix.
        self.edit(dialog.rotation.components[2], 0)
        np.testing.assert_array_equal(
            dialog.edited_recipe().features[0].orientation,
            body.cad.recipe.features[0].orientation,
        )
        self.assertFalse(dialog.apply_button.isEnabled())
        previous = dialog.preview_body
        dialog.tree.setCurrentItem(dialog.tree.topLevelItem(1))
        self.edit(dialog.dimensions[0], 1000)
        QTest.mouseClick(dialog.preview_button, Qt.MouseButton.LeftButton)
        self.wait(
            lambda: (
                dialog.controller.process is None
                and dialog.controller.failure_kind is not None
            )
        )
        self.assertEqual(dialog.controller.failure_kind, "operation_failed")
        self.assertIs(dialog.preview_body, previous)
        self.assertFalse(dialog.apply_button.isEnabled())
        self.assertEqual(dialog.project, project)
        self.edit(dialog.dimensions[0], 2)
        QTest.mouseClick(dialog.preview_button, Qt.MouseButton.LeftButton)
        self.wait(
            lambda: (
                dialog.controller.process is None and dialog.apply_button.isEnabled()
            )
        )
        QTest.mouseClick(dialog.apply_button, Qt.MouseButton.LeftButton)
        self.assertEqual(dialog.result(), QDialog.DialogCode.Accepted)
        self.assertEqual(dialog.result_body.id, body.id)

    @unittest.skipUnless(importlib.util.find_spec("build123d"), "Install OCCT 8")
    def test_editor_shortcut_commits_one_undoable_transaction(self):
        window = EditorWindow(
            QSettings(str(self.root / "settings.ini"), QSettings.Format.IniFormat)
        )
        window._discard_allowed = lambda: True
        self.addCleanup(window.close)
        window.show()
        body = Body(new_id(), "Stock", dimensions=(0.12, 0.075, 0.02), mass=1.404)
        project = Project(new_id(), "Edit stock", bodies=(body,))
        window._commit(project)
        window.select_object(body.id)
        before = window.project
        driven = []

        def interact():
            dialog = QApplication.activeModalWidget()
            try:
                self.assertIsInstance(dialog, CadHistoryDialog)
                self.edit(dialog.dimensions[0], 150)
                dialog.start_preview()
                self.wait(
                    lambda: (
                        dialog.controller.process is None
                        and dialog.preview_body is not None
                    )
                )
                self.assertEqual(window.project, before)
                dialog.accept()
                driven.append(True)
            finally:
                if dialog is not None and dialog.isVisible():
                    dialog.reject()

        window.activateWindow()
        QTest.qWait(80)
        QTimer.singleShot(50, interact)
        QTest.keyClick(
            window,
            Qt.Key.Key_H,
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier,
        )
        self.assertEqual(driven, [True])
        after = window.project
        self.assertEqual(after.revision, before.revision + 1)
        self.assertIsNotNone(after.bodies[0].cad.recipe)
        self.assertAlmostEqual(after.bodies[0].dimensions[0], 0.15, places=10)
        window.commands["undo"].trigger()
        self.assertEqual(window.project.bodies, before.bodies)
        window.commands["redo"].trigger()
        self.assertEqual(window.project.bodies, after.bodies)
