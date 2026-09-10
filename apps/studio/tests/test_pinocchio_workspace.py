"""Optional-worker lifecycle and captured operator validation in the GUI process."""

import csv
import json
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
from PySide6.QtCore import QSettings, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox
from shiboken6 import isValid
from vinkulum_studio.articulated_result import load_operators
from vinkulum_studio.articulated_window import ArticulatedWindow
from vinkulum_studio.document import Law, Load, new_id, pendulum
from vinkulum_studio.examples3d import double_pendulum
from vinkulum_studio.pinocchio_controller import PinocchioController


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1",
    "Requires the Studio desktop test environment",
)
class PinocchioWorkspace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.exceptions = patch("sys.excepthook")
        self.catch = self.exceptions.start()
        self.addCleanup(self.exceptions.stop)
        self.addCleanup(self.assert_no_exception)
        self.project = double_pendulum()

    def assert_no_exception(self):
        QTest.qWait(10)
        self.catch.assert_not_called()

    def wait_until(self, condition, seconds=10):
        deadline = time.monotonic() + seconds
        while not condition() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(
            condition(), "Pinocchio controller did not reach the expected state"
        )

    def script(self, name, source):
        path = self.root / name
        path.write_text(f"#!{sys.executable}\n" + source)
        path.chmod(0o755)
        return path

    def controller(self):
        c = PinocchioController()
        self.addCleanup(c.shutdown)
        return c

    def window(self):
        window = ArticulatedWindow(
            self.project,
            settings=QSettings(
                str(self.root / "settings.ini"), QSettings.Format.IniFormat
            ),
        )
        self.addCleanup(self.close_window, window)
        window.show()
        window.activateWindow()
        QTest.qWait(80)
        return window

    def close_window(self, window):
        if isValid(window):
            with patch.object(
                QMessageBox, "question", return_value=QMessageBox.StandardButton.Discard
            ):
                window.close()
            QTest.qWait(10)

    def edit(self, field, value):
        field.setFocus()
        field.selectAll()
        QTest.keyClicks(field, str(value))

    def test_unavailable_worker_invalid_edits_and_model_admission(self):
        window = self.window()
        window.interpreter.clear()
        self.assertFalse(window.run_action.isEnabled())
        self.assertIn("Execution settings", window.pending.text())
        window.interpreter.setText(sys.executable)
        self.assertTrue(window.run_action.isEnabled())
        self.edit(window.state_fields["q"][0], "invalid")
        self.assertFalse(window.run_action.isEnabled())
        with patch.object(
            QMessageBox, "question", return_value=QMessageBox.StandardButton.Cancel
        ) as question:
            window.open_example()
        question.assert_called_once()

        self.assertEqual(window.state_fields["q"][0].text(), "invalid")
        with patch.object(
            QMessageBox, "question", return_value=QMessageBox.StandardButton.Discard
        ):
            window.open_example()
        self.assertTrue(window.run_action.isEnabled())
        window.set_project(pendulum())  # G0's spherical joint is outside this adapter.
        self.assertFalse(window.run_action.isEnabled())
        self.assertEqual(window.input_table.rowCount(), 0)
        self.assertIn("revolute and prismatic", window.admission.text())
        self.edit(window.load_time, "invalid")
        with patch.object(
            QMessageBox, "question", return_value=QMessageBox.StandardButton.Cancel
        ) as question:
            self.assertFalse(window.close())
        question.assert_called_once()

    @unittest.skipUnless(
        os.environ.get("VINKULUM_PINOCCHIO_PYTHON"),
        "Set the separate qualified worker Python",
    )
    def test_loaded_derivative_tables_and_legacy_missing_channels(self):
        window = self.window()
        window.interpreter.setText(os.environ["VINKULUM_PINOCCHIO_PYTHON"])
        window.output.setText(str(self.root / "new-derivatives"))
        window.start()
        self.wait_until(lambda: window.controller.process is None, seconds=30)
        result = window.controller.last_result
        self.assertIsNotNone(result)
        for index, channel in (
            (6, "external_effort_derivatives"),
            (7, "loaded_inverse_derivatives"),
        ):
            window.channel.setCurrentIndex(index)
            np.testing.assert_array_equal(
                window.table.model().values, result.report[channel]["q"]
            )
            self.assertTrue(window.export_button.isEnabled())
            self.assertIn("row effort / column coordinate", window.table_note.text())
            self.assertIn("[rad]", window.table.model().headers[1])
            path = self.root / f"{channel}.csv"
            with patch.object(
                QFileDialog, "getSaveFileName", return_value=(str(path), "CSV")
            ):
                window.export_table()
            with path.open() as stream:
                rows = list(csv.reader(stream))
            np.testing.assert_array_equal(
                [[float(value) for value in row[1:]] for row in rows[1:]],
                result.report[channel]["q"],
            )
        old = (
            Path(__file__).resolve().parents[3]
            / "examples/studio/articulated/double-pendulum"
        )
        window._completed(load_operators(old))
        for index in (6, 7):
            window.channel.setCurrentIndex(index)
            self.assertEqual(window.table.model().rowCount(), 0)
            self.assertFalse(window.export_button.isEnabled())
            self.assertIn("does not contain", window.table_note.text())
        window.channel.setCurrentIndex(2)
        self.assertEqual(window.table.model().rowCount(), 2)
        self.assertTrue(window.export_button.isEnabled())

    def test_editor_shortcut_close_veto_and_real_cancel_button(self):
        from vinkulum_studio.editor import EditorWindow

        editor = EditorWindow(
            QSettings(str(self.root / "editor.ini"), QSettings.Format.IniFormat)
        )
        self.addCleanup(self.close_window, editor)
        editor.show()
        editor.load_example("Double pendulum")
        editor.activateWindow()
        QTest.qWait(80)
        QTest.keyClick(
            editor,
            Qt.Key.Key_P,
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier,
        )
        QTest.qWait(50)
        window = editor._articulated_window
        self.assertIsNotNone(window)
        self.edit(window.state_fields["q"][0], 0.2)
        editor.commands["articulated"].trigger()
        self.assertIs(editor._articulated_window, window)
        self.assertEqual(window.edited_state()["q"][0], 0.2)
        with (
            patch.object(editor, "_discard_allowed", return_value=True),
            patch.object(
                QMessageBox, "question", return_value=QMessageBox.StandardButton.Cancel
            ),
        ):
            self.assertFalse(editor.close())
        self.assertIs(editor._articulated_window, window)
        slow = self.script("slow-ui", "import time\ntime.sleep(30)\n")
        window.interpreter.setText(str(slow))
        window.output.setText(str(self.root / "runs"))
        QTest.mouseClick(window.run_button, Qt.MouseButton.LeftButton)
        self.wait_until(
            lambda: (
                window.controller.process is not None
                and window.controller.process.processId() != 0
            )
        )
        pid = window.controller.process.processId()
        self.assertFalse(window.input_table.isEnabled())
        self.assertFalse(window.open_result_action.isEnabled())
        self.assertTrue(window.cancel_button.isEnabled())
        QTest.mouseClick(window.cancel_button, Qt.MouseButton.LeftButton)
        self.wait_until(lambda: window.controller.process is None)
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)
        self.assertIn("cancelled", window.status.text())
        self.assertTrue(window.input_table.isEnabled())
        self.close_window(window)
        self.assertIsNone(editor._articulated_window)

    @unittest.skipUnless(
        os.environ.get("VINKULUM_PINOCCHIO_PYTHON"),
        "Set the separate qualified worker Python",
    )
    def test_workspace_evaluate_export_and_reopen_preserves_current_inputs(self):
        self.project = self.project.replace_object(
            Load(
                new_id(),
                "Ramp force",
                self.project.bodies[-1].id,
                point=(0.02, 0.0, 0.1),
                force=(Law("lineaire", (1.0, 4.0)), Law(), Law()),
            )
        )
        window = self.window()
        window.output.setText(str(self.root / "runs"))
        state = {
            "q": [0.2000000000001234, -0.3],
            "velocity": [0.1, 0.4],
            "acceleration": [0.7, 0.2],
            "effort": [0.9, 0.1],
        }
        for key, values in state.items():
            for field, value in zip(window.state_fields[key], values):
                self.edit(field, value)
        self.edit(window.load_time, 0.5)
        captured = window.edited_state()
        input_path = self.root / "input.json"
        with patch.object(
            QFileDialog, "getSaveFileName", return_value=(str(input_path), "")
        ):
            window.save_input()
        self.assertEqual(json.loads(input_path.read_text())["state"], captured)
        window.activateWindow()
        QTest.qWait(50)
        QTest.keyClick(window, Qt.Key.Key_Return, Qt.KeyboardModifier.ControlModifier)
        self.wait_until(
            lambda: window.controller.process is None and window.result is not None,
            seconds=30,
        )
        result = window.result
        self.assertEqual(result.state, captured)
        self.assertEqual(window.mode.currentText(), "Captured state")
        self.assertEqual(result.display_project.loads[0].force[0].value(0), 3.0)
        self.assertEqual(result.project.loads[0].force[0].kind, "lineaire")
        self.assertNotIn("pinocchio", sys.modules)
        for channel in range(6):
            window.channel.setCurrentIndex(channel)
            model = window.table.model()
            self.assertGreater(model.rowCount(), 0)
            self.assertFalse(model.values.flags.writeable)
            self.assertEqual(
                model.data(model.index(0, 1), Qt.ItemDataRole.ToolTipRole),
                repr(float(model.values[0, 0])),
            )
        window.body_choice.setCurrentIndex(1)
        np.testing.assert_array_equal(
            window.table.model().values, result.report["bodies"][1]["jacobian"]
        )
        csv_path = self.root / "jacobian.csv"
        with patch.object(
            QFileDialog, "getSaveFileName", return_value=(str(csv_path), "")
        ):
            window.export_table()
        with csv_path.open(newline="") as stream:
            rows = list(csv.reader(stream))
        self.assertEqual(rows[0], list(window.table.model().headers))
        np.testing.assert_array_equal(
            np.array([row[1:] for row in rows[1:]], dtype=float),
            window.table.model().values,
        )
        # A newly captured model and edited state never overwrite an old result.
        window.set_project(replace(self.project, name="Next mechanism"))
        self.assertIs(window.result, result)
        self.assertEqual(window.mode.currentText(), "Authored model")
        self.edit(window.state_fields["q"][0], 0.8)
        edited = window.edited_state()
        window.interpreter.clear()
        with patch.object(
            QFileDialog,
            "getOpenFileName",
            return_value=(str(result.directory / "result.json"), ""),
        ):
            window.open_result()
            self.wait_until(lambda: not window.archive.busy)
        self.assertEqual(window.project.name, "Next mechanism")
        self.assertEqual(window.edited_state(), edited)
        self.assertEqual(window.result.state, captured)
        self.assertIn("Inputs differ", window.pending.text())
        self.assertFalse(window.run_button.isEnabled())
        previous = window.result
        invalid = self.root / "invalid"
        shutil.copytree(result.directory, invalid)
        data = json.loads((invalid / "result.json").read_text())
        data["mass_matrix"][0][0] = -1.0
        (invalid / "result.json").write_text(json.dumps(data))
        with patch.object(
            QFileDialog,
            "getOpenFileName",
            return_value=(str(invalid / "result.json"), ""),
        ):
            window.open_result()
            self.wait_until(lambda: not window.archive.busy)
        self.assertIs(window.result, previous)
        self.assertIn("positive definite", window.status.text())
        with (
            patch.object(
                QFileDialog, "getOpenFileName", return_value=(str(input_path), "")
            ),
            patch.object(
                QMessageBox, "question", return_value=QMessageBox.StandardButton.Cancel
            ),
        ):
            window.open_input()
            self.wait_until(lambda: not window.archive.busy)
        self.assertEqual(window.edited_state(), edited)
        with (
            patch.object(
                QFileDialog, "getOpenFileName", return_value=(str(input_path), "")
            ),
            patch.object(
                QMessageBox, "question", return_value=QMessageBox.StandardButton.Discard
            ),
        ):
            window.open_input()
            self.wait_until(lambda: not window.archive.busy)
        self.assertEqual(window.project, self.project)
        self.assertEqual(window.edited_state(), captured)

    def test_cancel_timeout_start_failure_and_environment_isolation(self):
        c = self.controller()
        previous = object()
        c.last_result = previous
        errors = []
        c.problem.connect(errors.append)
        slow = self.script(
            "slow",
            "import os,time,json\nfrom pathlib import Path\nPath('pid.txt').write_text(str(os.getpid()))\nPath('env.json').write_text(json.dumps({key:os.environ.get(key) for key in ('PYTHONPATH','PYTHONHOME','VIRTUAL_ENV','OPENBLAS_NUM_THREADS')}))\ntime.sleep(30)\n",
        )
        root = self.root / "cancel"
        with patch.dict(
            os.environ,
            {
                "PYTHONPATH": "/tmp/irrelevant-gui-path",
                "PYTHONHOME": "/tmp/irrelevant-python-home",
            },
        ):
            c.start(self.project, {}, root, interpreter=slow)
        self.wait_until(lambda: (root / "pid.txt").is_file())
        pid = int((root / "pid.txt").read_text())
        self.assertEqual(
            json.loads((root / "env.json").read_text()),
            {
                "PYTHONPATH": None,
                "PYTHONHOME": None,
                "VIRTUAL_ENV": None,
                "OPENBLAS_NUM_THREADS": "1",
            },
        )
        c.cancel()
        self.wait_until(lambda: c.process is None)
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)
        self.assertIn("cancelled", errors[-1])
        self.assertIs(c.last_result, previous)
        c.start(self.project, {}, self.root / "timeout", interpreter=slow, timeout=0.05)
        self.wait_until(lambda: c.process is None)
        self.assertIn("timed out", errors[-1])
        broken = self.script("broken", "pass\n")
        broken.write_text("#!/vinkulum-missing-interpreter\n")
        c.start(self.project, {}, self.root / "failure", interpreter=broken)
        self.wait_until(lambda: c.process is None)
        self.assertIn("Could not start", errors[-1])
        self.assertIs(c.last_result, previous)
        root = self.root / "shutdown"
        c.start(self.project, {}, root, interpreter=slow)
        self.wait_until(lambda: (root / "pid.txt").is_file())
        pid = int((root / "pid.txt").read_text())
        c.shutdown()
        self.assertIsNone(c.process)
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)

    @unittest.skipUnless(
        os.environ.get("VINKULUM_PINOCCHIO_PYTHON"),
        "Set the separate qualified worker Python",
    )
    def test_real_worker_capture_and_rejection_of_inconsistent_results(self):
        c = self.controller()
        errors = []
        c.problem.connect(errors.append)
        state = {
            "q": [0.2, -0.3],
            "velocity": [0.1, 0.4],
            "acceleration": [0.7, 0.2],
            "effort": [0.9, 0.1],
        }
        original = json.loads(json.dumps(state))
        c.start(
            self.project,
            state,
            self.root / "valid",
            interpreter=os.environ["VINKULUM_PINOCCHIO_PYTHON"],
        )
        state["q"][0] = 99
        self.wait_until(lambda: c.process is None, seconds=30)
        self.assertFalse(errors, errors)
        result = c.last_result
        self.assertIsNotNone(result)
        self.assertEqual(result.state["q"], original["q"])
        self.assertEqual(result.project, self.project)
        self.assertEqual(result.report["scientific_status"], "NotAssessed")
        baseline = json.loads((result.directory / "result.json").read_text())
        before = {
            p.name: (p.read_bytes(), p.stat().st_mtime_ns)
            for p in result.directory.iterdir()
            if p.is_file()
        }
        reopened = load_operators(result.directory)
        self.assertEqual(reopened.report, result.report)
        self.assertEqual(
            before,
            {
                p.name: (p.read_bytes(), p.stat().st_mtime_ns)
                for p in result.directory.iterdir()
                if p.is_file()
            },
        )
        edits = (
            lambda r: r.update(scientific_status="Certified"),
            lambda r: r["coordinate_order"][0].update(coordinate_unit="deg"),
            lambda r: r["mass_matrix"][0].__setitem__(0, True),
            lambda r: r["inverse_effort"].__setitem__(0, r["inverse_effort"][0] + 1),
            lambda r: r["bodies"][0]["position_m"].__setitem__(0, 20),
            lambda r: r["conventions"].__setitem__("jacobians", "at the world origin"),
            lambda r: r["state"]["q"].__setitem__(0, 0.8),
            lambda r: r["intrinsic_inverse_derivatives"]["acceleration"][0].__setitem__(
                0, 0
            ),
            lambda r: r["mass_matrix"][0].__setitem__(0, -1.0),
            lambda r: r["checks"].update(relative_tolerance=1.0),
            lambda r: r.update(environment={}),
            lambda r: r.update(loads=[{}]),
        )
        for i, edit in enumerate(edits):
            folder = self.root / f"bad-{i}"
            shutil.copytree(result.directory, folder)
            altered = json.loads(json.dumps(baseline))
            edit(altered)
            (folder / "result.json").write_text(json.dumps(altered))
            with self.subTest(i=i), self.assertRaises((ValueError, TypeError)):
                load_operators(folder)
        invalid = self.script(
            "invalid",
            f"import shutil,sys,json\nfrom pathlib import Path\nroot=Path(sys.argv[sys.argv.index('--output')+1])\nshutil.copytree({str(result.directory)!r},root)\np=root/'result.json'\ndata=json.loads(p.read_text())\ndata['inverse_effort'][0]+=1\np.write_text(json.dumps(data))\n",
        )
        c.start(self.project, original, self.root / "bad-worker", interpreter=invalid)
        self.wait_until(lambda: c.process is None)
        self.assertIs(c.last_result, result)
        self.assertIn("Inconsistent Pinocchio", errors[-1])
