"""G0 requirements: real solver recipe, independent reference, supervised faults."""

import csv
from dataclasses import FrozenInstanceError, asdict, replace
import io
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("RAYON_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

from PySide6.QtCore import Qt, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from scipy.integrate import solve_ivp

from vinkulum_studio.controller import Controller
from vinkulum_studio.model import (G, INERTIA, Parameters, Result,
                                    load_parameters, read_json,
                                    save_parameters, write_json)
from vinkulum_studio.worker import simulate
from vinkulum_studio.window import MainWindow


class Contracts(unittest.TestCase):
    def test_GUI02_GUI07_strict_inputs_and_version(self):
        for change in ({"mass": True}, {"length": float("nan")}, {"step": 0.},
                       {"duration": 600., "step": .000001}, {"length": "1"},
                       {"mass": 10**400}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                Parameters(**(asdict(Parameters()) | change))
        with self.assertRaises(FrozenInstanceError):
            Parameters().mass = 1
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "params.json"
            save_parameters(path, Parameters())
            self.assertEqual(load_parameters(path), Parameters())
            original = read_json(path)
            for data in (original | {"schema_version": 2}, original | {"schema_version": True},
                         original | {"extra": "code"},
                         original | {"parameters": asdict(Parameters()) | {"script": "exit()"}}):
                write_json(path, data)
                with self.assertRaises(ValueError):
                    load_parameters(path)
            for content in ('{"x": 1, "x": 2}', '{"x": NaN}', '[' * 2000,
                            ' ' * 16_385, '{"x": Infinity}'):
                path.write_text(content)
                with self.assertRaises(ValueError):
                    load_parameters(path)

    def test_T10_failed_replace_preserves_previous_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "params.json"
            save_parameters(path, Parameters())
            before = path.read_bytes()
            with patch("vinkulum_studio.model.os.replace", side_effect=OSError("disk failure")):
                with self.assertRaises(OSError):
                    save_parameters(path, Parameters(length=2.))
            self.assertEqual(path.read_bytes(), before)
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_T07_reference_and_temporal_refinement(self):
        # Thresholds fixed before execution: < 4e-4 rad at h=.005;
        # refinement ratio > 3.5 on each pair (second order expected).
        p = Parameters()
        coefficient = p.mass * G * p.length / (p.mass * p.length**2 + INERTIA)
        reference = solve_ivp(lambda t, y: (y[1], -coefficient * math.sin(y[0])),
                              (0, p.duration), (math.radians(p.angle_deg), 0),
                              method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
        self.assertTrue(reference.success)
        errors = []
        for h in (.02, .01, .005):
            parameters = replace(p, step=h)
            result = Result.from_dict(simulate(parameters, "reference"), "reference", parameters)
            error = max(abs(row[3] - reference.sol(row[0])[0]) for row in result.samples)
            errors.append(float(error))
        ratios = [errors[i] / errors[i + 1] for i in range(2)]
        print("\nT07 errors(rad):", errors, "ratios:", ratios, flush=True)
        self.assertLess(errors[-1], 4e-4)
        self.assertGreater(min(ratios), 3.5)

    def test_GUI05_reject_mismatched_corrupt_and_partial_result(self):
        p = Parameters(duration=.1)
        data = simulate(p, "expected")
        for mutate in (lambda d: d.update(run_id="wrong"),
                       lambda d: d["samples"].pop(),
                       lambda d: d["samples"][-1].__setitem__(1, float("nan")),
                       lambda d: d["samples"][-1].__setitem__(3, 2.0),
                       lambda d: d["samples"][1].__setitem__(0, 0.)):
            copy = json.loads(json.dumps(data))
            mutate(copy)
            with self.assertRaises(ValueError):
                Result.from_dict(copy, "expected", p)

    def test_GUI07_export_records_captured_inputs_and_units(self):
        p = Parameters(duration=.1)
        result = Result.from_dict(simulate(p, "export"), "export", p)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "samples.csv"
            result.export_csv(path)
            first, rest = path.read_text().split("\n", 1)
            metadata = json.loads(first[2:])
            self.assertEqual(metadata["parameters"], asdict(p))
            self.assertEqual(metadata["manifest"]["units"]["angle"], "rad")
            rows = list(csv.reader(io.StringIO(rest)))
            self.assertEqual(len(rows) - 1, len(result.samples))
            self.assertEqual(tuple(map(float, rows[-1])), result.samples[-1])


# A controlled child makes cancellation/failure ordering reproducible. It announces
# readiness only after installing its signal handler; the tests then release it.
FAULT_WORKER = r'''
import json, os, signal, sys, time
from pathlib import Path
directory = Path(sys.argv[1])
signal.signal(signal.SIGTERM, signal.SIG_IGN)
print("diagnostic" * 8000, flush=True)
(directory / "ready").write_text("ready")
while not (directory / "release").exists():
    time.sleep(.01)
mode = (directory / "release").read_text()
if mode == "crash":
    os._exit(7)
if mode == "bad":
    (directory / "result.json").write_text('{"schema_version": 1}')
    sys.exit(0)
(directory / "result.json").write_text(json.dumps({"status": "failed", "message": "Non-convergence simulée"}))
sys.exit(1)
'''


class ControlledController(Controller):
    def _worker_command(self, directory):
        self.test_directory = directory
        return sys.executable, ["-c", FAULT_WORKER, str(directory)]


class GuiRecipe(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])
        cls.app.setStyle("Fusion")

    def wait_until(self, predicate, timeout=10):
        deadline = time.monotonic() + timeout
        while not predicate() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(predicate(), "Événement attendu absent avant le délai.")

    def test_T00_GUI01_to_GUI08_real_window_recipe(self):
        window = MainWindow()
        try:
            window.show()
            QTest.qWait(30)
            self.assertEqual(window.parameters(), Parameters())
            self.assertIsNone(window.result)
            QTest.mouseClick(window.run_button, Qt.MouseButton.LeftButton)
            self.assertFalse(window.run_button.isEnabled())
            self.assertTrue(window.stop_button.isEnabled())
            window.fields["length"].setValue(1.5)  # edit after capture, before collection
            self.wait_until(lambda: window.controller.process is None)
            self.assertIsNotNone(window.result, window.status.text())
            self.assertEqual(window.result.parameters.length, 1.)
            self.assertEqual(window.pendulum.length, 1.)
            self.assertIn("différents", window.result_label.text())
            QTest.mouseClick(window.play_button, Qt.MouseButton.LeftButton)
            QTest.qWait(120)
            self.assertGreater(window.slider.value(), 0)
            QTest.mouseClick(window.play_button, Qt.MouseButton.LeftButton)
            self.assertFalse(window.timer.isActive())
            window.slider.setValue(150)
            self.assertEqual(window.pendulum.sample, window.result.samples[150])
            first_id = window.result.run_id
            QTest.mouseClick(window.run_button, Qt.MouseButton.LeftButton)
            self.wait_until(lambda: window.controller.process is None)
            self.assertNotEqual(window.result.run_id, first_id)
            self.assertEqual(window.result.parameters.length, 1.5)
            previous = window.result
            QTest.mouseClick(window.run_button, Qt.MouseButton.LeftButton)
            QTest.mouseClick(window.stop_button, Qt.MouseButton.LeftButton)
            self.wait_until(lambda: window.controller.process is None)
            self.assertIs(window.result, previous)
            self.assertIn("précédent", window.result_label.text())
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "params.json"
                window.save(path)
                window.fields["length"].setValue(2.)
                window.load(path)
                self.assertEqual(window.parameters().length, 1.5)
                write_json(path, {"format": "vinkulum-studio-parameters", "schema_version": 1,
                                   "parameters": asdict(Parameters(length=1.000001))})
                with self.assertRaises(ValueError):
                    window.load(path)
                self.assertEqual(window.parameters().length, 1.5)
            screenshot = os.environ.get("VINKULUM_STUDIO_SCREENSHOT")
            if screenshot:
                window.fields["length"].setValue(1.)
                window.run()
                self.wait_until(lambda: window.controller.process is None)
                window.slider.setValue(80)
                QTest.qWait(30)
                self.assertTrue(window.grab().save(screenshot))
        finally:
            window.close()

    def test_GUI05_controlled_failure_cancel_and_crash_keep_result(self):
        p = Parameters(duration=.1)
        previous = Result.from_dict(simulate(p, "previous"), "previous", p)
        for mode in ("failed", "bad", "crash", "cancel", "close"):
            with self.subTest(mode=mode), patch("vinkulum_studio.window.Controller", ControlledController):
                window = MainWindow()
                try:
                    window.show()
                    window.controller.last_result = previous
                    window._completed(previous)
                    ticks = []
                    heartbeat = QTimer(window)
                    heartbeat.setInterval(10)
                    heartbeat.timeout.connect(lambda: ticks.append(1))
                    heartbeat.start()
                    window.run()
                    directory = window.controller.test_directory
                    self.wait_until(lambda: (directory / "ready").exists())
                    window.fields["mass"].setValue(.5)
                    QTest.qWait(50)
                    if mode == "cancel":
                        window.stop()
                    elif mode == "close":
                        window.close()
                    elif mode == "crash":
                        window.controller.process.kill()
                    else:
                        (directory / "release").write_text(mode)
                    self.wait_until(lambda: window.controller.process is None)
                    self.assertIs(window.result, previous)
                    self.assertIs(window.controller.last_result, previous)
                    self.assertIn("précédent", window.result_label.text())
                    self.assertGreaterEqual(len(ticks), 3)
                    self.assertLessEqual(len(window.controller._log), 32_768)
                    self.assertFalse(directory.exists())
                    if mode == "failed":
                        self.assertIn("Non-convergence simulée", window.status.text())
                finally:
                    window.close()

    def test_GUI05_worker_unavailable_reenables_launch(self):
        class MissingController(Controller):
            def _worker_command(self, directory):
                return "/nonexistent/vinkulum-python", []
        with patch("vinkulum_studio.window.Controller", MissingController):
            window = MainWindow()
            try:
                window.run()
                self.wait_until(lambda: window.controller.process is None)
                self.assertIn("impossible", window.status.text())
                self.assertTrue(window.run_button.isEnabled())
                self.assertIsNone(window.result)
            finally:
                window.close()


if __name__ == "__main__":
    unittest.main()
