"""CPU allocation contracts, real solver equivalence and cross-engine admission."""

import json
import os
import shutil
import sys
import tempfile
import threading
import time
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from test_meshing import request
from vinkulum_studio.calculix import load_static_result, run_static
from vinkulum_studio.engine_threads import calculix_threads, thread_environment
from vinkulum_studio.execution import ExecutionPlan
from vinkulum_studio.meshing import HxtMeshRequest, MeshRequest, parse_info
from vinkulum_studio.static_window import tension_example


class NativeEngineContracts(unittest.TestCase):
    def test_native_policies_override_inherited_oversubscription(self):
        for engine in ("occt", "gmsh", "calculix", "pinocchio"):
            values = thread_environment(2, engine=engine)
            self.assertEqual(values["OMP_NUM_THREADS"], "2")
            self.assertEqual(values["OMP_THREAD_LIMIT"], "2")
            self.assertEqual(values["OMP_MAX_ACTIVE_LEVELS"], "1")
            self.assertEqual(values["OPENBLAS_NUM_THREADS"], "1")
            self.assertEqual(values["MKL_NUM_THREADS"], "1")
        ccx = thread_environment(4, engine="calculix")
        for key in (
            "NUMBER_OF_CPUS",
            "CCX_NPROC_STIFFNESS",
            "CCX_NPROC_RESULTS",
            "CCX_NPROC_EQUATION_SOLVER",
        ):
            self.assertEqual(ccx[key], "4")
        for value in (True, 0, -1, 1.0, "2"):
            with self.assertRaises((ValueError, TypeError)):
                thread_environment(value, engine="calculix")
        for text in ("", "Using up to 3 cpu(s) for spooles."):
            with self.assertRaises(ValueError):
                calculix_threads(text, ExecutionPlan(2, 2))

    def test_hxt_requires_parallel_build_and_legacy_request_remains_exact(self):
        legacy = request()
        self.assertEqual(set(asdict(legacy)), {"body", "size_mm", "order", "threads"})
        self.assertIn("Mesh.Algorithm3D = 1;", legacy.geo_script())
        self.assertEqual(MeshRequest.from_dict(asdict(legacy)), legacy)
        hxt = HxtMeshRequest(**asdict(legacy) | {"body": legacy.body})
        self.assertEqual(MeshRequest.from_dict(asdict(hxt)), hxt)
        self.assertIn("Mesh.Algorithm3D = 10;", hxt.geo_script())
        info = "Version : 5.0\nOCC version : 8.0.1\n"
        self.assertEqual(parse_info(info, legacy), ("5.0", "8.0.1"))
        for options in ("", "Hxt", "OpenMP"):
            with self.assertRaises(ValueError):
                parse_info(info + "Build options : " + options, hxt)
        self.assertEqual(
            parse_info(info + "Build options : Hxt OpenMP", hxt), ("5.0", "8.0.1")
        )

    @unittest.skipUnless(shutil.which("ccx"), "Requires the real CalculiX executable")
    def test_real_calculix_one_two_four_threads_preserve_physics_and_strict_archive(
        self,
    ):
        with (
            tempfile.TemporaryDirectory() as directory,
            patch.dict(
                os.environ, {"NUMBER_OF_CPUS": "128", "OPENBLAS_NUM_THREADS": "128"}
            ),
        ):
            reports = []
            for threads in (1, 2, 4):
                root = Path(directory) / str(threads)
                result = run_static(tension_example(), root, threads=threads)
                self.assertEqual(result["requested_threads"], threads)
                phases = result["execution"]["reported_phases"]
                self.assertTrue(phases)
                self.assertTrue(all(1 <= p["threads"] <= threads for p in phases))
                self.assertEqual(load_static_result(root)[1], result)
                reports.append(result)
            for result in reports[1:]:
                for key in ("displacements", "stress", "reactions", "strain_energy_J"):
                    np.testing.assert_allclose(
                        result[key], reports[0][key], rtol=1e-12, atol=1e-15
                    )
            root = Path(directory) / "1"
            path = root / "result.json"
            report = json.loads(path.read_text())
            report["execution"]["allocation"]["threads"] = True
            path.write_text(json.dumps(report))
            with self.assertRaises(ValueError):
                load_static_result(root)


@unittest.skipUnless(os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires Qt")
class CrossEngineAdmission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        from vinkulum_studio.cpu_scheduler import CpuScheduler

        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.scheduler = CpuScheduler(2)
        caught = patch("sys.excepthook")
        self.errors = caught.start()
        self.addCleanup(caught.stop)
        self.addCleanup(self.errors.assert_not_called)

    def wait_for(self, predicate):
        deadline = time.monotonic() + 10
        while not predicate() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(predicate(), "CPU/worker lifecycle did not settle")

    def engine(self, name):
        path = self.root / name
        path.write_text(
            f"#!{sys.executable}\nimport os,json\nfrom pathlib import Path\nPath({str(self.root / (name + '.json'))!r}).write_text(json.dumps(dict(os.environ)))\nraise SystemExit(2)\n"
        )
        path.chmod(0o755)
        return path

    def controllers(self):
        from vinkulum_studio.cad_controller import CadController
        from vinkulum_studio.examples3d import double_pendulum
        from vinkulum_studio.meshing_controller import MeshingController
        from vinkulum_studio.pinocchio_controller import PinocchioController
        from vinkulum_studio.static_controller import StaticController

        cad = CadController(scheduler=self.scheduler)
        mesh = MeshingController(scheduler=self.scheduler)
        ccx = StaticController(scheduler=self.scheduler)
        pin = PinocchioController(scheduler=self.scheduler)
        jobs = (
            (
                cad,
                lambda: cad.start(
                    {"operation": "self_check"},
                    executable=self.engine("cad"),
                    arguments=[],
                    threads=1,
                ),
            ),
            (
                mesh,
                lambda: mesh.start(
                    replace(request(), threads=1),
                    self.root / "mesh-run",
                    executable=self.engine("mesh"),
                ),
            ),
            (
                ccx,
                lambda: ccx.start(
                    tension_example(),
                    self.root / "ccx-run",
                    executable=self.engine("ccx"),
                    threads=1,
                ),
            ),
            (
                pin,
                lambda: pin.start(
                    double_pendulum(),
                    {},
                    self.root / "pin-run",
                    interpreter=self.engine("pin"),
                ),
            ),
        )
        for controller, _ in jobs:
            controller.last_result = self
            self.addCleanup(controller.shutdown)
        return jobs

    def test_queue_cancellation_and_failures_release_a_shared_budget(self):
        blocker = self.scheduler.request(2)
        self.wait_for(lambda: self.scheduler.allocated == 2)
        jobs = self.controllers()
        for controller, start in jobs:
            start()
            self.assertFalse(controller.timer.isActive())
        QTest.qWait(30)
        self.assertFalse(list(self.root.glob("*.json")))
        jobs[0][0].cancel()
        self.assertIsNone(jobs[0][0].process)
        blocker.release()
        self.wait_for(
            lambda: (
                self.scheduler.allocated == 0
                and all(c._ticket is None for c, _ in jobs)
            )
        )
        self.assertFalse((self.root / "cad.json").exists())
        for name in ("mesh", "ccx", "pin"):
            environment = json.loads((self.root / (name + ".json")).read_text())
            self.assertEqual(environment["OMP_NUM_THREADS"], "1")
            self.assertEqual(environment["OPENBLAS_NUM_THREADS"], "1")
        for controller, _ in jobs:
            self.assertIs(controller.last_result, self)

    def test_cancelling_from_queue_signal_never_starts_any_engine(self):
        for controller, start in self.controllers():
            controller.stage_changed.connect(controller.cancel)
            start()
        QTest.qWait(50)
        self.assertEqual(self.scheduler.allocated, 0)
        self.assertFalse(list(self.root.glob("*.json")))

    @unittest.skipUnless(shutil.which("ccx"), "Requires real CalculiX")
    def test_static_checks_keep_gui_live_and_lease_until_cancelled_reader_finishes(
        self,
    ):
        from vinkulum_studio.static_controller import StaticController, finish_run

        entered, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)
        identities = []

        def validate(*args, **kwargs):
            identities.append(threading.get_ident())
            entered.set()
            if not release.wait(5):
                raise RuntimeError("Test did not release validator")
            return finish_run(*args, **kwargs)

        controller = StaticController(scheduler=self.scheduler)
        controller.last_result = self
        self.addCleanup(controller.shutdown)
        self.addCleanup(release.set)
        ticks = []
        heartbeat = QTimer()
        heartbeat.timeout.connect(lambda: ticks.append(True))
        heartbeat.start(5)
        with patch(
            "vinkulum_studio.static_controller.finish_run", side_effect=validate
        ):
            controller.start(tension_example(), self.root / "static", threads=2)
            self.wait_for(entered.is_set)
            self.assertNotEqual(identities, [threading.get_ident()])
            self.assertEqual(self.scheduler.allocated, 2)
            before = len(ticks)
            QTest.qWait(40)
            self.assertGreater(len(ticks), before)
            controller.cancel()
            self.assertEqual(self.scheduler.allocated, 2)
            release.set()
            self.wait_for(lambda: controller.process is None)
        heartbeat.stop()
        self.assertEqual(self.scheduler.allocated, 0)
        self.assertIs(controller.last_result, self)
        self.assertFalse((self.root / "static/result.json").exists())
        self.assertFalse((self.root / "static/validated-result.json").exists())

    def test_pinocchio_reader_holds_lease_and_cancellation_preserves_previous_result(
        self,
    ):
        from vinkulum_studio.examples3d import double_pendulum
        from vinkulum_studio.pinocchio_controller import PinocchioController

        interpreter = self.root / "operator-fixture"
        interpreter.write_text(f"#!{sys.executable}\nraise SystemExit(0)\n")
        interpreter.chmod(0o755)
        entered, release = threading.Event(), threading.Event()
        identities = []

        def validate(*args, **kwargs):
            identities.append(threading.get_ident())
            entered.set()
            if not release.wait(5):
                raise RuntimeError("Test did not release operator reader")
            return object()

        controller = PinocchioController(scheduler=self.scheduler)
        controller.last_result = self
        self.addCleanup(controller.shutdown)
        self.addCleanup(release.set)
        with patch(
            "vinkulum_studio.pinocchio_controller.load_operators", side_effect=validate
        ):
            controller.start(
                double_pendulum(), {}, self.root / "operators", interpreter=interpreter
            )
            self.wait_for(entered.is_set)
            self.assertNotEqual(identities, [threading.get_ident()])
            self.assertEqual(self.scheduler.allocated, 1)
            controller.cancel()
            self.assertEqual(self.scheduler.allocated, 1)
            release.set()
            self.wait_for(lambda: controller.process is None)
        self.assertIs(controller.last_result, self)
        self.assertEqual(self.scheduler.allocated, 0)

    @unittest.skipUnless(shutil.which("ccx"), "Requires real CalculiX")
    def test_static_window_close_defers_destruction_during_result_checks(self):
        from shiboken6 import isValid
        from vinkulum_studio.static_controller import finish_run
        from vinkulum_studio.static_window import StaticWindow

        entered, release = threading.Event(), threading.Event()

        def validate(*args, **kwargs):
            entered.set()
            if not release.wait(5):
                raise RuntimeError("Test did not release result checks")
            return finish_run(*args, **kwargs)

        window = StaticWindow()
        window._discard_allowed = lambda: True
        self.addCleanup(lambda: window.close() if isValid(window) else None)
        self.addCleanup(release.set)
        with patch(
            "vinkulum_studio.static_controller.finish_run", side_effect=validate
        ):
            window.controller.start(tension_example(), self.root / "closing", threads=1)
            self.wait_for(entered.is_set)
            started = time.monotonic()
            window.close()
            self.assertLess(time.monotonic() - started, 0.25)
            self.assertTrue(isValid(window))
            self.assertTrue(window._close_pending)
            release.set()
            self.wait_for(lambda: not isValid(window))
        self.assertFalse((self.root / "closing/result.json").exists())


if __name__ == "__main__":
    unittest.main()
