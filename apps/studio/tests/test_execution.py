"""CPU admission, cancellation races and real native worker provenance."""

import json
import os
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("VINKULUM_STUDIO_CPUS", "2")

from PySide6.QtCore import QProcess, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio.controller import Controller, ResultValidation
from vinkulum_studio.cpu_scheduler import CpuScheduler
from vinkulum_studio.execution import ExecutionPlan, default_budget
from vinkulum_studio.model import Parameters, Result, write_json
from vinkulum_studio.worker import simulate


class CpuContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def wait_until(self, predicate, timeout=15):
        deadline = time.monotonic() + timeout
        while not predicate() and time.monotonic() < deadline:
            QTest.qWait(5)
        self.assertTrue(predicate(), "CPU scheduling did not settle in time.")

    def test_cpu_inputs_and_captured_provenance(self):
        for data in (
            {"threads": True, "budget": 2},
            {"threads": 3, "budget": 2},
            {"threads": 0, "budget": 2},
            {"threads": 1},
            {"threads": 1, "budget": 2, "extra": 0},
        ):
            with self.subTest(data=data), self.assertRaises(ValueError):
                ExecutionPlan.from_dict(data)
        with patch("vinkulum_studio.execution.available_cpus", return_value=4):
            with patch.dict(os.environ, {"VINKULUM_STUDIO_CPUS": "3"}):
                self.assertEqual(default_budget(), 3)
            with (
                patch.dict(os.environ, {"VINKULUM_STUDIO_CPUS": "5"}),
                self.assertRaises(ValueError),
            ):
                default_budget()
        plan = ExecutionPlan(2, 2)
        output = simulate(Parameters(duration=0.01), "allocation", plan)
        plan.check_manifest(output["manifest"]["execution"])
        output["manifest"]["execution"]["threads"] = 1
        with self.assertRaises(ValueError):
            plan.check_manifest(output["manifest"]["execution"])

    def test_fifo_bound_release_and_cancel_before_dispatch(self):
        scheduler = CpuScheduler(3)
        observed = []
        first, second, third = [scheduler.request(n) for n in (2, 2, 1)]
        for name, ticket in zip(("first", "second", "third"), (first, second, third)):
            ticket.ready.connect(
                lambda plan, n=name: observed.append((n, scheduler.allocated))
            )
        QTest.qWait(20)
        self.assertEqual(observed, [("first", 2)])
        self.assertEqual(second.state, "queued")
        first.release()
        QTest.qWait(20)
        self.assertEqual(observed, [("first", 2), ("second", 2), ("third", 3)])
        second.release()
        third.release()
        self.assertEqual(scheduler.allocated, 0)
        cancelled = scheduler.request(3)
        cancelled.ready.connect(lambda _: self.fail("Cancelled job dispatched"))
        cancelled.release()
        cancelled.release()
        QTest.qWait(10)
        self.assertEqual(scheduler.allocated, 0)

    def test_synchronous_cancel_never_starts_a_process(self):
        scheduler = CpuScheduler(2)
        controller = Controller(scheduler=scheduler)
        busy = []
        controller.busy_changed.connect(busy.append)
        controller.busy_changed.connect(
            lambda state: controller.cancel() if state else None
        )
        with patch.object(
            QProcess, "start", side_effect=AssertionError("cancelled process started")
        ):
            controller.start(Parameters(duration=0.01))
            QTest.qWait(20)
        self.assertIsNone(controller.process)
        self.assertEqual(busy, [True, False])
        self.assertEqual(scheduler.allocated, 0)

    def test_real_workers_queue_capture_threads_and_release(self):
        scheduler = CpuScheduler(2)
        first, second = Controller(scheduler=scheduler), Controller(scheduler=scheduler)
        errors, stages = [], []
        for controller in (first, second):
            controller.problem.connect(errors.append)
        second.stage_changed.connect(stages.append)
        try:
            first.start(Parameters(duration=0.05), threads=2)
            second.start(Parameters(duration=0.02), threads=1)
            QTest.qWait(10)
            self.assertEqual(second._ticket.state, "queued")
            self.wait_until(lambda: first.process is None and second.process is None)
            self.assertEqual(errors, [])
            for controller, threads in ((first, 2), (second, 1)):
                manifest = json.loads(controller.last_result.manifest_json)
                self.assertEqual(manifest["execution"]["threads"], threads)
                self.assertEqual(manifest["execution"]["budget"], 2)
            self.assertEqual(scheduler.allocated, 0)
            self.assertTrue(any("1 native CPU" in stage for stage in stages))
        finally:
            first.shutdown()
            second.shutdown()

    def test_failed_start_releases_cpu_for_queued_calculation(self):
        scheduler = CpuScheduler(2)
        failed, good = Controller(scheduler=scheduler), Controller(scheduler=scheduler)
        errors = []
        failed.problem.connect(errors.append)
        failed._worker_command = lambda _: ("/nonexistent/vinkulum-worker", [])
        try:
            failed.start(Parameters(duration=0.01))
            good.start(Parameters(duration=0.01))
            self.wait_until(lambda: failed.process is None and good.process is None)
            self.assertTrue(errors)
            self.assertIsNotNone(good.last_result)
            self.assertEqual(scheduler.allocated, 0)
        finally:
            failed.shutdown()
            good.shutdown()

    def test_result_reader_runs_off_gui_thread_and_timer_remains_live(self):
        main_thread = threading.get_ident()
        reader_threads, ticks = [], []
        parameters, plan = Parameters(duration=0.01), ExecutionPlan(1, 2)
        output = simulate(parameters, "reader", plan)
        original = Result.from_dict

        def slow_reader(*args):
            reader_threads.append(threading.get_ident())
            time.sleep(0.15)
            return original(*args)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.json"
            write_json(path, output)
            reader = ResultValidation(path, "reader", parameters, plan, 0, None)
            timer = QTimer()
            timer.setInterval(5)
            timer.timeout.connect(lambda: ticks.append(1))
            with patch.object(Result, "from_dict", side_effect=slow_reader):
                timer.start()
                reader.start()
                self.wait_until(lambda: reader.isFinished())
                timer.stop()
                reader.wait()
            self.assertIsNotNone(reader.result, reader.message)
            self.assertNotEqual(reader_threads, [main_thread])
            self.assertGreater(len(ticks), 5)
            reader.deleteLater()

    def test_cancel_and_close_during_validation_preserve_previous_result_and_files(
        self,
    ):
        parameters = Parameters(duration=0.01)
        previous = Result.from_dict(
            simulate(parameters, "previous"), "previous", parameters
        )
        original = Result.from_dict
        for close in (False, True):
            with self.subTest(close=close):
                entered, release = threading.Event(), threading.Event()
                scheduler = CpuScheduler(2)
                controller = Controller(scheduler=scheduler)
                controller.last_result = previous
                published = []
                controller.completed.connect(published.append)

                def delayed_reader(*args, entered=entered, release=release):
                    entered.set()
                    if not release.wait(3):
                        raise ValueError("Reader test timed out.")
                    return original(*args)

                with patch.object(Result, "from_dict", side_effect=delayed_reader):
                    try:
                        controller.start(parameters)
                        directory = Path(controller._temporary.name)
                        self.wait_until(entered.is_set)
                        controller.cancel()
                        self.assertTrue(directory.is_dir())
                        self.assertEqual(scheduler.allocated, 2)
                        release.set()
                        if close:
                            controller.shutdown()
                        else:
                            self.wait_until(lambda c=controller: c.process is None)
                        self.assertIs(controller.last_result, previous)
                        self.assertEqual(published, [])
                        self.assertFalse(directory.exists())
                        self.assertEqual(scheduler.allocated, 0)
                    finally:
                        release.set()
                        controller.shutdown()

    def test_cancel_or_close_from_validation_stage_never_starts_reader(self):
        for close in (False, True):
            with self.subTest(close=close):
                scheduler = CpuScheduler(2)
                controller = Controller(scheduler=scheduler)
                problems = []
                controller.problem.connect(problems.append)

                def at_stage(message, c=controller, close=close):
                    if message.startswith("Checking captured"):
                        (c.shutdown if close else c.cancel)()

                controller.stage_changed.connect(at_stage)
                try:
                    with patch.object(
                        ResultValidation,
                        "start",
                        side_effect=AssertionError("cancelled reader started"),
                    ):
                        controller.start(Parameters(duration=0.01))
                        self.wait_until(lambda c=controller: c.process is None)
                    self.assertTrue(problems)
                    self.assertIsNone(controller.last_result)
                    self.assertEqual(scheduler.allocated, 0)
                finally:
                    controller.shutdown()


if __name__ == "__main__":
    unittest.main()
