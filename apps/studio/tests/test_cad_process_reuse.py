"""Real OCCT reuse with CPU admission, bounded retention and transaction isolation."""

import importlib.util
import json
import os
import sys
import tempfile
import textwrap
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QCoreApplication, QEvent, QProcess, QSettings, QTimer, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from shiboken6 import isValid

from vinkulum_studio.cad_controller import CadController
from vinkulum_studio.cad_process_cache import CadProcessCache, process_key
from vinkulum_studio.cpu_scheduler import CpuScheduler

BOX = {
    "operation": "box",
    "name": "Reuse plate",
    "dimensions_mm": [100, 60, 20],
    "density": 7800,
}
CAD = importlib.util.find_spec("build123d") is not None


@unittest.skipUnless(
    CAD and os.environ.get("VINKULUM_3D_TESTS"), "Requires desktop and OCCT"
)
class CadProcessReuse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.exceptions = self.enterContext(patch("sys.excepthook"))
        self.addCleanup(self.exceptions.assert_not_called)
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.scheduler = CpuScheduler(2)
        self.cache = CadProcessCache(idle_ms=5000)
        self.addCleanup(self.cache.shutdown)

    def tearDown(self):
        self.exceptions.assert_not_called()

    def controller(self):
        controller = CadController(scheduler=self.scheduler, cache=self.cache)
        self.addCleanup(controller.shutdown)
        return controller

    def wait(self, condition, seconds=15):
        deadline = time.monotonic() + seconds
        while not condition() and time.monotonic() < deadline:
            QTest.qWait(5)
        self.assertTrue(condition(), "CAD lifecycle did not settle")

    def run_request(self, controller, request=BOX, *, wrapper=None, **kwargs):
        controller.start(request, **kwargs)
        if wrapper is not None:
            # Instrument the real service before its deferred dispatch. Both
            # commands and cache identity include the actual wrapper path.
            controller._command = (sys.executable, [str(wrapper)])
            controller._cache_key = process_key(
                *controller._command,
                controller.process.processEnvironment(),
                controller._plan.threads,
            )

    def success(self, controller):
        self.wait(lambda: controller.process is None)
        self.assertIsNone(controller.failure_kind, controller.last_log)
        self.assertEqual(self.scheduler.allocated, 0)
        self.assertIsNotNone(self.cache.idle_process)
        self.assertEqual(self.cache.idle_process.state(), QProcess.ProcessState.Running)
        return self.cache.idle_process.processId()

    def wrapper(self, name, code):
        path = self.root / f"{name}.py"
        path.write_text(textwrap.dedent(code))
        return path

    def paused_service(self):
        return self.wrapper(
            "pause",
            """
            import json, pathlib, sys, time
            from vinkulum_studio import cad_worker
            original = cad_worker.evaluate
            def evaluate(source, output, **kwargs):
                code = original(source, output, **kwargs)
                if json.loads(pathlib.Path(source).read_text())["operation"] == "fillet":
                    pathlib.Path(__file__).with_suffix(".entered").touch()
                    time.sleep(120)
                return code
            cad_worker.evaluate = evaluate
            raise SystemExit(cad_worker.service())
        """,
        )

    def test_two_controllers_reuse_one_worker_and_keep_distinct_transactions(self):
        first = self.controller()
        self.run_request(first)
        pid = self.success(first)
        result = first.last_result
        first.shutdown()
        second = self.controller()
        self.run_request(second, {**BOX, "dimensions_mm": [10, 20, 30]})
        self.assertEqual(self.success(second), pid)
        self.assertAlmostEqual(second.last_result["body"]["mass"], 0.0468, places=13)
        self.assertNotEqual(
            result["service_request_id"], second.last_result["service_request_id"]
        )
        self.assertNotEqual(
            result["request_sha256"], second.last_result["request_sha256"]
        )
        self.assertEqual(first.last_result, result)
        self.assertEqual(self.cache.idle_process.parent(), self.cache)

    def test_editor_dialogs_are_destroyed_while_the_shared_worker_survives(self):
        from vinkulum_studio.cad_dialog import CadDialog
        from vinkulum_studio.cad_history_dialog import CadHistoryDialog
        from vinkulum_studio.cad_process_cache import application_cad_cache
        from vinkulum_studio.editor import EditorWindow

        cache = application_cad_cache()
        cache.expire()
        self.addCleanup(cache.expire)
        window = EditorWindow(
            QSettings(str(self.root / "window.ini"), QSettings.Format.IniFormat)
        )
        window._discard_allowed = lambda: True
        self.addCleanup(window.close)
        window.new_project()
        window.show()
        pids, errors = [], []
        for action in ("cancel", "create", "history"):

            def interact():
                dialog = self.app.activeModalWidget()
                try:
                    if action == "cancel":
                        dialog.reject()
                    elif action == "create":
                        self.assertIsInstance(dialog, CadDialog)
                        dialog.start()
                        self.wait(lambda: dialog.controller.process is None)
                    else:
                        self.assertIsInstance(dialog, CadHistoryDialog)
                        field = dialog.dimensions[0]
                        field.setFocus()
                        QTest.keyClick(
                            field, Qt.Key.Key_A, Qt.KeyboardModifier.ControlModifier
                        )
                        QTest.keyClicks(field, "110")
                        dialog.start_preview()
                        self.wait(
                            lambda: dialog.controller.process is None
                            and dialog.preview_body is not None
                        )
                        dialog.accept()
                except BaseException as error:
                    errors.append(error)
                    dialog.reject()

            QTimer.singleShot(0, interact)
            if action == "history":
                window.open_cad_history()
            else:
                window.open_cad()
            self.assertFalse(errors, errors)
            QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
            self.assertEqual(window.findChildren(CadDialog), [])
            self.assertEqual(window.findChildren(CadHistoryDialog), [])
            if action != "cancel":
                pids.append(cache.idle_process.processId())
        self.assertEqual(pids[0], pids[1])
        self.assertEqual(len(window.project.bodies), 1)
        self.assertAlmostEqual(window.project.bodies[0].dimensions[0], 0.11, places=12)

    def test_pool_size_and_environment_changes_replace_the_idle_process(self):
        controller = self.controller()
        self.run_request(controller, threads=2)
        first = self.success(controller)
        self.run_request(controller, threads=1)
        second = self.success(controller)
        self.assertNotEqual(first, second)
        self.assertEqual(controller.last_result["execution"]["threads"], 1)
        with patch.dict(os.environ, {"VINKULUM_CAD_CACHE_TEST": "changed"}):
            self.run_request(controller, threads=1)
            self.assertNotEqual(second, self.success(controller))

    def test_concurrent_admission_runs_two_workers_and_retains_only_one(self):
        first, second = self.controller(), self.controller()
        pids = []
        first.completed.connect(
            lambda result: pids.append(self.cache.idle_process.processId())
        )
        second.completed.connect(
            lambda result: pids.append(self.cache.idle_process.processId())
        )
        self.run_request(first, threads=1)
        self.run_request(second, {**BOX, "dimensions_mm": [10, 20, 30]}, threads=1)
        self.wait(
            lambda: first.process is not None
            and second.process is not None
            and first.process.processId() > 0
            and second.process.processId() > 0
        )
        self.assertEqual(self.scheduler.allocated, 2)
        self.wait(lambda: first.process is None and second.process is None)
        self.assertIsNone(first.failure_kind, first.last_log)
        self.assertIsNone(second.failure_kind, second.last_log)
        self.assertEqual(len(set(pids)), 2)
        self.assertAlmostEqual(first.last_result["body"]["mass"], 0.936, places=12)
        self.assertAlmostEqual(second.last_result["body"]["mass"], 0.0468, places=12)
        self.assertEqual(self.scheduler.allocated, 0)
        self.wait(lambda: not self.cache._retiring)
        self.assertIsNotNone(self.cache.idle_process)

    def test_idle_expiry_reaps_the_process_without_reserving_cpu(self):
        self.cache.idle_ms = 30
        controller = self.controller()
        self.run_request(controller)
        self.success(controller)
        process = self.cache.idle_process
        self.assertEqual(self.scheduler.allocated, 0)
        self.wait(lambda: self.cache.idle_process is None and not self.cache._retiring)
        self.wait(lambda: not isValid(process))

    def test_queued_cancellation_preserves_the_idle_worker_and_fifo(self):
        controller = self.controller()
        self.run_request(controller)
        pid = self.success(controller)
        previous = controller.last_result
        blocker = self.scheduler.request(2)
        self.addCleanup(blocker.release)
        self.wait(lambda: self.scheduler.allocated == 2)
        self.run_request(controller)
        self.assertEqual(controller.process.state(), QProcess.ProcessState.NotRunning)
        self.assertFalse(controller.timer.isActive())
        following = self.scheduler.request(2)
        self.addCleanup(following.release)
        controller.cancel()
        self.assertIsNone(controller.process)
        self.assertEqual(self.cache.idle_process.processId(), pid)
        self.assertIs(controller.last_result, previous)
        blocker.release()
        self.wait(lambda: following.state == "running")
        following.release()
        self.run_request(controller)
        self.assertEqual(self.success(controller), pid)

    def test_native_operation_error_retires_the_worker_before_releasing_its_lease(self):
        controller = self.controller()
        self.run_request(controller)
        pid = self.success(controller)
        previous = controller.last_result
        request = {
            "operation": "fillet",
            "a": previous["body"],
            "radius_mm": 1000,
            "density": 7800,
        }
        self.run_request(controller, request)
        self.wait(lambda: controller._acknowledged)
        if controller.process is not None:
            self.assertEqual(self.scheduler.allocated, 2)
        self.wait(lambda: controller.process is None)
        self.assertEqual(controller.failure_kind, "operation_failed")
        self.assertIs(controller.last_result, previous)
        self.assertIsNone(self.cache.idle_process)
        self.assertEqual(self.scheduler.allocated, 0)
        self.run_request(controller)
        self.assertNotEqual(pid, self.success(controller))

    def test_cancellation_discards_a_real_unpublished_result(self):
        self._fault_case("cancel")

    def test_cancel_during_validation_keeps_lease_and_files_until_reader_finishes(self):
        self._validation_fault("cancel")

    def test_process_death_during_validation_discards_the_prepared_document(self):
        self._validation_fault("death")

    def test_partial_extra_response_during_validation_is_rejected(self):
        self._validation_fault("extra")

    def test_nonzero_exit_during_validation_discards_the_prepared_document(self):
        self._validation_fault("error_exit")

    def _validation_fault(self, fault):
        from vinkulum_studio.cad_admission import admit_cad_response

        wrapper = None
        if fault in ("extra", "error_exit"):
            action = (
                'os.write(protocol, b"incomplete extra response")'
                if fault == "extra"
                else "os._exit(7)"
            )
            wrapper = self.wrapper(
                fault,
                f"""
                import json, os, pathlib, threading, time
                from vinkulum_studio import cad_worker
                protocol = os.dup(1)
                original = cad_worker.evaluate
                def noise():
                    marker = pathlib.Path(__file__).with_suffix(".noise")
                    while not marker.exists():
                        time.sleep(.005)
                    {action}
                def evaluate(source, output, **kwargs):
                    code = original(source, output, **kwargs)
                    if json.loads(pathlib.Path(source).read_text())["dimensions_mm"][0] == 10:
                        threading.Thread(target=noise, daemon=True).start()
                    return code
                cad_worker.evaluate = evaluate
                raise SystemExit(cad_worker.service())
            """,
            )
        controller = self.controller()
        self.run_request(controller, wrapper=wrapper)
        pid = self.success(controller)
        previous = controller.last_result
        entered, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)
        published = []
        controller.completed.connect(published.append)

        def admit(*args):
            result = admit_cad_response(*args)
            entered.set()
            if not release.wait(5):
                raise RuntimeError("Test did not release the CAD validator")
            return result

        with patch(
            "vinkulum_studio.cad_controller.admit_cad_response", side_effect=admit
        ):
            self.run_request(
                controller, {**BOX, "dimensions_mm": [10, 20, 30]}, wrapper=wrapper
            )
            self.wait(entered.is_set)
            process = controller.process
            temporary = Path(controller.temporary.name)
            self.assertEqual(process.processId(), pid)
            if fault == "cancel":
                controller.cancel()
            elif fault in ("extra", "error_exit"):
                wrapper.with_suffix(".noise").touch()
            else:
                process.kill()
            self.wait(lambda: process.state() == QProcess.ProcessState.NotRunning)
            if fault == "error_exit":
                self.assertEqual(process.exitCode(), 7)
                self.assertEqual(process.exitStatus(), QProcess.ExitStatus.NormalExit)
            self.assertEqual(self.scheduler.allocated, 2)
            self.assertTrue(temporary.exists())
            self.assertIs(controller.process, process)
            self.assertIsNone(self.cache.idle_process)
            self.assertEqual(published, [])
            release.set()
            self.wait(lambda: controller.process is None)
        self.assertEqual(
            controller.failure_kind,
            {
                "cancel": "cancelled",
                "death": "crashed",
                "extra": "invalid_result",
                "error_exit": "worker_failed",
            }[fault],
        )
        self.assertIs(controller.last_result, previous)
        self.assertFalse(temporary.exists())
        self.assertEqual(self.scheduler.allocated, 0)
        self.assertEqual(published, [])
        self.run_request(controller)
        self.assertNotEqual(self.success(controller), pid)

    def test_timeout_discards_a_real_unpublished_result(self):
        self._fault_case("timeout")

    def test_process_death_discards_a_real_unpublished_result(self):
        self._fault_case("death")

    def _fault_case(self, fault):
        self.cache.expire()
        wrapper = self.paused_service()
        marker = wrapper.with_suffix(".entered")
        marker.unlink(missing_ok=True)
        controller = self.controller()
        self.run_request(controller, wrapper=wrapper)
        pid = self.success(controller)
        previous = controller.last_result
        published, ticks = [], []
        controller.completed.connect(published.append)
        request = {
            "operation": "fillet",
            "a": previous["body"],
            "radius_mm": 1,
            "density": 7800,
        }
        self.run_request(
            controller,
            request,
            wrapper=wrapper,
            timeout=1 if fault == "timeout" else 10,
        )
        self.wait(marker.exists)
        self.assertEqual(controller.process.processId(), pid)
        self.assertTrue((Path(controller.temporary.name) / "output.json").is_file())
        self.assertIs(controller.last_result, previous)
        self.assertEqual(self.scheduler.allocated, 2)
        timer = QTimer()
        timer.setInterval(5)
        timer.timeout.connect(lambda: ticks.append(True))
        timer.start()
        try:
            self.wait(lambda: len(ticks) >= 3, seconds=0.3)
        finally:
            timer.stop()
        if fault == "cancel":
            controller.cancel()
            self.assertEqual(self.scheduler.allocated, 2)
        elif fault == "death":
            controller.process.kill()
        self.wait(lambda: controller.process is None)
        self.assertEqual(
            controller.failure_kind,
            {"cancel": "cancelled", "timeout": "timed_out", "death": "crashed"}[fault],
        )
        self.assertEqual(published, [])
        self.assertIs(controller.last_result, previous)
        self.assertIsNone(self.cache.idle_process)
        self.assertEqual(self.scheduler.allocated, 0)
        self.run_request(controller)
        self.assertNotEqual(pid, self.success(controller))

    def test_stale_acknowledgement_and_stale_result_identity_are_rejected(self):
        for location in ("ack", "file"):
            with self.subTest(location=location):
                self.cache.expire()
                wrapper = self.wrapper(
                    "stale",
                    """
                    import json, pathlib, sys
                    from vinkulum_studio import cad_worker
                    original = cad_worker.evaluate
                    def evaluate(source, output, **kwargs):
                        code = original(source, output, **kwargs)
                        result = json.loads(pathlib.Path(output).read_text())
                        result["service_request_id"] = "0" * 32
                        pathlib.Path(output).write_text(json.dumps(result))
                        return code
                    cad_worker.evaluate = evaluate
                    raise SystemExit(cad_worker.service())
                """,
                )
                if location == "ack":
                    wrapper.write_text(textwrap.dedent("""
                        import json, sys
                        from vinkulum_studio import cad_worker
                        command = json.loads(sys.stdin.readline())
                        code = cad_worker.evaluate(command["input"], command["output"], request_id=command["id"])
                        print(json.dumps({"protocol": 1, "id": "0" * 32, "code": code}), flush=True)
                        import time
                        time.sleep(120)
                    """))
                controller = self.controller()
                self.run_request(controller, wrapper=wrapper)
                self.wait(lambda: controller.process is None)
                self.assertEqual(
                    controller.failure_kind, "invalid_result", controller.last_log
                )
                self.assertIsNone(controller.last_result)
                self.assertIsNone(self.cache.idle_process)
                self.assertEqual(self.scheduler.allocated, 0)


if __name__ == "__main__":
    unittest.main()
