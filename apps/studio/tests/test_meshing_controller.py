"""Direct-process cancellation and threaded mesh admission through real Qt."""

import os
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from test_meshing import request, tetra_msh
from vinkulum_studio.meshing import finish_mesh, load_mesh
from vinkulum_studio.meshing_controller import MeshingController


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires a desktop Qt context"
)
class MeshingLifecycle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.controller = MeshingController()
        self.addCleanup(self.controller.shutdown)
        self.exceptions = patch("sys.excepthook").start()
        self.addCleanup(patch.stopall)
        self.addCleanup(self.exceptions.assert_not_called)

    def wait_until(self, condition, seconds=10):
        deadline = time.monotonic() + seconds
        while not condition() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(condition(), "Mesh worker did not reach the expected state.")

    def executable(self, name, body, occ="8.0.1"):
        path = self.root / name
        path.write_text(
            f"#!{sys.executable}\nimport sys, os, time\nfrom pathlib import Path\n"
            f'if "-info" in sys.argv:\n print("Version : 5.0.0-test\\nOCC version : {occ}", flush=True)\n sys.exit(0)\n{body}\n'
        )
        path.chmod(0o755)
        return path

    def successful_engine(self):
        return self.executable(
            "fixture",
            f'Path("mesh.msh").write_bytes({tetra_msh()!r})\nprint("Synthetic controller fixture", flush=True)',
        )

    def test_success_reopens_then_failure_preserves_previous_mesh(self):
        controller = self.controller
        output = self.root / "successful"
        controller.start(request(), output, executable=self.successful_engine())
        self.wait_until(lambda: not controller.busy)
        self.assertIsNone(controller.failure_kind)
        previous = controller.last_result
        self.assertEqual(load_mesh(output)[0], previous[0])
        cases = (
            ("old-occ", "", "7.6.3", "invalid_version"),
            ("failed", "sys.exit(2)", "8.0.1", "mesher_failed"),
            ("invalid", "print('Error : deliberate failure')", "8.0.1", "invalid_mesh"),
            (
                "crash",
                "import signal\nos.kill(os.getpid(), signal.SIGKILL)",
                "8.0.1",
                "crashed",
            ),
        )
        for name, body, occ, expected in cases:
            executable = self.executable(name, body, occ)
            root = self.root / (name + "-output")
            controller.start(request(), root, executable=executable)
            self.wait_until(lambda: not controller.busy)
            self.assertEqual(controller.failure_kind, expected)
            self.assertIs(controller.last_result, previous)
            self.assertFalse((root / "result.json").exists())

    def test_cancellation_kills_direct_process_and_timeout_is_distinct(self):
        controller = self.controller
        engine = self.executable(
            "hanging", 'Path("pid.txt").write_text(str(os.getpid()))\ntime.sleep(30)'
        )
        root = self.root / "cancelled"
        controller.start(request(), root, executable=engine)
        self.wait_until(lambda: (root / "pid.txt").exists())
        pid = int((root / "pid.txt").read_text())
        controller.cancel()
        self.wait_until(lambda: not controller.busy)
        self.assertEqual(controller.failure_kind, "cancelled")
        with self.assertRaises(ProcessLookupError):
            os.kill(pid, 0)
        controller.start(
            request(), self.root / "timeout", executable=engine, timeout=0.05
        )
        self.wait_until(lambda: not controller.busy)
        self.assertEqual(controller.failure_kind, "timed_out")

    def test_cancellation_from_stage_signal_never_starts_the_engine(self):
        controller = self.controller
        marker = self.root / "must-not-be-written"
        executable = self.root / "never-start"
        executable.write_text(
            f'#!{sys.executable}\nfrom pathlib import Path\nPath({str(marker)!r}).write_text("bad")\n'
        )
        executable.chmod(0o755)
        controller.stage_changed.connect(lambda _: controller.cancel())
        controller.start(request(), self.root / "cancelled", executable=executable)
        self.assertFalse(controller.busy)
        self.assertEqual(controller.failure_kind, "cancelled")
        QTest.qWait(30)
        self.assertFalse(marker.exists())

    def test_validation_is_off_thread_and_cancellable_without_publication(self):
        started = threading.Event()
        heartbeats = []
        timer = QTimer()
        timer.timeout.connect(lambda: heartbeats.append(1))
        timer.start(5)
        self.addCleanup(timer.stop)

        def slow_validation(run, code, *, publish, cancellation):
            started.set()
            deadline = time.monotonic() + 5
            while not cancellation() and time.monotonic() < deadline:
                time.sleep(0.005)
            if cancellation():
                raise InterruptedError("Cancelled during numerical admission.")
            return finish_mesh(run, code, publish=publish, cancellation=cancellation)

        controller = self.controller
        root = self.root / "validation"
        with patch(
            "vinkulum_studio.meshing_controller.finish_mesh",
            side_effect=slow_validation,
        ):
            controller.start(request(), root, executable=self.successful_engine())
            self.wait_until(started.is_set)
            before = len(heartbeats)
            self.wait_until(lambda: len(heartbeats) > before + 2)
            controller.cancel()
            self.wait_until(lambda: not controller.busy)
        self.assertEqual(controller.failure_kind, "cancelled")
        self.assertFalse((root / "result.json").exists())
        self.assertIsNone(controller.last_result)


if __name__ == "__main__":
    unittest.main()
