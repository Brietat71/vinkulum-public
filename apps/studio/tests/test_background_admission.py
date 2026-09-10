"""Real Qt cancellation and publication boundaries for CAD and stored results."""

import importlib.util
import os
import sys
import tempfile
import threading
import time
import unittest
import zipfile
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QDialog, QFileDialog
from shiboken6 import isValid
from vinkulum_studio.archive_controller import ArchiveController
from vinkulum_studio.cpu_scheduler import CpuScheduler


@unittest.skipUnless(os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires Qt")
class BackgroundAdmission(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        caught = patch("sys.excepthook")
        self.errors = caught.start()
        self.addCleanup(caught.stop)
        self.addCleanup(self.check_errors)

    def check_errors(self):
        QTest.qWait(10)
        self.errors.assert_not_called()

    def wait_for(self, predicate, seconds=15):
        deadline = time.monotonic() + seconds
        while not predicate() and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertTrue(
            predicate(), "Background operation did not reach its expected state"
        )

    def test_archive_queue_cancel_and_close_from_stage_never_read_files(self):
        scheduler = CpuScheduler(1)
        blocker = scheduler.request(1)
        self.wait_for(lambda: scheduler.allocated == 1)
        controller = ArchiveController(scheduler=scheduler)
        self.addCleanup(controller.shutdown)

        def read(path):
            with open(path) as stream:
                return stream.read()

        with patch(
            "builtins.open",
            side_effect=AssertionError("Queued archive must not read files"),
        ):
            controller.start(read, "not-read")
            QTest.qWait(20)
            self.assertIsNone(controller.reader)
            controller.cancel()
        self.assertFalse(controller.busy)
        self.assertEqual(scheduler.allocated, 1)
        blocker.release()
        called = []

        def close_at_read(message):
            if message.startswith("Reading archive"):
                controller.shutdown()

        controller.stage_changed.connect(close_at_read)
        controller.start(lambda: called.append(True))
        self.wait_for(lambda: not controller.busy)
        self.assertFalse(called)
        self.assertEqual(scheduler.allocated, 0)

    def test_static_archive_reader_keeps_gui_alive_and_cancel_preserves_inputs(self):
        from vinkulum_studio.calculix import load_static_result
        from vinkulum_studio.static_window import StaticWindow

        source = (
            Path(__file__).resolve().parents[3]
            / "docs/bancs/studio-c3d8-refinement-2026/mesh-20-2-1.zip"
        )
        archive = self.root / "archive"
        archive.mkdir()
        with zipfile.ZipFile(source) as z:
            for name in z.namelist():
                (archive / name).write_bytes(z.read(name))
        previous = load_static_result(archive)
        window = StaticWindow()
        self.addCleanup(window.close)
        window._discard_allowed = lambda: True
        window.show()
        window._present_result(previous, archived=True)
        window.young.setText("123456789.01234567")
        pending = window.edited_study()
        positions = window.viewport.display_positions.copy()
        entered, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)
        identities, ticks = [], []

        def load(path):
            identities.append(threading.get_ident())
            entered.set()
            if not release.wait(5):
                raise RuntimeError("Test did not release archive reader")
            return load_static_result(path)

        heartbeat = QTimer()
        heartbeat.timeout.connect(lambda: ticks.append(True))
        heartbeat.start(5)
        with (
            patch("vinkulum_studio.static_window.load_static_result", side_effect=load),
            patch.object(
                QFileDialog,
                "getOpenFileName",
                return_value=(str(archive / "result.json"), ""),
            ),
        ):
            window.open_result_dialog()
            self.wait_for(entered.is_set)
            self.assertNotEqual(identities, [threading.get_ident()])
            before = len(ticks)
            QTest.qWait(40)
            self.assertGreater(len(ticks), before)
            window.cancel_work()
            self.assertTrue(window.archive.busy)
            release.set()
            self.wait_for(lambda: not window.archive.busy)
        heartbeat.stop()
        self.assertIs(window.result, previous)
        self.assertEqual(window.edited_study(), pending)
        self.assertEqual(window.viewport.display_positions.tolist(), positions.tolist())
        self.assertIn("cancelled", window.status.text())

    def test_window_close_defers_archive_reader_destruction(self):
        from vinkulum_studio.articulated_window import ArticulatedWindow
        from vinkulum_studio.examples3d import double_pendulum
        from vinkulum_studio.static_window import StaticWindow

        for window_type, arguments in (
            (StaticWindow, ()),
            (ArticulatedWindow, (double_pendulum(),)),
        ):
            with self.subTest(window=window_type.__name__):
                self.check_window_close(window_type, arguments)

    def check_window_close(self, window_type, arguments):
        window = window_type(*arguments)
        if hasattr(window, "interpreter"):
            window.interpreter.setText(sys.executable)
        self.assertTrue(window.run_button.isEnabled())
        window._discard_allowed = lambda: True
        self.addCleanup(lambda: window.close() if isValid(window) else None)
        entered, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)

        def load():
            entered.set()
            if not release.wait(5):
                raise RuntimeError("Test did not release archive reader")
            return object()

        window.archive.start(load, context=("result", None))
        self.wait_for(entered.is_set)
        self.assertFalse(window.run_button.isEnabled())
        before = time.monotonic()
        window.close()
        self.assertLess(time.monotonic() - before, 0.25)
        self.assertTrue(isValid(window))
        release.set()
        self.wait_for(lambda: not isValid(window))

    @unittest.skipUnless(importlib.util.find_spec("build123d"), "Requires OCCT 8")
    def test_editor_uses_the_body_and_document_prepared_off_gui(self):
        from vinkulum_studio.cad_dialog import CadDialog
        from vinkulum_studio.document import Body
        from vinkulum_studio.editor import EditorWindow

        window = EditorWindow()
        window._discard_allowed = lambda: True
        self.addCleanup(window.close)
        window.new_project()
        window.show()
        before = window.project
        original = Body.from_dict
        identities = []

        def decode(data):
            identities.append(threading.get_ident())
            if threading.get_ident() == threading.main_thread().ident:
                raise ValueError("CAD body was reconstructed on the GUI thread")
            return original(data)

        def create():
            dialog = self.app.activeModalWidget()
            if not isinstance(dialog, CadDialog):
                self.fail("CAD command did not open the expected dialog")
            dialog.name.setText("Background-admitted solid")
            dialog.start()

        timer = QTimer()
        timer.setSingleShot(True)
        timer.timeout.connect(create)
        timer.start(50)
        timed_out = []

        def timeout():
            dialog = self.app.activeModalWidget()
            timed_out.append(
                dialog.status.text()
                if isinstance(dialog, CadDialog)
                else "CAD dialog did not complete"
            )
            if dialog is not None:
                dialog.reject()

        watchdog = QTimer()
        watchdog.setSingleShot(True)
        watchdog.timeout.connect(timeout)
        watchdog.start(15000)
        try:
            with patch.object(Body, "from_dict", side_effect=decode):
                window.open_cad()
        finally:
            timer.stop()
            watchdog.stop()
        self.assertFalse(timed_out, timed_out)
        self.assertTrue(identities)
        self.assertNotIn(threading.get_ident(), identities)
        self.assertEqual(len(window.project.bodies), 1, window.status.text())
        accepted = window.project
        window.undo()
        self.assertEqual(
            window.project, replace(before, revision=window.project.revision)
        )
        self.assertGreater(window.project.revision, accepted.revision)
        window.redo()
        self.assertEqual(
            window.project, replace(accepted, revision=window.project.revision)
        )

    @unittest.skipUnless(importlib.util.find_spec("build123d"), "Requires OCCT 8")
    def test_cad_dialog_cancel_keeps_reader_and_files_alive_until_it_finishes(self):
        from vinkulum_studio.cad_dialog import CadDialog
        from vinkulum_studio.cad_history_dialog import CadHistoryDialog
        from vinkulum_studio.document import Body, Project, new_id

        project = Project(new_id(), "CAD cancellation")
        body = Body(new_id(), "Stock", dimensions=(0.12, 0.075, 0.02), mass=1.404)
        prepared = replace(project, bodies=(body,))
        for history in (False, True):
            with self.subTest(history=history):
                dialog = (
                    CadHistoryDialog(prepared, body.id)
                    if history
                    else CadDialog(project, None)
                )
                self.check_cad_cancel(dialog, history)

    def check_cad_cancel(self, dialog, history):
        from vinkulum_studio.cad_admission import admit_cad_response

        project = dialog.project
        self.addCleanup(dialog.reject)
        entered, release = threading.Event(), threading.Event()
        self.addCleanup(release.set)

        def admit(*args):
            entered.set()
            if not release.wait(5):
                raise RuntimeError("Test did not release CAD validator")
            return admit_cad_response(*args)

        dialog.show()
        with patch(
            "vinkulum_studio.cad_controller.admit_cad_response", side_effect=admit
        ):
            (dialog.start_preview if history else dialog.start)()
            self.wait_for(entered.is_set)
            temporary = Path(dialog.controller.temporary.name)
            before = time.monotonic()
            dialog.reject()
            self.assertLess(time.monotonic() - before, 0.25)
            self.assertTrue(dialog.isVisible())
            self.assertTrue(temporary.exists())
            self.assertIsNone(dialog.result_project)
            release.set()
            self.wait_for(lambda: not dialog.isVisible())
        self.assertIsNone(dialog.controller.process)
        self.assertIsNone(dialog.result_project)
        self.assertFalse(temporary.exists())
        self.assertEqual(dialog.result(), QDialog.DialogCode.Rejected)
        self.assertIs(dialog.project, project)


if __name__ == "__main__":
    unittest.main()
