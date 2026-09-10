"""A supervised CAD transaction, with distinct failure causes and captured input."""

import hashlib
import json
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path

from PySide6.QtCore import (
    QObject,
    QProcess,
    QProcessEnvironment,
    QThread,
    QTimer,
    Signal,
    Slot,
)

from .cad_admission import CadAdmissionError, admit_cad_response
from .cpu_scheduler import application_scheduler
from .document import Project
from .engine_threads import thread_environment
from .execution import ExecutionPlan
from .model import finite_number, write_json


class CadValidation(QThread):
    def __init__(self, path, request, input_sha256, plan, code, project, parent):
        super().__init__(parent)
        self.arguments = (path, request, input_sha256, plan, code, project)
        self.result, self.error, self.kind = None, None, None

    def run(self):
        try:
            result = admit_cad_response(*self.arguments)
            if not self.isInterruptionRequested():
                self.result = result
        except CadAdmissionError as error:
            self.error, self.kind = str(error), error.kind
        except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
            self.error, self.kind = str(error), "invalid_result"


class CadController(QObject):
    completed = Signal(object)
    failed = Signal(str, str)
    busy_changed = Signal(bool)
    stage_changed = Signal(str)

    def __init__(self, parent=None, *, scheduler=None, project=None):
        super().__init__(parent)
        if project is not None and not isinstance(project, Project):
            raise TypeError("CAD admission requires an immutable project snapshot.")
        self.scheduler = scheduler if scheduler is not None else application_scheduler()
        self._ticket = None
        self.process = None
        self.temporary = None
        self.last_result = None
        self.last_admission = None
        self.validation = None
        self._context_project = project
        self.last_request = None
        self.last_log = ""
        self.failure_kind = None
        self._reason = None
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self._timed_out)

    @Slot()
    def _timed_out(self):
        self.cancel("timed_out")

    def start(
        self, request, *, timeout=60, executable=None, arguments=None, threads=None
    ):
        if self.process is not None:
            raise RuntimeError("A CAD operation is already running.")
        if not finite_number(timeout) or not 0 < timeout <= 120:
            raise ValueError("CAD timeout must be between 0 and 120 seconds.")
        captured = json.loads(json.dumps(request, allow_nan=False))
        if not isinstance(captured, dict):
            raise TypeError("A CAD operation request is required.")
        self._plan = ExecutionPlan(
            min(4, self.scheduler.capacity) if threads is None else threads,
            self.scheduler.capacity,
        )
        captured["execution"] = asdict(self._plan)
        self.temporary = tempfile.TemporaryDirectory(prefix="vinkulum-cad-")
        root = Path(self.temporary.name)
        try:
            write_json(root / "input.json", captured)
        except OSError, ValueError, TypeError:
            self.temporary.cleanup()
            self.temporary = None
            raise
        self._input_hash = hashlib.sha256(
            (root / "input.json").read_bytes()
        ).hexdigest()
        self.last_request = captured
        self.last_log, self.failure_kind, self._reason = "", None, None
        process = QProcess(self)
        self.process = process
        environment = QProcessEnvironment.systemEnvironment()
        for key, value in thread_environment(self._plan.threads, engine="occt").items():
            environment.insert(key, value)
        process.setProcessEnvironment(environment)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.setStandardOutputFile(str(root / "worker.log"))
        process.finished.connect(self._process_finished)
        process.errorOccurred.connect(self._process_error)
        if arguments is None:
            arguments = (
                ["--cad-worker"]
                if getattr(sys, "frozen", False)
                else ["-m", "vinkulum_studio.cad_worker"]
            )
        self._command = (
            str(executable or sys.executable),
            [*arguments, str(root / "input.json"), str(root / "output.json")],
        )
        self._timeout = timeout
        self._ticket = self.scheduler.request(self._plan.threads)
        self._ticket.ready.connect(self._allocated)
        self.busy_changed.emit(True)
        if self.process is process:
            self.stage_changed.emit("Waiting for CAD CPU resources…")

    @Slot(object)
    def _allocated(self, plan):
        process = self.process
        if process is None or self._reason:
            return
        self.timer.start(max(1, round(self._timeout * 1000)))
        self.stage_changed.emit(
            f"Building CAD geometry · {plan.threads} native CPU threads…"
        )
        if process is self.process and not self._reason:
            process.start(*self._command)

    @Slot(int, QProcess.ExitStatus)
    def _process_finished(self, code, status):
        self._finished(self.sender(), code, status)

    @Slot(QProcess.ProcessError)
    def _process_error(self, error):
        self._error(self.sender(), error)

    def cancel(self, reason="cancelled"):
        if self.process is not None:
            self._reason = reason
            self.timer.stop()
            if self.validation is not None:
                self.validation.requestInterruption()
            elif self.process.state() == QProcess.ProcessState.NotRunning:
                self._settle(self.process, reason, "CAD cancelled before execution.")
            else:
                self.process.kill()

    def _error(self, process, error):
        if (
            process is not None
            and process is self.process
            and error == QProcess.ProcessError.FailedToStart
        ):
            self._settle(
                process,
                "start_failed",
                "The CAD worker could not start: " + process.errorString(),
            )

    def _read_log(self):
        path = Path(self.temporary.name) / "worker.log"
        if path.is_file():
            with path.open("rb") as stream:
                stream.seek(0, 2)
                stream.seek(max(0, stream.tell() - 8192))
                self.last_log = stream.read(8192).decode("utf-8", errors="replace")

    def _finished(self, process, code, status):
        if process is None or process is not self.process:
            return
        self.timer.stop()
        try:
            self._read_log()
        except OSError as error:
            self._settle(process, "invalid_result", str(error))
            return
        if self._reason:
            self._settle(
                process, self._reason, "CAD stopped. The previous model is preserved."
            )
            return
        if status == QProcess.ExitStatus.CrashExit:
            self._settle(
                process,
                "crashed",
                "The CAD worker crashed. The previous model is preserved.",
            )
            return
        worker = CadValidation(
            Path(self.temporary.name) / "output.json",
            self.last_request,
            self._input_hash,
            self._plan,
            code,
            self._context_project,
            self,
        )
        self.validation = worker
        worker.finished.connect(self._reader_finished)
        self.stage_changed.emit(
            "Checking CAD geometry, mass properties and document attachments…"
        )
        if self.validation is not worker:
            return
        if self._reason:
            self._validated(worker)
        else:
            worker.start()

    @Slot()
    def _reader_finished(self):
        self._validated(self.sender())

    def _validated(self, worker):
        if worker is None or worker is not self.validation:
            return
        self.validation = None
        worker.deleteLater()
        if self._reason:
            self._settle(
                self.process,
                self._reason,
                "CAD validation cancelled. The previous model is preserved.",
            )
            return
        if worker.result is None:
            self._settle(
                self.process,
                worker.kind or "invalid_result",
                worker.error or "No validated CAD response returned.",
            )
            return
        admission = worker.result
        self.last_admission = admission
        self.last_result = admission.payload
        self._settle(self.process)
        self.completed.emit(admission.payload)

    def _settle(self, process, kind=None, message=""):
        if process is None or process is not self.process:
            return
        self.timer.stop()
        self.process = None
        process.deleteLater()
        if self._ticket is not None:
            self._ticket.release()
            self._ticket = None
        self.temporary.cleanup()
        self.temporary = None
        self.failure_kind = kind
        self.busy_changed.emit(False)
        if kind:
            self.failed.emit(kind, message)

    def shutdown(self):
        process = self.process
        if process is not None:
            self.cancel()
            if self.validation is not None:
                worker = self.validation
                if not worker.wait(2000):
                    raise RuntimeError(
                        "CAD validation is still stopping; keep its owner alive."
                    )
                self._validated(worker)
                return
            process.waitForFinished(2000)
            if (
                process is self.process
                and process.state() == QProcess.ProcessState.NotRunning
            ):
                self._settle(
                    process, "cancelled", "CAD stopped when its window closed."
                )
