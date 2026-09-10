"""Qt process supervision. A subprocess is isolation, not a security sandbox."""

import sys
import tempfile
import uuid
import zipfile
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

from .cpu_scheduler import application_scheduler
from .document import MAX_PROJECT_BYTES, Project
from .execution import ExecutionPlan
from .mechanism import MechanicalResult
from .model import MAX_RESULT_BYTES, Result, read_json, write_json


class ResultValidation(QThread):
    """Own no widgets; validate captured inputs and arrays off the GUI thread."""

    def __init__(self, output, run_id, parameters, plan, code, parent):
        super().__init__(parent)
        self.output, self.run_id, self.parameters = output, run_id, parameters
        self.plan, self.code = plan, code
        self.result, self.message = None, None

    def run(self):
        try:
            data = read_json(self.output, max(MAX_RESULT_BYTES, MAX_PROJECT_BYTES * 2))
            if self.code != 0:
                self.message = f"Calculation failed: {str(data.get('message', 'Error without diagnostic'))[:2000]}"
                return
            self.plan.check_manifest(data.get("manifest", {}).get("execution"))
            self.result = (
                MechanicalResult.read(
                    data, self.run_id, self.parameters, self.output.parent
                )
                if isinstance(self.parameters, Project)
                else Result.from_dict(data, self.run_id, self.parameters)
            )
        except (
            OSError,
            ValueError,
            TypeError,
            KeyError,
            AttributeError,
            EOFError,
            zipfile.BadZipFile,
        ) as exc:
            self.message = f"Worker output rejected: {exc}"


class Controller(QObject):
    busy_changed = Signal(bool)
    completed = Signal(object)
    problem = Signal(str)
    stage_changed = Signal(str)

    def __init__(self, parent=None, *, scheduler=None):
        super().__init__(parent)
        self.process = None
        self.last_result = None
        self._temporary = None
        self._cancelled = False
        self._log = b""
        self.scheduler = scheduler if scheduler is not None else application_scheduler()
        self._ticket = None
        self._validation = None
        self.threads = self.scheduler.capacity

    def start(self, parameters, *, threads=None):
        if self.process is not None:
            raise RuntimeError("A calculation is already running.")
        if isinstance(parameters, Project) and parameters.diagnostics():
            raise ValueError("\n".join(d.message for d in parameters.diagnostics()))
        self._plan = ExecutionPlan(
            self.threads if threads is None else threads, self.scheduler.capacity
        )
        self._temporary = tempfile.TemporaryDirectory(prefix="vinkulum-run-")
        directory = Path(self._temporary.name)
        self._run_id = str(uuid.uuid4())
        self._parameters = parameters
        self._output = directory / "result.json"
        try:
            write_json(
                directory / "input.json",
                {
                    "run_id": self._run_id,
                    "parameters": asdict(parameters),
                    "execution": asdict(self._plan),
                    "kind": "mechanism"
                    if isinstance(parameters, Project)
                    else "pendulum",
                },
            )
            self._command = self._worker_command(directory)
        except Exception:
            self._temporary.cleanup()
            self._temporary = None
            raise
        process = QProcess(self)
        self.process = process
        self._cancelled = False
        self._log = b""
        environment = QProcessEnvironment.systemEnvironment()
        for key, value in {
            "RAYON_NUM_THREADS": str(self._plan.threads),
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "PYTHONUNBUFFERED": "1",
        }.items():
            environment.insert(key, value)
        process.setProcessEnvironment(environment)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.readyReadStandardOutput.connect(self._process_output)
        process.finished.connect(self._process_finished)
        process.errorOccurred.connect(self._process_error)
        self._ticket = self.scheduler.request(self._plan.threads)
        self._ticket.ready.connect(self._allocated)
        self.busy_changed.emit(True)
        if self.process is process:
            self.stage_changed.emit("Waiting for native CPU resources…")

    @Slot(object)
    def _allocated(self, plan):
        process = self.process
        if process is None or self._cancelled:
            return
        program, arguments = self._command
        self.stage_changed.emit(f"Running with {plan.threads} native CPU threads…")
        if process is self.process and not self._cancelled:
            process.start(program, arguments)

    @Slot()
    def _process_output(self):
        self._drain(self.sender())

    @Slot(int, QProcess.ExitStatus)
    def _process_finished(self, code, status):
        self._finish(self.sender(), code, status)

    @Slot(QProcess.ProcessError)
    def _process_error(self, error):
        self._error(self.sender(), error)

    def _worker_command(self, directory):
        if getattr(sys, "frozen", False):
            return sys.executable, [
                "--worker",
                str(directory / "input.json"),
                str(self._output),
            ]
        return sys.executable, [
            "-m",
            "vinkulum_studio.worker",
            str(directory / "input.json"),
            str(self._output),
        ]

    def _drain(self, process):
        if process is self.process:
            self._log = (self._log + bytes(process.readAllStandardOutput()))[-32_768:]

    def cancel(self):
        process = self.process
        if process is None:
            return
        self._cancelled = True
        if self._validation is not None:
            return  # Keep files alive until the reader finishes; discard its output.
        if process.state() == QProcess.ProcessState.NotRunning:
            self._cleanup(process)
            self.problem.emit("Calculation stopped. No new result published.")
            return
        process.terminate()
        QTimer.singleShot(
            1000, lambda: process.kill() if self.process is process else None
        )

    def _error(self, process, error):
        if process is self.process and error == QProcess.ProcessError.FailedToStart:
            message = process.errorString()
            self._cleanup(process)
            self.problem.emit(f"Could not start worker: {message}")

    def _finish(self, process, code, status):
        if process is not self.process:
            return
        self._drain(process)
        if self._cancelled or status == QProcess.ExitStatus.CrashExit:
            message = (
                "Calculation stopped. No new result published."
                if self._cancelled
                else f"The calculation process exited abnormally (code {code})."
            )
            self._cleanup(process)
            self.problem.emit(message)
            return
        validation = ResultValidation(
            self._output, self._run_id, self._parameters, self._plan, code, self
        )
        self._validation = validation
        validation.finished.connect(self._reader_finished)
        # Retain the CPU lease through validation (which also consumes CPU).
        self.stage_changed.emit("Checking captured inputs and trajectory…")
        if self._validation is not validation:
            return  # A synchronous close handler already settled this reader.
        if self._cancelled:
            self._validated(process, validation)
            return
        validation.start()

    @Slot()
    def _reader_finished(self):
        validation = self.sender()
        if validation is not None and validation is self._validation:
            self._validated(self.process, validation)

    def _validated(self, process, validation):
        if process is not self.process or validation is not self._validation:
            return
        result, message = validation.result, validation.message
        if self._cancelled:
            result, message = None, "Calculation stopped. No new result published."
        self._validation = None
        validation.deleteLater()
        self._cleanup(process)
        if result is not None:
            self.last_result = result
            self.completed.emit(result)
        else:
            self.problem.emit(message or "Calculation failed.")

    def _cleanup(self, process):
        self.process = None
        if self._ticket is not None:
            self._ticket.release()
            self._ticket = None
        if self._temporary is not None:
            self._temporary.cleanup()
            self._temporary = None
        process.deleteLater()
        self.busy_changed.emit(False)

    def shutdown(self):
        process = self.process
        if process is not None:
            self.cancel()
            if self._validation is not None:
                validation = self._validation
                validation.wait()
                self._validated(process, validation)
            elif self.process is process:
                process.kill()
                process.waitForFinished(1000)
