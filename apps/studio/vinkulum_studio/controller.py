"""Qt process supervision. A subprocess is isolation, not a security sandbox."""

import sys
import tempfile
import uuid
import zipfile
from dataclasses import asdict
from pathlib import Path

from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, QTimer, Signal

from .document import MAX_PROJECT_BYTES, Project
from .mechanism import MechanicalResult
from .model import MAX_RESULT_BYTES, Result, read_json, write_json


class Controller(QObject):
    busy_changed = Signal(bool)
    completed = Signal(object)
    problem = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.process = None
        self.last_result = None
        self._temporary = None
        self._cancelled = False
        self._log = b""

    def start(self, parameters):
        if self.process is not None:
            raise RuntimeError("A calculation is already running.")
        if isinstance(parameters, Project) and parameters.diagnostics():
            raise ValueError("\n".join(d.message for d in parameters.diagnostics()))
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
                    "kind": "mechanism"
                    if isinstance(parameters, Project)
                    else "pendulum",
                },
            )
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
            "RAYON_NUM_THREADS": "2",
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "PYTHONUNBUFFERED": "1",
        }.items():
            environment.insert(key, value)
        process.setProcessEnvironment(environment)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.readyReadStandardOutput.connect(lambda: self._drain(process))
        process.finished.connect(
            lambda code, status: self._finish(process, code, status)
        )
        process.errorOccurred.connect(lambda error: self._error(process, error))
        self.busy_changed.emit(True)
        program, arguments = self._worker_command(directory)
        process.start(program, arguments)

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
        result = None
        message = None
        try:
            if self._cancelled:
                message = "Calculation stopped. No new result published."
            elif status == QProcess.ExitStatus.CrashExit:
                message = f"The calculation process exited abnormally (code {code})."
            else:
                data = read_json(
                    self._output, max(MAX_RESULT_BYTES, MAX_PROJECT_BYTES * 2)
                )
                if code != 0:
                    detail = str(data.get("message", "Error without diagnostic"))[:2000]
                    message = f"Calculation failed: {detail}"
                else:
                    result = (
                        MechanicalResult.read(
                            data, self._run_id, self._parameters, self._output.parent
                        )
                        if isinstance(self._parameters, Project)
                        else Result.from_dict(data, self._run_id, self._parameters)
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
            message = f"Worker output rejected: {exc}"
        self._cleanup(process)
        if result is not None:
            self.last_result = result
            self.completed.emit(result)
        else:
            self.problem.emit(message or "Calculation failed.")

    def _cleanup(self, process):
        self.process = None
        if self._temporary is not None:
            self._temporary.cleanup()
            self._temporary = None
        process.deleteLater()
        self.busy_changed.emit(False)

    def shutdown(self):
        if self.process is not None:
            self._cancelled = True
            self.process.kill()
            self.process.waitForFinished(1000)
