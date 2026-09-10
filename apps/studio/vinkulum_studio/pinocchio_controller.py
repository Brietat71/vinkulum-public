"""Supervise the optional Pinocchio Python worker without importing its engine."""

import os
from pathlib import Path

from PySide6.QtCore import QObject, QProcess, QThread, QTimer, Signal, Slot

from .articulated import state_vector, tree_links
from .articulated_result import load_operators
from .cpu_scheduler import application_scheduler
from .document import save_project
from .engine_environment import external_engine_environment
from .model import finite_number, write_json


class OperatorValidation(QThread):
    def __init__(self, directory, project, state, parent):
        super().__init__(parent)
        self.directory, self.project, self.state = directory, project, state
        self.result, self.error = None, None

    def run(self):
        try:
            result = load_operators(
                self.directory / "operators",
                expected_project=self.project,
                expected_state=self.state,
            )
            if not self.isInterruptionRequested():
                self.result = result
        except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
            self.error = str(error)


class PinocchioController(QObject):
    busy_changed = Signal(bool)
    completed = Signal(object)
    problem = Signal(str)
    stage_changed = Signal(str)

    def __init__(self, parent=None, *, scheduler=None):
        super().__init__(parent)
        self.scheduler = scheduler if scheduler is not None else application_scheduler()
        self._ticket = None
        self.validation = None
        self.process = None
        self.last_result = None
        self._reason = None
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self._timed_out)

    @Slot()
    def _timed_out(self):
        self.cancel("Pinocchio worker timed out.")

    def start(self, project, state, directory, *, interpreter, timeout=30):
        if self.process is not None:
            raise RuntimeError("A Pinocchio analysis is already running.")
        links = tree_links(project)
        if not isinstance(state, dict) or not set(state) <= {
            "q",
            "velocity",
            "acceleration",
            "effort",
            "time_s",
        }:
            raise ValueError("Invalid articulated state.")
        captured = {
            key: state_vector(state.get(key), links, key, initial=key == "q").tolist()
            for key in ("q", "velocity", "acceleration", "effort")
        }
        captured["time_s"] = state.get("time_s", 0.0)
        if (
            not finite_number(captured["time_s"])
            or not 0 <= captured["time_s"] <= project.duration
        ):
            raise ValueError("Load time must lie within the project duration.")
        if not finite_number(timeout) or not 0 < timeout <= 120:
            raise ValueError("Timeout must be between 0 and 120 seconds.")
        # Resolving a venv's Python symlink would lose that environment.
        python = Path(interpreter).expanduser().absolute()
        if not python.is_file() or not os.access(python, os.X_OK):
            raise ValueError(
                "Select the Python executable in the separate Pinocchio environment."
            )
        root = Path(directory).resolve()
        root.mkdir(parents=True, exist_ok=False)
        save_project(root / "request-project.json", project)
        write_json(root / "request-state.json", captured)
        self._project, self._state, self._directory = project, captured, root
        self._reason = None
        process = QProcess(self)
        self.process = process
        # Single-state recursive algorithms are sequential. Reserve one CPU;
        # more OpenMP threads would not parallelise this operator evaluation.
        environment = external_engine_environment(threads=1, engine="pinocchio")
        process.setProcessEnvironment(environment)
        process.setWorkingDirectory(str(root))
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.setStandardOutputFile(str(root / "worker.log"))
        process.finished.connect(self._process_finished)
        process.errorOccurred.connect(self._process_error)
        self._python, self._timeout = python, timeout
        self._ticket = self.scheduler.request(1)
        self._ticket.ready.connect(self._allocated)
        self.busy_changed.emit(True)
        if self.process is not None:
            self.stage_changed.emit("Waiting for Pinocchio CPU resources…")

    @Slot(object)
    def _allocated(self, plan):
        process = self.process
        if process is None or self._reason:
            return
        self.timer.start(max(1, round(self._timeout * 1000)))
        self.stage_changed.emit("Evaluating captured articulated operators…")
        if process is not self.process or self._reason:
            return
        root = self._directory
        process.start(
            str(self._python),
            [
                "-m",
                "vinkulum_studio.pinocchio_backend",
                str(root / "request-project.json"),
                "--state",
                str(root / "request-state.json"),
                "--output",
                str(root / "operators"),
            ],
        )

    def cancel(self, reason="Pinocchio analysis cancelled. Previous result preserved."):
        if self.process is not None:
            self._reason = reason
            self.timer.stop()
            if self.validation is not None:
                self.validation.requestInterruption()
            elif self.process.state() == QProcess.ProcessState.NotRunning:
                self._settle(self.process, reason)
            else:
                self.process.kill()

    @Slot(QProcess.ProcessError)
    def _process_error(self, error):
        self._error(self.sender(), error)

    @Slot(int, QProcess.ExitStatus)
    def _process_finished(self, code, status):
        self._finished(self.sender(), code, status)

    def _error(self, process, error):
        if (
            process is not None
            and process is self.process
            and error == QProcess.ProcessError.FailedToStart
        ):
            self._settle(
                process, f"Could not start Pinocchio worker: {process.errorString()}"
            )

    def _finished(self, process, code, status):
        if process is None or process is not self.process:
            return
        self.timer.stop()
        try:
            if self._reason:
                raise ValueError(self._reason)
            if status != QProcess.ExitStatus.NormalExit or code != 0:
                with (self._directory / "worker.log").open("rb") as stream:
                    stream.seek(0, 2)
                    size = stream.tell()
                    stream.seek(max(0, size - 4096))
                    detail = stream.read(4096).decode("utf-8", errors="replace")
                raise ValueError(f"Pinocchio worker failed (code {code}).\n{detail}")
        except (OSError, ValueError, TypeError, KeyError) as error:
            self._settle(process, str(error))
            return
        worker = OperatorValidation(self._directory, self._project, self._state, self)
        self.validation = worker
        worker.finished.connect(self._reader_finished)
        self.stage_changed.emit("Checking captured inputs and operator identities…")
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
        error = self._reason or worker.error
        if error or worker.result is None:
            self._settle(self.process, error or "No validated operators returned.")
            return
        self.last_result = worker.result
        self._settle(self.process)
        self.completed.emit(self.last_result)

    def _settle(self, process, message=None):
        if process is None or process is not self.process:
            return
        self.timer.stop()
        self.process = None
        process.deleteLater()
        if self._ticket is not None:
            self._ticket.release()
            self._ticket = None
        self.busy_changed.emit(False)
        if message:
            self.problem.emit(f"{message}\nInputs and diagnostics: {self._directory}")

    def shutdown(self):
        process = self.process
        if process is not None:
            self.cancel("Window closed. Pinocchio worker stopped.")
            if self.validation is not None:
                worker = self.validation
                if not worker.wait(2000):
                    raise RuntimeError(
                        "Operator validation is still stopping; keep its owner alive."
                    )
                self._validated(worker)
                return
            process.waitForFinished(2000)
            if (
                process is self.process
                and process.state() == QProcess.ProcessState.NotRunning
            ):
                self._settle(process, self._reason)
