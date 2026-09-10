"""Supervise the optional Pinocchio Python worker without importing its engine."""

import os
from pathlib import Path

from PySide6.QtCore import QObject, QProcess, QTimer, Signal

from .articulated import state_vector, tree_links
from .articulated_result import load_operators
from .document import save_project
from .engine_environment import external_engine_environment
from .model import finite_number, write_json


class PinocchioController(QObject):
    busy_changed = Signal(bool)
    completed = Signal(object)
    problem = Signal(str)
    stage_changed = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.process = None
        self.last_result = None
        self._reason = None
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(lambda: self.cancel("Pinocchio worker timed out."))

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
        environment = external_engine_environment()
        for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
            environment.insert(key, "1")
        environment.insert("PYTHONUNBUFFERED", "1")
        process.setProcessEnvironment(environment)
        process.setWorkingDirectory(str(root))
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.setStandardOutputFile(str(root / "worker.log"))
        process.finished.connect(
            lambda code, status: self._finished(process, code, status)
        )
        process.errorOccurred.connect(lambda error: self._error(process, error))
        self.timer.start(max(1, round(timeout * 1000)))
        self.busy_changed.emit(True)
        self.stage_changed.emit("Evaluating captured articulated operators…")
        process.start(
            str(python),
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
            self.process.kill()

    def _error(self, process, error):
        if process is self.process and error == QProcess.ProcessError.FailedToStart:
            self._settle(
                process, f"Could not start Pinocchio worker: {process.errorString()}"
            )

    def _finished(self, process, code, status):
        if process is not self.process:
            return
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
            self.stage_changed.emit("Checking captured inputs and operator identities…")
            result = load_operators(
                self._directory / "operators",
                expected_project=self._project,
                expected_state=self._state,
            )
        except (OSError, ValueError, TypeError, KeyError) as error:
            self._settle(process, str(error))
            return
        self.last_result = result
        self._settle(process)
        self.completed.emit(result)

    def _settle(self, process, message=None):
        if process is not self.process:
            return
        self.timer.stop()
        self.process = None
        process.deleteLater()
        self.busy_changed.emit(False)
        if message:
            self.problem.emit(f"{message}\nInputs and diagnostics: {self._directory}")

    def shutdown(self):
        process = self.process
        if process is not None:
            self.cancel("Window closed. Pinocchio worker stopped.")
            process.waitForFinished(2000)
            if (
                process is self.process
                and process.state() == QProcess.ProcessState.NotRunning
            ):
                self._settle(process, self._reason)
