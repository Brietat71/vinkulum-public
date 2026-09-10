"""Direct, cancellable Qt supervision of the optional CalculiX executable."""

from PySide6.QtCore import QObject, QProcess, QTimer, Signal

from .calculix import (
    StaticStudy,
    finish_run,
    identify_engine,
    parse_version,
    prepare_run,
)
from .engine_environment import external_engine_environment
from .model import finite_number


class StaticController(QObject):
    busy_changed = Signal(bool)
    stage_changed = Signal(str)
    completed = Signal(object)
    problem = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.process = None
        self.last_result = None
        self.run = None
        self._reason = None
        self._version_output = b""
        self._stage = None
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(lambda: self.cancel("CalculiX timed out."))

    def start(self, study, directory, *, executable="ccx", timeout=60):
        if self.process is not None:
            raise RuntimeError("A CalculiX calculation is already running.")
        if not isinstance(study, StaticStudy):
            raise TypeError("A validated static study is required.")
        if not finite_number(timeout) or not 0 < timeout <= 300:
            raise ValueError("Timeout must be between 0 and 300 seconds.")
        self._engine, self._engine_hash = identify_engine(executable)
        self._study, self._directory, self._timeout = study, directory, timeout
        self.run = None
        self._reason = None
        self._version_output = b""
        self._launch("version", ["-v"], 10)

    def _launch(self, stage, arguments, timeout):
        process = QProcess(self)
        self.process, self._stage = process, stage
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        environment = external_engine_environment()
        for name in (
            "OMP_NUM_THREADS",
            "CCX_NPROC_RESULTS",
            "CCX_NPROC_EQUATION_SOLVER",
        ):
            environment.insert(name, "1")
        process.setProcessEnvironment(environment)
        if stage == "solve":
            process.setWorkingDirectory(str(self.run.root))
            process.setStandardOutputFile(str(self.run.root / "solver.log"))
        else:
            process.readyReadStandardOutput.connect(lambda: self._drain(process))
        process.finished.connect(
            lambda code, status: self._finished(process, code, status)
        )
        process.errorOccurred.connect(lambda error: self._error(process, error))
        self.timer.start(max(1, round(timeout * 1000)))
        if stage == "version":
            self.busy_changed.emit(True)
        self.stage_changed.emit(
            "Identifying CalculiX…" if stage == "version" else "Solving linear statics…"
        )
        process.start(str(self._engine), arguments)

    def _drain(self, process):
        if process is self.process and self._stage == "version":
            self._version_output += bytes(process.readAllStandardOutput())
            if len(self._version_output) > 32768:
                self._version_output = self._version_output[:32768]
                self.cancel("CalculiX version output exceeded its budget.")

    def cancel(self, reason="Calculation cancelled. Previous result preserved."):
        if self.process is not None:
            self._reason = reason
            self.timer.stop()
            # This is ccx itself, not a Python wrapper which could leave ccx alive.
            self.process.kill()

    def _error(self, process, error):
        if process is self.process and error == QProcess.ProcessError.FailedToStart:
            self._settle(process, f"Could not start CalculiX: {process.errorString()}")

    def _finished(self, process, code, status):
        if process is not self.process:
            return
        self.timer.stop()
        self._drain(process)
        try:
            if self._reason:
                raise ValueError(self._reason)
            if status != QProcess.ExitStatus.NormalExit:
                raise ValueError(f"CalculiX exited abnormally (code {code}).")
            if self._stage == "version":
                # CalculiX 2.21 returns 201 for a successful -v probe.
                version = parse_version(
                    self._version_output.decode("utf-8", errors="replace")
                )
                self.run = prepare_run(
                    self._study,
                    self._directory,
                    self._engine,
                    self._engine_hash,
                    version,
                )
                self.process = None
                process.deleteLater()
                self._launch("solve", ["-i", "study"], self._timeout)
                return
            self.stage_changed.emit("Checking forces, moments and energy…")
            report = finish_run(self.run, code)
            result = (self.run.study, report, self.run.root)
        except (OSError, ValueError, TypeError, RuntimeError) as error:
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
            if self.run is not None:
                message += f"\nInput and diagnostics: {self.run.root}"
            self.problem.emit(message)

    def shutdown(self):
        process = self.process
        if process is not None:
            self.cancel("Window closed. Calculation stopped.")
            process.waitForFinished(2000)
            if (
                process is self.process
                and process.state() == QProcess.ProcessState.NotRunning
            ):
                self._settle(process, self._reason)
