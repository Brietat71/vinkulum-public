"""Cancellable CalculiX supervision with CPU admission and off-GUI validation."""

import os

from PySide6.QtCore import QObject, QProcess, QThread, QTimer, Signal, Slot

from .calculix import (
    StaticStudy,
    finish_run,
    identify_engine,
    parse_version,
    prepare_run,
)
from .cpu_scheduler import application_scheduler
from .engine_environment import external_engine_environment
from .execution import ExecutionPlan
from .model import finite_number, write_json


class StaticValidation(QThread):
    def __init__(self, run, code, parent):
        super().__init__(parent)
        self.calculation, self.code = run, code
        self.result, self.error = None, None

    def run(self):
        try:
            report = finish_run(self.calculation, self.code, publish=False)
            if self.isInterruptionRequested():
                return
            # Encode and write large reports here; the GUI only commits a rename.
            write_json(self.calculation.root / "validated-result.json", report)
            self.result = (self.calculation.study, report, self.calculation.root)
        except (OSError, ValueError, TypeError, RuntimeError) as error:
            self.error = str(error)


class StaticController(QObject):
    busy_changed = Signal(bool)
    stage_changed = Signal(str)
    completed = Signal(object)
    problem = Signal(str)

    def __init__(self, parent=None, *, scheduler=None):
        super().__init__(parent)
        self.scheduler = scheduler if scheduler is not None else application_scheduler()
        self._ticket = None
        self.process, self.validation, self.last_result, self.run = (
            None,
            None,
            None,
            None,
        )
        self._reason = None
        self._version_output = b""
        self._stage = None
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self._timed_out)

    @Slot()
    def _timed_out(self):
        self.cancel("CalculiX timed out.")

    def start(self, study, directory, *, executable="ccx", timeout=60, threads=None):
        if self.process is not None:
            raise RuntimeError("A CalculiX calculation is already running.")
        if not isinstance(study, StaticStudy):
            raise TypeError("A validated static study is required.")
        if not finite_number(timeout) or not 0 < timeout <= 300:
            raise ValueError("Timeout must be between 0 and 300 seconds.")
        self._plan = ExecutionPlan(
            min(4, self.scheduler.capacity) if threads is None else threads,
            self.scheduler.capacity,
        )
        self._engine, self._engine_hash = identify_engine(executable)
        self._study, self._directory, self._timeout = study, directory, timeout
        self.run, self._reason, self._version_output = None, None, b""
        self._stage = "queued"
        self.process = QProcess(self)  # Busy sentinel, including the queue wait.
        self._ticket = self.scheduler.request(self._plan.threads)
        self._ticket.ready.connect(self._allocated)
        self.busy_changed.emit(True)
        if self.process is not None:
            self.stage_changed.emit("Waiting for CalculiX CPU resources…")

    @Slot(object)
    def _allocated(self, plan):
        if self.process is None or self._reason:
            return
        self.process.deleteLater()
        self.process = None
        self._launch("version", ["-v"], 10)

    def _launch(self, stage, arguments, timeout):
        process = QProcess(self)
        self.process, self._stage = process, stage
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.setProcessEnvironment(
            external_engine_environment(threads=self._plan.threads, engine="calculix")
        )
        if stage == "solve":
            process.setWorkingDirectory(str(self.run.root))
            process.setStandardOutputFile(str(self.run.root / "solver.log"))
        else:
            process.readyReadStandardOutput.connect(self._output_ready)
        process.finished.connect(self._process_finished)
        process.errorOccurred.connect(self._process_error)
        self.timer.start(max(1, round(timeout * 1000)))
        self.stage_changed.emit(
            "Identifying CalculiX…"
            if stage == "version"
            else f"Solving linear statics · {self._plan.threads} CPU threads…"
        )
        if process is self.process and not self._reason:
            process.start(str(self._engine), arguments)

    @Slot()
    def _output_ready(self):
        self._drain(self.sender())

    def _drain(self, process):
        if process is not None and process is self.process and self._stage == "version":
            self._version_output += bytes(process.readAllStandardOutput())
            if len(self._version_output) > 32768:
                self._version_output = self._version_output[:32768]
                self.cancel("CalculiX version output exceeded its budget.")

    def cancel(self, reason="Calculation cancelled. Previous result preserved."):
        process = self.process
        if process is None:
            return
        self._reason = reason
        self.timer.stop()
        if self.validation is not None:
            self.validation.requestInterruption()
        elif process.state() == QProcess.ProcessState.NotRunning:
            self._settle(process, reason)
        else:
            process.kill()

    @Slot(QProcess.ProcessError)
    def _process_error(self, error):
        process = self.sender()
        if (
            process is not None
            and process is self.process
            and error == QProcess.ProcessError.FailedToStart
        ):
            self._settle(process, f"Could not start CalculiX: {process.errorString()}")

    @Slot(int, QProcess.ExitStatus)
    def _process_finished(self, code, status):
        process = self.sender()
        if process is None or process is not self.process:
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
                    execution=self._plan,
                )
                self.process = None
                process.deleteLater()
                self._launch("solve", ["-i", "study"], self._timeout)
                return
        except (OSError, ValueError, TypeError, RuntimeError) as error:
            self._settle(process, str(error))
            return
        worker = StaticValidation(self.run, code, self)
        self.validation = worker
        worker.finished.connect(self._reader_finished)
        self.stage_changed.emit("Checking forces, moments and energy…")
        if self.validation is not worker:
            return
        if self._reason:
            self._validated(worker)
        else:
            worker.start()

    @Slot()
    def _reader_finished(self):
        worker = self.sender()
        if worker is self.validation:
            self._validated(worker)

    def _validated(self, worker):
        if worker is None or worker is not self.validation:
            return
        self.validation = None
        worker.deleteLater()
        message = self._reason or worker.error
        if message or worker.result is None:
            self._settle(
                self.process, message or "No validated CalculiX result returned."
            )
            return
        try:
            os.replace(
                self.run.root / "validated-result.json", self.run.root / "result.json"
            )
        except OSError as error:
            self._settle(self.process, str(error))
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
        if self.run is not None:
            try:
                (self.run.root / "validated-result.json").unlink(missing_ok=True)
            except OSError as error:
                message = f"{message or 'Calculation settled.'}\nCould not remove staged result: {error}"
        self.busy_changed.emit(False)
        if message:
            if self.run is not None:
                message += f"\nInput and diagnostics: {self.run.root}"
            self.problem.emit(message)

    def shutdown(self):
        process = self.process
        if process is not None:
            self.cancel("Window closed. Calculation stopped.")
            if self.validation is not None:
                worker = self.validation
                if not worker.wait(2000):
                    raise RuntimeError(
                        "CalculiX result checks are stopping; keep their owner alive."
                    )
                self._validated(worker)
            elif self.process is process:
                process.waitForFinished(2000)
                if (
                    process is self.process
                    and process.state() == QProcess.ProcessState.NotRunning
                ):
                    self._settle(process, self._reason)
