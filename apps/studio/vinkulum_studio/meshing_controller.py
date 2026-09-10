"""Direct Gmsh supervision; numerical mesh validation runs off the GUI thread."""

import os

from PySide6.QtCore import QObject, QProcess, QThread, QTimer, Signal, Slot

from .cpu_scheduler import application_scheduler
from .engine_environment import external_engine_environment
from .execution import ExecutionPlan
from .meshing import (
    MeshRequest,
    finish_mesh,
    identify_mesher,
    parse_info,
    prepare_mesh,
    save_mesh_result,
)
from .model import finite_number


class MeshValidation(QThread):
    def __init__(self, run, parent=None):
        super().__init__(parent)
        self.mesh_run = run
        self.result = None
        self.error = None

    def run(self):
        try:
            self.result = finish_mesh(
                self.mesh_run,
                0,
                publish=False,
                cancellation=self.isInterruptionRequested,
            )
            if not self.isInterruptionRequested():
                save_mesh_result(
                    self.mesh_run, self.result, filename="validated-result.json"
                )
        except (OSError, ValueError, TypeError, RuntimeError, KeyError) as error:
            self.error = str(error)


class MeshingController(QObject):
    busy_changed = Signal(bool)
    stage_changed = Signal(str)
    completed = Signal(object)
    failed = Signal(str, str)

    def __init__(self, parent=None, *, scheduler=None):
        super().__init__(parent)
        self.scheduler = scheduler if scheduler is not None else application_scheduler()
        self._ticket = None
        self.process = None
        self.validation = None
        self.busy = False
        self.run = None
        self.last_result = None
        self.failure_kind = None
        self._reason = None
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self._timed_out)

    @Slot()
    def _timed_out(self):
        self.cancel("timed_out")

    def start(self, request, directory, *, executable="gmsh", timeout=120):
        if self.busy:
            raise RuntimeError("A mesh operation is already running.")
        if not isinstance(request, MeshRequest):
            raise TypeError("A captured MeshRequest is required.")
        if not finite_number(timeout) or not 0 < timeout <= 600:
            raise ValueError("Meshing timeout must be between 0 and 600 seconds.")
        self._plan = ExecutionPlan(request.threads, self.scheduler.capacity)
        self._engine, self._identity = identify_mesher(executable)
        self._request, self._directory, self._timeout = request, directory, timeout
        self.run = None
        self._reason = None
        self.failure_kind = None
        self._output = b""
        self.busy = True
        self._ticket = self.scheduler.request(request.threads)
        self._ticket.ready.connect(self._allocated)
        self.busy_changed.emit(True)
        if self.busy:
            self.stage_changed.emit("Waiting for Gmsh CPU resources…")

    @Slot(object)
    def _allocated(self, plan):
        if not self.busy or self._reason:
            return
        self._launch("version", ["-info"], 10)

    def _launch(self, stage, arguments, timeout):
        if self._reason:
            self._settle(self._reason, "Meshing cancelled before execution.")
            return
        process = QProcess(self)
        self.process = process
        self._stage = stage
        process.setProcessEnvironment(
            external_engine_environment(threads=self._plan.threads, engine="gmsh")
        )
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        if stage == "mesh":
            process.setWorkingDirectory(str(self.run.root))
            process.setStandardOutputFile(str(self.run.root / "mesher.log"))
        else:
            process.readyReadStandardOutput.connect(self._output_ready)
        process.finished.connect(self._process_finished)
        process.errorOccurred.connect(self._process_error)
        self.timer.start(max(1, round(timeout * 1000)))
        self.stage_changed.emit(
            "Identifying Gmsh and OCCT…"
            if stage == "version"
            else f"Meshing the captured CAD solid · {self._plan.threads} CPU threads…"
        )
        # A signal handler can cancel synchronously while the process has not
        # started yet. Do not launch it after such a cancellation.
        if process is not self.process or not self.busy:
            return
        if self._reason:
            self._drop_process(process)
            self._settle(self._reason, "Meshing cancelled before execution.")
            return
        process.start(str(self._engine), arguments)

    @Slot()
    def _output_ready(self):
        self._drain(self.sender())

    @Slot(int, QProcess.ExitStatus)
    def _process_finished(self, code, status):
        self._finished(self.sender(), code, status)

    @Slot(QProcess.ProcessError)
    def _process_error(self, error):
        self._error(self.sender(), error)

    def _drain(self, process):
        if process is not None and process is self.process and self._stage == "version":
            self._output += bytes(process.readAllStandardOutput())
            if len(self._output) > 32768:
                self._output = self._output[:32768]
                self.cancel("invalid_version")

    def cancel(self, kind="cancelled"):
        if not self.busy:
            return
        self._reason = kind
        self.timer.stop()
        if self.process is not None:
            if self.process.state() == QProcess.ProcessState.NotRunning:
                self._drop_process(self.process)
                self._settle(kind, "Meshing cancelled before execution.")
            else:
                self.process.kill()
        if self.validation is not None:
            self.validation.requestInterruption()
        elif self.process is None and self.busy:
            self._settle(kind, "Meshing cancelled before execution.")

    def _error(self, process, error):
        if (
            process is not None
            and process is self.process
            and error == QProcess.ProcessError.FailedToStart
        ):
            self._drop_process(process)
            self._settle(
                "start_failed", "Could not start Gmsh: " + process.errorString()
            )

    def _drop_process(self, process):
        if self.process is process:
            self.process = None
        process.deleteLater()

    def _finished(self, process, code, status):
        if process is None or process is not self.process:
            return
        self._drain(process)
        self._drop_process(process)
        if self._reason:
            self._settle(
                self._reason, "Meshing stopped. The previous mesh is preserved."
            )
            return
        if status == QProcess.ExitStatus.CrashExit:
            self._settle("crashed", "Gmsh crashed. The previous mesh is preserved.")
            return
        if code != 0:
            detail = (
                self._output.decode(errors="replace")
                if self._stage == "version"
                else f"Inspect {self.run.root / 'mesher.log'}"
            )
            self._settle("mesher_failed", f"Gmsh exited with code {code}.\n{detail}")
            return
        if self._stage == "version":
            info = self._output.decode(errors="replace")
            try:
                parse_info(info, self._request)
            except ValueError as error:
                self._settle("invalid_version", str(error))
                return
            try:
                self.run = prepare_mesh(
                    self._request,
                    self._directory,
                    self._engine,
                    self._identity,
                    info,
                )
            except (OSError, ValueError, TypeError) as error:
                self._settle("capture_failed", str(error))
                return
            self._launch(
                "mesh",
                [
                    "-nt",
                    str(self._plan.threads),
                    "mesh.geo",
                    "-3",
                    "-format",
                    "msh2",
                    "-o",
                    "mesh.msh",
                ],
                self._timeout,
            )
            return
        self.stage_changed.emit("Checking element Jacobians and boundary groups…")
        if self._reason:
            self._settle(self._reason, "Mesh validation cancelled before execution.")
            return
        worker = MeshValidation(self.run, self)
        self.validation = worker
        worker.finished.connect(self._reader_finished)
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
                self._reason, "Meshing stopped. The previous mesh is preserved."
            )
            return
        if worker.error or worker.result is None:
            self._settle("invalid_mesh", worker.error or "No validated mesh returned.")
            return
        try:
            os.replace(
                self.run.root / "validated-result.json", self.run.root / "result.json"
            )
        except (OSError, ValueError, TypeError) as error:
            self._settle("write_failed", str(error))
            return
        self.last_result = (worker.result, self.run.root)
        self._settle()
        self.completed.emit(self.last_result)

    def _settle(self, kind=None, message=""):
        if not self.busy:
            return
        self.timer.stop()
        self.busy = False
        self.failure_kind = kind
        if self._ticket is not None:
            self._ticket.release()
            self._ticket = None
        if self.run is not None:
            try:
                (self.run.root / "validated-result.json").unlink(missing_ok=True)
            except OSError as error:
                kind = kind or "cleanup_failed"
                self.failure_kind = kind
                message += f"\nCould not remove staged result: {error}"
        self.busy_changed.emit(False)
        if kind:
            self.failed.emit(kind, message)

    def shutdown(self):
        if self.busy:
            self.cancel()
        if self.process is not None:
            self.process.waitForFinished(2000)
        if self.validation is not None and not self.validation.wait(2000):
            raise RuntimeError(
                "Mesh validation is still stopping; keep its owner alive."
            )
        if self.validation is not None:
            self._validated(self.validation)
