"""Direct Gmsh supervision; numerical mesh validation runs off the GUI thread."""

from PySide6.QtCore import QObject, QProcess, QThread, QTimer, Signal

from .engine_environment import external_engine_environment
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
        except (OSError, ValueError, TypeError, RuntimeError, KeyError) as error:
            self.error = str(error)


class MeshingController(QObject):
    busy_changed = Signal(bool)
    stage_changed = Signal(str)
    completed = Signal(object)
    failed = Signal(str, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.process = None
        self.validation = None
        self.busy = False
        self.run = None
        self.last_result = None
        self.failure_kind = None
        self._reason = None
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(lambda: self.cancel("timed_out"))

    def start(self, request, directory, *, executable="gmsh", timeout=120):
        if self.busy:
            raise RuntimeError("A mesh operation is already running.")
        if not isinstance(request, MeshRequest):
            raise TypeError("A captured MeshRequest is required.")
        if not finite_number(timeout) or not 0 < timeout <= 600:
            raise ValueError("Meshing timeout must be between 0 and 600 seconds.")
        self._engine, self._identity = identify_mesher(executable)
        self._request, self._directory, self._timeout = request, directory, timeout
        self.run = None
        self._reason = None
        self.failure_kind = None
        self._output = b""
        self.busy = True
        self.busy_changed.emit(True)
        self._launch("version", ["-info"], 10)

    def _launch(self, stage, arguments, timeout):
        if self._reason:
            self._settle(self._reason, "Meshing cancelled before execution.")
            return
        process = QProcess(self)
        self.process = process
        self._stage = stage
        process.setProcessEnvironment(external_engine_environment())
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        if stage == "mesh":
            process.setWorkingDirectory(str(self.run.root))
            process.setStandardOutputFile(str(self.run.root / "mesher.log"))
        else:
            process.readyReadStandardOutput.connect(lambda: self._drain(process))
        process.finished.connect(
            lambda code, status: self._finished(process, code, status)
        )
        process.errorOccurred.connect(lambda error: self._error(process, error))
        self.timer.start(max(1, round(timeout * 1000)))
        self.stage_changed.emit(
            "Identifying Gmsh and OCCT…"
            if stage == "version"
            else "Meshing the captured CAD solid…"
        )
        # A signal handler can cancel synchronously while the process has not
        # started yet. Do not launch it after such a cancellation.
        if self._reason:
            self._drop_process(process)
            self._settle(self._reason, "Meshing cancelled before execution.")
            return
        process.start(str(self._engine), arguments)

    def _drain(self, process):
        if process is self.process and self._stage == "version":
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
            self.process.kill()
        if self.validation is not None:
            self.validation.requestInterruption()

    def _error(self, process, error):
        if process is self.process and error == QProcess.ProcessError.FailedToStart:
            self._drop_process(process)
            self._settle(
                "start_failed", "Could not start Gmsh: " + process.errorString()
            )

    def _drop_process(self, process):
        if self.process is process:
            self.process = None
        process.deleteLater()

    def _finished(self, process, code, status):
        if process is not self.process:
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
                parse_info(info)
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
                ["mesh.geo", "-3", "-format", "msh2", "-o", "mesh.msh"],
                self._timeout,
            )
            return
        self.stage_changed.emit("Checking element Jacobians and boundary groups…")
        if self._reason:
            self._settle(self._reason, "Mesh validation cancelled before execution.")
            return
        worker = MeshValidation(self.run, self)
        self.validation = worker
        worker.finished.connect(lambda: self._validated(worker))
        worker.start()

    def _validated(self, worker):
        if worker is not self.validation:
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
            save_mesh_result(self.run, worker.result)
        except (OSError, ValueError, TypeError) as error:
            self._settle("write_failed", str(error))
            return
        self.last_result = (worker.result, self.run.root)
        self._settle()
        self.completed.emit(self.last_result)

    def _settle(self, kind=None, message=""):
        self.timer.stop()
        self.busy = False
        self.failure_kind = kind
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
