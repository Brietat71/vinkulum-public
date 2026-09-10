"""One short-lived idle CAD process per application; active CPU leases stay separate."""

import hashlib
import json
import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from PySide6.QtCore import QCoreApplication, QObject, QProcess, QTimer, Slot


def process_key(executable, arguments, environment, threads):
    """Invalidate reuse after interpreter, source, package or environment changes."""

    def stamp(path):
        try:
            stat = path.stat()
            return str(path), stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns
        except OSError:
            return str(path), None

    packages = {}
    for name in ("vinkulum-studio", "build123d", "cadquery-ocp-novtk"):
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = None  # Frozen bundles are identified by their executable.
    root = Path(__file__).parent
    identity = [
        executable,
        arguments,
        str(Path.cwd()),
        sorted(environment.toStringList()),
        threads,
        packages,
        stamp(Path(executable)),
        [stamp(path) for path in sorted(root.glob("*.py"))],
    ]
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()


class CadProcessCache(QObject):
    def __init__(self, parent=None, *, idle_ms=15000):
        super().__init__(parent)
        if type(idle_ms) is not int or not 0 <= idle_ms <= 30000:
            raise ValueError("CAD idle retention must be between 0 and 30000 ms.")
        self.idle_ms = idle_ms
        self.idle_process = None
        self._key = None
        self._retiring = set()
        self._closing = False
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.expire)

    def take(self, key, parent):
        process = self.idle_process
        if process is None:
            return None
        if (
            self._closing
            or self._key != key
            or process.state() != QProcess.ProcessState.Running
            or process.bytesAvailable()
        ):
            self.expire()
            return None
        self.timer.stop()
        self.idle_process, self._key = None, None
        self._disconnect(process)
        process.readAllStandardError()
        process.setParent(parent)
        return process

    def retain(self, process, key):
        self.expire()
        process.setParent(self)
        process.finished.connect(self._ended)
        process.errorOccurred.connect(self._error)
        process.readyReadStandardOutput.connect(self._noise)
        process.readyReadStandardError.connect(self._stderr)
        if (
            self._closing
            or not self.idle_ms
            or process.state() != QProcess.ProcessState.Running
        ):
            self._retire(process)
            return
        self.idle_process, self._key = process, key
        self.timer.start(self.idle_ms)

    def _disconnect(self, process):
        process.finished.disconnect(self._ended)
        process.errorOccurred.disconnect(self._error)
        process.readyReadStandardOutput.disconnect(self._noise)
        process.readyReadStandardError.disconnect(self._stderr)

    def _retire(self, process):
        self._retiring.add(process)
        if process.state() == QProcess.ProcessState.NotRunning:
            self._remove(process)
        else:
            process.kill()

    def _remove(self, process):
        self._disconnect(process)
        self._retiring.discard(process)
        process.deleteLater()

    @Slot()
    def expire(self):
        self.timer.stop()
        process = self.idle_process
        self.idle_process, self._key = None, None
        if process is not None:
            self._retire(process)

    @Slot(int, QProcess.ExitStatus)
    def _ended(self, code, status):
        process = self.sender()
        if process is self.idle_process:
            self.timer.stop()
            self.idle_process, self._key = None, None
        self._remove(process)

    @Slot()
    def _noise(self):
        process = self.sender()
        process.readAllStandardOutput()
        if process is self.idle_process:
            self.expire()  # A quiescent worker must not emit unsolicited responses.

    @Slot()
    def _stderr(self):
        self.sender().readAllStandardError()

    @Slot(QProcess.ProcessError)
    def _error(self, error):
        if self.sender() is self.idle_process:
            self.expire()

    @Slot()
    def shutdown(self):
        if self._closing and not self._retiring:
            return
        self._closing = True
        self.expire()
        for process in tuple(self._retiring):
            process.kill()
            process.waitForFinished(2000)


def application_cad_cache():
    app = QCoreApplication.instance()
    if app is None:
        raise RuntimeError("CAD process reuse requires a Qt application.")
    if not hasattr(app, "_vinkulum_cad_cache"):
        try:
            seconds = int(os.environ.get("VINKULUM_CAD_IDLE_SECONDS", "15"))
            if not 0 <= seconds <= 30:
                raise ValueError
        except ValueError:
            raise ValueError(
                "VINKULUM_CAD_IDLE_SECONDS must be an integer from 0 to 30."
            ) from None
        app._vinkulum_cad_cache = CadProcessCache(app, idle_ms=1000 * seconds)
        app.aboutToQuit.connect(app._vinkulum_cad_cache.shutdown)
        app.destroyed.connect(app._vinkulum_cad_cache.shutdown)
    return app._vinkulum_cad_cache
