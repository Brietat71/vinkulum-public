"""Linux process-group lifetime for the external CAD/mesh/static worker.

SPDX-License-Identifier: Apache-2.0
"""

import os
import signal
from pathlib import Path

from PySide6.QtCore import QProcess

from .bridge import Job


class StaticJob(Job):
    def __init__(self, source, directory, request, interpreter, gmsh, ccx):
        super().__init__(source, directory, request, interpreter)
        self.group = None
        self.process.setUnixProcessParameters(QProcess.UnixProcessFlag.CreateNewSession)
        self.process.setArguments(
            [
                str(Path(__file__).with_name("static_guard.py")),
                str(os.getpid()),
                str(directory),
                "--gmsh",
                str(gmsh),
                "--ccx",
                str(ccx),
            ]
        )
        self.process.started.connect(self.started)

    def started(self):
        self.group = int(self.process.processId())
        if self.cancelled:
            self.kill_group()

    def start(self):
        super().start()
        self.timer.start(240000)

    def kill_group(self):
        if self.group is not None:
            try:
                os.killpg(self.group, signal.SIGKILL)
            except ProcessLookupError:
                pass
        elif self.process.state() == QProcess.ProcessState.Running:
            self.process.kill()

    def cancel(self):
        if self.terminal:
            return
        self.cancelled = True
        self.timer.stop()
        self.kill_group()
        if self.process.state() == QProcess.ProcessState.NotRunning:
            self._settle(error="Calculation cancelled")

    def _timeout(self):
        self.timed_out = True
        self.kill_group()

    def _finished(self, code, status):
        # The worker owns all Gmsh and CalculiX descendants in this session.
        self.kill_group()
        self.group = None
        self._read()
        try:
            (self.directory / "worker.log").write_bytes(self.diagnostics)
        except OSError as error:
            self._settle(error=str(error))
            return
        if self.cancelled:
            self._settle(error="Calculation cancelled")
        elif self.timed_out or code or status != QProcess.ExitStatus.NormalExit:
            self._settle(
                error="Static worker failed or timed out: "
                + self.diagnostics.decode(errors="replace")
            )
        else:
            self._settle(result=self.directory)
