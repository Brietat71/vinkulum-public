"""A supervised CAD transaction, with distinct failure causes and captured input."""

import hashlib
import json
import math
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
from PySide6.QtCore import QObject, QProcess, QTimer, Signal

from .cad_history import CadRecipe
from .document import MAX_PROJECT_BYTES, Body
from .model import finite_number, read_json, write_json
from .sketch import Sketch


class CadController(QObject):
    completed = Signal(object)
    failed = Signal(str, str)
    busy_changed = Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.process = None
        self.temporary = None
        self.last_result = None
        self.last_request = None
        self.last_log = ""
        self.failure_kind = None
        self._reason = None
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(lambda: self.cancel("timed_out"))

    def start(self, request, *, timeout=60, executable=None, arguments=None):
        if self.process is not None:
            raise RuntimeError("A CAD operation is already running.")
        if not finite_number(timeout) or not 0 < timeout <= 120:
            raise ValueError("CAD timeout must be between 0 and 120 seconds.")
        captured = json.loads(json.dumps(request, allow_nan=False))
        if not isinstance(captured, dict):
            raise TypeError("A CAD operation request is required.")
        self.temporary = tempfile.TemporaryDirectory(prefix="vinkulum-cad-")
        root = Path(self.temporary.name)
        try:
            write_json(root / "input.json", captured)
        except OSError, ValueError, TypeError:
            self.temporary.cleanup()
            self.temporary = None
            raise
        self._input_hash = hashlib.sha256(
            (root / "input.json").read_bytes()
        ).hexdigest()
        self.last_request = captured
        self.last_log, self.failure_kind, self._reason = "", None, None
        process = QProcess(self)
        self.process = process
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.setStandardOutputFile(str(root / "worker.log"))
        process.finished.connect(
            lambda code, status: self._finished(process, code, status)
        )
        process.errorOccurred.connect(lambda error: self._error(process, error))
        self.busy_changed.emit(True)
        self.timer.start(max(1, round(timeout * 1000)))
        if arguments is None:
            arguments = (
                ["--cad-worker"]
                if getattr(sys, "frozen", False)
                else ["-m", "vinkulum_studio.cad_worker"]
            )
        process.start(
            str(executable or sys.executable),
            [*arguments, str(root / "input.json"), str(root / "output.json")],
        )

    def cancel(self, reason="cancelled"):
        if self.process is not None:
            self._reason = reason
            self.timer.stop()
            self.process.kill()

    def _error(self, process, error):
        if process is self.process and error == QProcess.ProcessError.FailedToStart:
            self._settle(
                process,
                "start_failed",
                "The CAD worker could not start: " + process.errorString(),
            )

    def _read_log(self):
        path = Path(self.temporary.name) / "worker.log"
        if path.is_file():
            with path.open("rb") as stream:
                stream.seek(0, 2)
                stream.seek(max(0, stream.tell() - 8192))
                self.last_log = stream.read(8192).decode("utf-8", errors="replace")

    def _finished(self, process, code, status):
        if process is not self.process:
            return
        try:
            self._read_log()
            if self._reason:
                message = (
                    "CAD timed out. The previous model is preserved."
                    if self._reason == "timed_out"
                    else "CAD cancelled. The previous model is preserved."
                )
                self._settle(process, self._reason, message)
                return
            if status == QProcess.ExitStatus.CrashExit:
                self._settle(
                    process,
                    "crashed",
                    "The CAD worker crashed. The previous model is preserved.",
                )
                return
            path = Path(self.temporary.name) / "output.json"
            result = read_json(path, MAX_PROJECT_BYTES) if path.is_file() else {}
            if not isinstance(result, dict):
                raise TypeError("Invalid CAD worker response.")
            if result.get("status") == "failed":
                self._settle(
                    process,
                    "operation_failed",
                    "CAD operation failed: "
                    + str(result.get("error", "No detail returned.")),
                )
                return
            if code != 0:
                self._settle(
                    process,
                    "worker_failed",
                    f"CAD worker exited with code {code}. The previous model is preserved.",
                )
                return
            if (
                result.get("status") != "completed"
                or result.get("request_sha256") != self._input_hash
            ):
                raise ValueError("CAD response does not match the captured operation.")
            if "body" in result:
                body = Body.from_dict(result["body"])
                if body.cad is None:
                    raise ValueError("The CAD worker returned no solid geometry.")
                density = self.last_request.get("density", 7800)
                if body.inertia_mode != "homogeneous" or not math.isclose(
                    body.mass, body.cad.volume_m3 * density, rel_tol=1e-12, abs_tol=0
                ):
                    raise ValueError(
                        "CAD mass properties do not match the captured density."
                    )
                operation = self.last_request.get("operation")
                if operation == "extrude_sketch":
                    features = body.cad.recipe.features if body.cad.recipe else ()
                    if (
                        len(features) != 1
                        or features[0].kind != "extrude_sketch"
                        or features[0].profile
                        != Sketch.from_dict(self.last_request["profile"])
                        or features[0].dimensions_mm
                        != tuple(self.last_request["dimensions_mm"])
                    ):
                        raise ValueError(
                            "The CAD worker changed the captured sketch or extrusion height."
                        )
                if (
                    operation in ("regenerate", "fillet", "cut", "fuse", "common")
                    and body.id != self.last_request["a"]["id"]
                ):
                    raise ValueError("The CAD worker changed the target body identity.")
                if operation == "regenerate":
                    expected = CadRecipe.from_dict(self.last_request["recipe"])
                    if (
                        body.cad.recipe is None
                        or replace(
                            body.cad.recipe, origin_in_body_m=expected.origin_in_body_m
                        )
                        != expected
                    ):
                        raise ValueError(
                            "The CAD worker returned a different feature graph."
                        )
                    source = Body.from_dict(self.last_request["a"])
                    if (
                        body.name != source.name
                        or body.orientation != source.orientation
                    ):
                        raise ValueError(
                            "CAD regeneration changed the captured body frame."
                        )
                    expected_position = np.array(source.position) + np.array(
                        source.orientation
                    ).reshape(3, 3) @ (
                        np.array(expected.origin_in_body_m)
                        - body.cad.recipe.origin_in_body_m
                    )
                    tolerance = 1e-10 + 16 * np.spacing(
                        np.maximum(np.abs(expected_position), np.abs(source.position))
                    )
                    if not np.all(
                        np.abs(np.array(body.position) - expected_position) <= tolerance
                    ):
                        raise ValueError(
                            "CAD regeneration changed the design frame's world origin."
                        )
            elif self.last_request.get("operation") == "export_step":
                if not isinstance(result.get("step"), str) or not result[
                    "step"
                ].startswith("ISO-10303-21;"):
                    raise ValueError("Invalid CAD STEP output.")
            elif self.last_request.get("operation") != "self_check":
                raise ValueError("The CAD worker returned no part.")
        except (OSError, ValueError, TypeError, KeyError) as error:
            self._settle(process, "invalid_result", str(error))
            return
        self.last_result = result
        self._settle(process)
        self.completed.emit(result)

    def _settle(self, process, kind=None, message=""):
        if process is not self.process:
            return
        self.timer.stop()
        self.process = None
        process.deleteLater()
        self.temporary.cleanup()
        self.temporary = None
        self.failure_kind = kind
        self.busy_changed.emit(False)
        if kind:
            self.failed.emit(kind, message)

    def shutdown(self):
        process = self.process
        if process is not None:
            self.cancel()
            process.waitForFinished(2000)
            if (
                process is self.process
                and process.state() == QProcess.ProcessState.NotRunning
            ):
                self._settle(
                    process, "cancelled", "CAD stopped when its window closed."
                )
