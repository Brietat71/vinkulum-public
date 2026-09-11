"""Experimental single-solid FreeCAD capture and native-motion playback.

SPDX-License-Identifier: Apache-2.0
Loaded inside FreeCAD; never imports Studio, numpy or a second OCCT binding.
"""

import copy
import hashlib
import json
import math
import uuid
from pathlib import Path

import FreeCAD as App
import FreeCADGui as Gui
import Part
from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, QTimer, Signal


def sha(data):
    return hashlib.sha256(data).hexdigest()


def request_bytes(request):
    return (json.dumps(request, indent=2, allow_nan=False) + "\n").encode("ascii")


def matrix_values(matrix):
    return [getattr(matrix, f"A{i}{j}") for i in range(1, 4) for j in range(1, 4)]


def single_solid(shape):
    """Unwrap containers without silently discarding extra geometry."""
    if not shape.isValid():
        raise ValueError("One valid positive-volume solid is required")
    solid = shape
    while solid.ShapeType in ("Compound", "CompSolid"):
        children = solid.childShapes()
        if len(children) != 1:
            raise ValueError("One solid without additional geometry is required")
        solid = children[0]
    if solid.ShapeType != "Solid" or solid.Volume <= 0:
        raise ValueError("One valid positive-volume solid is required")
    return solid


def signature(source):
    placement = source.getGlobalPlacement().toMatrix()
    pose = [getattr(placement, f"A{i}{j}") for i in range(1, 5) for j in range(1, 5)]
    return sha(
        source.Shape.exportBrepToString().encode("ascii")
        + json.dumps(pose, allow_nan=False).encode("ascii")
    )


def capture(
    source,
    directory,
    *,
    density,
    body_id,
    project_id,
    duration=2,
    step=0.005,
    pivot_mm=(0, 0, 0),
    axis_world=(0, 1, 0),
    threads=2,
):
    """Capture one top-level solid and an explicit world-frame revolute joint."""
    if not math.isfinite(density) or density <= 0:
        raise ValueError("A finite positive density in kg/m³ is required")
    if not (0 < step <= duration <= 10 and math.ceil(duration / step) <= 2000):
        raise ValueError("Prototype limit: 10 seconds and 2001 native samples")
    for vector in (pivot_mm, axis_world):
        if len(vector) != 3 or not all(math.isfinite(v) for v in vector):
            raise ValueError("The pivot and axis require three finite coordinates")
    norm = math.hypot(*axis_world)
    if norm < 1e-12:
        raise ValueError("The revolute axis must be nonzero")
    if type(threads) is not int or not 1 <= threads <= 64:
        raise ValueError("Select between 1 and 64 engine threads")
    local = source.Placement.toMatrix()
    global_placement = source.getGlobalPlacement().toMatrix()
    if any(
        abs(getattr(local, f"A{i}{j}") - getattr(global_placement, f"A{i}{j}")) > 1e-12
        for i in range(1, 5)
        for j in range(1, 5)
    ):
        raise ValueError("Select a top-level solid or PartDesign Body")
    shape = source.Shape
    solid = single_solid(shape)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    step_path = directory / "part.step"
    shape.exportStep(str(step_path))
    centre = list(solid.CenterOfMass)
    request = {
        "format": "vinkulum-freecad-prototype-1",
        "run_id": str(uuid.uuid4()),
        "body_id": body_id,
        "project_id": project_id,
        "label": source.Label,
        "source_name": source.Name,
        "source_geometry_sha256": signature(source),
        "step_sha256": sha(step_path.read_bytes()),
        "freecad_version": App.Version(),
        "freecad_occt": Part.OCC_VERSION,
        "density_kg_m3": density,
        "properties_si": {
            "volume_m3": solid.Volume * 1e-9,
            "mass_kg": density * solid.Volume * 1e-9,
            "centre_m": [v * 1e-3 for v in centre],
            "inertia_kg_m2": [
                v * density * 1e-15 for v in matrix_values(solid.MatrixOfInertia)
            ],
        },
        "characteristic_length_m": shape.BoundBox.DiagonalLength * 1e-3,
        "pivot_m": [v * 1e-3 for v in pivot_mm],
        "axis_world": [v / norm for v in axis_world],
        "threads": threads,
        "duration_s": duration,
        "step_s": step,
    }
    (directory / "request.json").write_bytes(request_bytes(request))
    return request


def admit(source, directory, request):
    """Bind the result to its capture, reject stale geometry and malformed poses."""
    captured = request_bytes(request)
    if (directory / "request.json").read_bytes() != captured:
        raise ValueError("Captured request changed during calculation")
    if signature(source) != request["source_geometry_sha256"]:
        raise ValueError(
            "Geometry changed during calculation; recapture before playback"
        )
    path = directory / "result/playback.json"
    if path.stat().st_size > 2_000_000:
        raise ValueError("Playback exceeds the prototype budget")
    result = json.loads(path.read_text())
    if (
        result.get("format") != "vinkulum-freecad-playback-1"
        or result.get("run_id") != request["run_id"]
        or result.get("body_id") != request["body_id"]
        or result.get("request_sha256") != sha(captured)
    ):
        raise ValueError("Playback does not identify the captured request")
    times, positions, rotations = (
        result[key] for key in ("time_s", "position_m", "rotation")
    )
    if not (2 <= len(times) <= 2001 and len(times) == len(positions) == len(rotations)):
        raise ValueError("Invalid native sample count")
    if times[0] != 0 or abs(times[-1] - request["duration_s"]) > 1e-9:
        raise ValueError("Incomplete native trajectory")
    for index, (t, position, rotation) in enumerate(zip(times, positions, rotations)):
        if len(position) != 3 or len(rotation) != 9:
            raise ValueError("Invalid pose dimensions")
        if not all(
            type(v) in (float, int) and math.isfinite(v)
            for v in [t, *position, *rotation]
        ):
            raise ValueError("Nonfinite native pose")
        if index and times[index - 1] >= t:
            raise ValueError("Native times must increase")
        for i in range(3):
            for j in range(3):
                dot = sum(rotation[3 * k + i] * rotation[3 * k + j] for k in range(3))
                if abs(dot - (i == j)) > 1e-8:
                    raise ValueError("Native rotation is not orthogonal")
        determinant = sum(
            rotation[i]
            * (
                rotation[3 + (i + 1) % 3] * rotation[6 + (i + 2) % 3]
                - rotation[3 + (i + 2) % 3] * rotation[6 + (i + 1) % 3]
            )
            for i in range(3)
        )
        if abs(determinant - 1) > 1e-8:
            raise ValueError("Native rotation is not proper")
    if any(
        abs(a - b) > 1e-9
        for a, b in zip(positions[0], request["properties_si"]["centre_m"])
    ):
        raise ValueError("Initial centre differs from the captured solid")
    if any(
        abs(a - b) > 1e-9 for a, b in zip(rotations[0], (1, 0, 0, 0, 1, 0, 0, 0, 1))
    ):
        raise ValueError("Initial body axes differ from the captured solid")
    return result


class Job(QObject):
    completed = Signal(object)
    failed = Signal(str)
    settled = Signal()

    def __init__(
        self, source, directory, request, interpreter, parent=None, *, adapter=None
    ):
        super().__init__(parent or Gui.getMainWindow())
        self.source, self.directory = source, Path(directory)
        self.request = copy.deepcopy(request)
        self.admit_result = adapter.admit if adapter is not None else admit
        worker = "assembly_worker.py" if adapter is not None else "worker.py"
        self.process = QProcess(self)
        self.diagnostics = bytearray()
        self.timed_out = False
        self.cancelled = False
        self.terminal = False
        environment = QProcessEnvironment.systemEnvironment()
        # The official AppImage exports PYTHONHOME for its own Python 3.11.
        # The explicitly selected Vinkulum interpreter must use its own runtime.
        for key in (
            "PYTHONHOME",
            "PYTHONPATH",
            "LD_LIBRARY_PATH",
            "QT_PLUGIN_PATH",
            "QT_QPA_PLATFORM_PLUGIN_PATH",
        ):
            environment.remove(key)
        environment.insert("OPENBLAS_NUM_THREADS", "1")
        environment.insert("OMP_NUM_THREADS", "1")
        self.process.setProcessEnvironment(environment)
        self.process.setProgram(str(interpreter))
        self.process.setArguments(
            [str(Path(__file__).with_name(worker)), str(self.directory)]
        )
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read)
        self.process.finished.connect(self._finished)
        self.process.errorOccurred.connect(self._error)
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self._timeout)

    def start(self):
        if self.terminal:
            raise RuntimeError("A settled job cannot be restarted")
        self.timer.start(60000)
        self.process.start()

    def cancel(self):
        if self.terminal:
            return
        self.cancelled = True
        self.timer.stop()
        if self.process.state() == QProcess.ProcessState.NotRunning:
            self._settle(error="Calculation cancelled")
        else:
            self.process.kill()

    def _settle(self, *, error=None, result=None):
        if self.terminal:
            return
        self.terminal = True
        self.timer.stop()
        try:
            if error is not None:
                self.failed.emit(error)
            else:
                self.completed.emit(result)
        finally:
            self.settled.emit()

    def _read(self):
        self.diagnostics.extend(bytes(self.process.readAllStandardOutput()))
        del self.diagnostics[:-8192]

    def _error(self, error):
        if error == QProcess.ProcessError.FailedToStart:
            self._settle(error=self.process.errorString())

    def _timeout(self):
        self.timed_out = True
        self.process.kill()

    def _finished(self, code, status):
        if self.terminal:
            return
        self.timer.stop()
        self._read()
        try:
            (self.directory / "worker.log").write_bytes(self.diagnostics)
        except OSError as error:
            self._settle(error=f"Cannot retain worker diagnostics: {error}")
            return
        if self.cancelled:
            self._settle(error="Calculation cancelled")
            return
        if self.timed_out or code != 0 or status != QProcess.ExitStatus.NormalExit:
            self._settle(
                error="Worker failed or timed out: "
                + self.diagnostics.decode(errors="replace")
            )
            return
        try:
            result = self.admit_result(self.source, self.directory, self.request)
        except Exception as error:  # noqa: BLE001 - report native errors at the process/UI boundary
            self._settle(error=str(error))
            return
        self._settle(result=result)


def set_pose(playback, original_placement, original_centre_m, position_m, rotation):
    """Apply the native body-to-world transform to a separate captured shape."""
    matrix = App.Matrix()
    for i in range(3):
        for j in range(3):
            setattr(matrix, f"A{i + 1}{j + 1}", rotation[3 * i + j])
        offset = position_m[i] - sum(
            rotation[3 * i + j] * original_centre_m[j] for j in range(3)
        )
        setattr(matrix, f"A{i + 1}4", 1000 * offset)
    playback.Placement = App.Placement(matrix).multiply(original_placement)
