"""CAD operation dialog. Native geometry runs outside the desktop process."""

import sys
import tempfile
from dataclasses import asdict
from pathlib import Path

from PySide6.QtCore import QProcess, QTimer, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
)

from .document import MAX_PROJECT_BYTES
from .model import read_json, write_json


class CadDialog(QDialog):
    completed = Signal(object)

    def __init__(self, project, selection, parent=None):
        super().__init__(parent)
        self.setWindowTitle("CAD design · OCCT 8 / build123d")
        self.setMinimumWidth(430)
        self.project = project
        self.process = None
        self.temporary = None
        self._log = b""
        self.input_path = None
        self.output_path = None
        self.result_data = None
        layout = QVBoxLayout(self)
        title = QLabel("Create a part or modify a solid")
        title.setObjectName("section_title")
        layout.addWidget(title)
        self.form = QFormLayout()
        layout.addLayout(self.form)
        self.operation = QComboBox()
        for label, key in (
            ("Box", "box"),
            ("Cylinder", "cylinder"),
            ("Sphere", "sphere"),
            ("Extrude · XY rectangle", "extrude_rectangle"),
            ("Extrude · XY disk", "extrude_circle"),
            ("Subtract A − B", "cut"),
            ("Union A + B", "fuse"),
            ("Intersect A ∩ B", "common"),
            ("Fillet · all edges", "fillet"),
            ("Import STEP part", "import_step"),
            ("Export A to STEP", "export_step"),
        ):
            self.operation.addItem(label, key)
        self.form.addRow("Operation", self.operation)
        self.name = QLineEdit("CAD part")
        self.form.addRow("Name", self.name)
        self.a, self.b = QComboBox(), QComboBox()
        for body in project.bodies:
            self.a.addItem(body.name, body.id)
            self.b.addItem(body.name, body.id)
        index = self.a.findData(selection)
        if index >= 0:
            self.a.setCurrentIndex(index)
        if self.b.count() > 1:
            self.b.setCurrentIndex((self.a.currentIndex() + 1) % self.b.count())
        self.form.addRow("Body A · modified", self.a)
        self.form.addRow("Body B · retained", self.b)
        self.dimensions = [self.number(v, 0.001, 1e6, " mm") for v in (100, 60, 20)]
        for label, field in zip(("Length", "Width", "Height"), self.dimensions):
            self.form.addRow(label, field)
        self.position = [self.number(0, -1e6, 1e6, " mm") for _ in range(3)]
        for axis, field in zip("XYZ", self.position):
            self.form.addRow("Position " + axis, field)
        self.density = self.number(7800, 0.001, 1e6, " kg/m³")
        self.form.addRow("Density", self.density)
        self.radius = self.number(1, 0.001, 1e6, " mm")
        self.form.addRow("Fillet radius", self.radius)
        self.explanation = QLabel("")
        self.explanation.setWordWrap(True)
        self.explanation.setObjectName("muted")
        layout.addWidget(self.explanation)
        self.status = QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel)
        self.apply_button = QPushButton("Create part")
        self.apply_button.setObjectName("primary")
        self.buttons.addButton(
            self.apply_button, QDialogButtonBox.ButtonRole.ActionRole
        )
        self.apply_button.clicked.connect(self.start)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self._timeout)
        self.operation.currentIndexChanged.connect(self._fields)
        self.a.currentIndexChanged.connect(self._density_from_body)
        self._fields()

    @staticmethod
    def number(value, minimum, maximum, suffix):
        field = QDoubleSpinBox()
        field.setDecimals(4)
        field.setRange(minimum, maximum)
        field.setValue(value)
        field.setSuffix(suffix)
        field.setKeyboardTracking(False)
        return field

    def _density_from_body(self):
        body = next(
            (b for b in self.project.bodies if b.id == self.a.currentData()), None
        )
        if body and body.cad:
            self.density.setValue(body.mass / body.cad.volume_m3)

    def _fields(self):
        operation = self.operation.currentData()
        edit = operation in ("cut", "fuse", "common", "fillet", "export_step")
        primitive = operation in (
            "box",
            "sphere",
            "cylinder",
            "extrude_rectangle",
            "extrude_circle",
        )
        labels = (
            ("Radius", "Height")
            if operation in ("sphere", "cylinder", "extrude_circle")
            else ("Length", "Width", "Height")
        )
        count = {
            "box": 3,
            "sphere": 1,
            "cylinder": 2,
            "extrude_rectangle": 3,
            "extrude_circle": 2,
        }.get(operation, 0)
        for index, field in enumerate(self.dimensions):
            self.form.setRowVisible(field, index < count)
            if index < count:
                self.form.labelForField(field).setText(labels[index])
        for field in self.position:
            self.form.setRowVisible(field, primitive)
        self.form.setRowVisible(self.name, not edit)
        self.form.setRowVisible(self.a, edit)
        self.form.setRowVisible(self.b, operation in ("cut", "fuse", "common"))
        self.form.setRowVisible(self.radius, operation == "fillet")
        self.form.setRowVisible(self.density, operation != "export_step")
        self.apply_button.setText(
            "Export…"
            if operation == "export_step"
            else "Import…"
            if operation == "import_step"
            else "Apply to A"
            if edit
            else "Create part"
        )
        self.explanation.setText(
            "Body B is retained. The result must be a single solid. Attachments keep their world positions."
            if operation in ("cut", "fuse", "common")
            else "XY profile extruded toward +Z. Position specifies the centre of the starting profile."
            if operation.startswith("extrude")
            else "OCCT reads STEP units. Import is limited to one solid part and 8 MB."
            if operation == "import_step"
            else "Geometry is stored as BREP; mass and inertia are computed from the exact volume, then converted to SI."
        )
        if edit:
            self._density_from_body()

    def request(self):
        operation = self.operation.currentData()
        count = {
            "box": 3,
            "sphere": 1,
            "cylinder": 2,
            "extrude_rectangle": 3,
            "extrude_circle": 2,
        }.get(operation, 0)
        result = {
            "operation": operation,
            "name": self.name.text(),
            "density": self.density.value(),
            "dimensions_mm": [f.value() for f in self.dimensions[:count]],
            "position_mm": [f.value() for f in self.position],
            "radius_mm": self.radius.value(),
        }
        for key in ("a", "b"):
            body = next(
                (
                    body
                    for body in self.project.bodies
                    if body.id == getattr(self, key).currentData()
                ),
                None,
            )
            result[key] = asdict(body) if body else None
        return result

    def start(self):
        if self.process is not None:
            return
        request = self.request()
        operation = request["operation"]
        if operation == "import_step":
            path, _ = QFileDialog.getOpenFileName(
                self, "Import STEP part", "", "STEP (*.step *.stp)"
            )
            if not path:
                return
            request["path"] = path
        elif operation == "export_step":
            path, _ = QFileDialog.getSaveFileName(
                self, "Export STEP part", "", "STEP (*.step *.stp)"
            )
            if not path:
                return
            self.output_path = Path(path)
        self.temporary = tempfile.TemporaryDirectory(prefix="vinkulum-cad-")
        directory = Path(self.temporary.name)
        write_json(directory / "input.json", request)
        self.process = process = QProcess(self)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        process.readyReadStandardOutput.connect(self._drain)
        process.finished.connect(self._finish)
        process.errorOccurred.connect(self._error)
        arguments = (
            ["--cad-worker"]
            if getattr(sys, "frozen", False)
            else ["-m", "vinkulum_studio.cad_worker"]
        )
        self.apply_button.setEnabled(False)
        self.operation.setEnabled(False)
        self.status.setText("CAD operation in progress…")
        self.timer.start(60000)
        process.start(
            sys.executable,
            [*arguments, str(directory / "input.json"), str(directory / "output.json")],
        )

    def _drain(self):
        if self.process:
            self._log = (self._log + bytes(self.process.readAllStandardOutput()))[
                -8192:
            ]

    def _error(self, error):
        if error == QProcess.ProcessError.FailedToStart:
            self._finish(-1, QProcess.ExitStatus.CrashExit)

    def _timeout(self):
        if self.process:
            self.status.setText("CAD timeout (60 s). Document preserved.")
            self.process.kill()

    def _finish(self, code, status):
        if self.process is None:
            return
        self.timer.stop()
        process, self.process = self.process, None
        try:
            path = Path(self.temporary.name) / "output.json"
            result = read_json(path, MAX_PROJECT_BYTES) if path.exists() else {}
            if (
                code != 0
                or status != QProcess.ExitStatus.NormalExit
                or result.get("status") != "completed"
            ):
                raise ValueError(
                    result.get("error")
                    or "CAD operation interrupted. Document preserved."
                )
            self.result_data = result
            self.completed.emit(result)
            self.accept()
        except (OSError, ValueError) as error:
            self.status.setText(str(error))
        finally:
            process.deleteLater()
            self.temporary.cleanup()
            self.temporary = None
            self.apply_button.setEnabled(True)
            self.operation.setEnabled(True)

    def reject(self):
        if self.process:
            self.process.kill()
            self.process.waitForFinished(2000)
        super().reject()
