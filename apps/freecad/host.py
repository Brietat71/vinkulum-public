"""Native FreeCAD task panel for the isolated Vinkulum mechanics process.

SPDX-License-Identifier: Apache-2.0
"""

import json
import os
import time
import uuid
from pathlib import Path

import FreeCAD as App
import FreeCADGui as Gui
from PySide6 import QtCore, QtGui, QtWidgets

from . import bridge

_active_job = None
_installed = False
PREFERENCES = "User parameter:BaseApp/Preferences/Mod/Vinkulum"


def selected_solid():
    selection = Gui.Selection.getSelection()
    if len(selection) != 1 or not hasattr(selection[0], "Shape"):
        raise ValueError("Select one solid or PartDesign Body in the document first.")
    return selection[0]


class CompactNumber(QtWidgets.QDoubleSpinBox):
    def textFromValue(self, value):
        text = super().textFromValue(value)
        if self.locale().decimalPoint() in text:
            text = text.rstrip(self.locale().zeroDigit()).removesuffix(
                self.locale().decimalPoint()
            )
        return text


def number(value, minimum, maximum, decimals=6):
    field = CompactNumber()
    field.setDecimals(decimals)
    field.setRange(minimum, maximum)
    field.setValue(value)
    field.setKeyboardTracking(False)
    field.setMinimumWidth(45)
    field.setSizePolicy(
        QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Fixed
    )
    return field


class MotionPanel:
    """Inputs for the next calculation and an explicitly captured native result."""

    def __init__(self, source):
        self.source = source
        self.source_name = source.Name
        self.document = source.Document
        self.document_name = self.document.Name
        self.closed = False
        self.job = None
        self.result = None
        self.request = None
        self.directory = None
        self.preview = None
        self.source_visibility = None
        self.geometry_dirty = True
        self.form = QtWidgets.QWidget()
        self.form.setWindowTitle("Vinkulum · rigid motion")
        layout = QtWidgets.QVBoxLayout(self.form)
        title = QtWidgets.QLabel(source.Label)
        title.setWordWrap(True)
        layout.addWidget(title)
        scope = QtWidgets.QLabel(
            "One rigid solid, one revolute joint, gravity −Z. "
            "Pivot coordinates use the FreeCAD world frame in mm."
        )
        scope.setWordWrap(True)
        layout.addWidget(scope)
        self.inputs = QtWidgets.QGroupBox("Inputs for next run")
        fields = QtWidgets.QFormLayout(self.inputs)
        fields.setRowWrapPolicy(QtWidgets.QFormLayout.RowWrapPolicy.WrapLongRows)
        self.density = number(7800, 0.000001, 1e8)
        self.duration = number(2, 0.000001, 10)
        self.step = number(0.005, 0.000001, 10)
        fields.addRow("Density [kg/m³]", self.density)
        fields.addRow("Duration [s]", self.duration)
        fields.addRow("Native step [s]", self.step)
        self.pivot = [number(0, -1e9, 1e9) for _ in range(3)]
        self.axis = [number(v, -1e6, 1e6) for v in (0, 1, 0)]
        for label, controls in (
            ("Pivot [mm]", self.pivot),
            ("Axis direction", self.axis),
        ):
            row = QtWidgets.QWidget()
            boxes = QtWidgets.QHBoxLayout(row)
            boxes.setContentsMargins(0, 0, 0, 0)
            for axis, control in zip("XYZ", controls):
                control.setPrefix(axis + " ")
                boxes.addWidget(control)
            fields.addRow(label, row)
        self.threads = QtWidgets.QSpinBox()
        self.threads.setRange(1, min(64, os.cpu_count() or 1))
        self.threads.setValue(min(2, self.threads.maximum()))
        fields.addRow("Engine threads", self.threads)
        layout.addWidget(self.inputs)
        settings = QtWidgets.QGroupBox("Vinkulum engine and saved calculations")
        settings_layout = QtWidgets.QFormLayout(settings)
        settings_layout.setRowWrapPolicy(
            QtWidgets.QFormLayout.RowWrapPolicy.WrapLongRows
        )
        prefs = App.ParamGet(PREFERENCES)
        self.interpreter = QtWidgets.QLineEdit(
            os.environ.get(
                "VINKULUM_FREECAD_PYTHON", prefs.GetString("EnginePython", "")
            )
        )
        self.output = QtWidgets.QLineEdit(
            os.environ.get(
                "VINKULUM_FREECAD_OUTPUT_ROOT",
                prefs.GetString(
                    "OutputRoot", str(Path(App.getUserAppDataDir()) / "VinkulumRuns")
                ),
            )
        )
        for label, control, callback in (
            ("Engine Python", self.interpreter, self.choose_engine),
            ("Save runs in", self.output, self.choose_output),
        ):
            row = QtWidgets.QWidget()
            boxes = QtWidgets.QHBoxLayout(row)
            boxes.setContentsMargins(0, 0, 0, 0)
            boxes.addWidget(control)
            button = QtWidgets.QPushButton("Browse…")
            button.clicked.connect(callback)
            boxes.addWidget(button)
            settings_layout.addRow(label, row)
        layout.addWidget(settings)
        buttons = QtWidgets.QHBoxLayout()
        self.run_button = QtWidgets.QPushButton("Run motion")
        self.run_button.clicked.connect(self.run)
        self.cancel_button = QtWidgets.QPushButton("Cancel calculation")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.cancel_button)
        layout.addLayout(buttons)
        self.open_button = QtWidgets.QPushButton("Open saved calculation…")
        self.open_button.clicked.connect(self.open_capture)
        layout.addWidget(self.open_button)
        self.play_button = QtWidgets.QPushButton("Play captured motion")
        self.play_button.setEnabled(False)
        self.play_button.clicked.connect(self.toggle_playback)
        layout.addWidget(self.play_button)
        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self.display_sample)
        layout.addWidget(self.slider)
        self.time_label = QtWidgets.QLabel("No captured motion")
        layout.addWidget(self.time_label)
        self.status = QtWidgets.QLabel(
            "Ready. Select the Vinkulum Python environment before running."
        )
        self.status.setWordWrap(True)
        self.status.setMinimumWidth(0)
        self.status.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Ignored, QtWidgets.QSizePolicy.Policy.Preferred
        )
        self.status.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        layout.addWidget(self.status)
        self.timer = QtCore.QTimer(self.form)
        self.timer.setInterval(16)
        self.timer.timeout.connect(self.tick)
        App.addDocumentObserver(self)

    def getStandardButtons(self):
        return QtWidgets.QDialogButtonBox.StandardButton.Close.value

    def isAllowedAlterDocument(self):
        return True

    def isAllowedAlterSelection(self):
        return True

    def choose_engine(self):
        name, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.form, "Vinkulum engine Python", self.interpreter.text()
        )
        if name:
            self.interpreter.setText(name)

    def choose_output(self):
        name = QtWidgets.QFileDialog.getExistingDirectory(
            self.form, "Save calculations in", self.output.text()
        )
        if name:
            self.output.setText(name)

    def source_is_current(self):
        document = App.listDocuments().get(self.document_name)
        return (
            document is not None and document.getObject(self.source_name) is self.source
        )

    def run(self):
        global _active_job
        if self.closed:
            return
        if _active_job is not None:
            self.status.setText(
                "A Vinkulum calculation is already active. Wait or cancel it first."
            )
            return
        try:
            if not self.source_is_current():
                raise ValueError("The selected solid no longer exists.")
            interpreter = Path(self.interpreter.text()).expanduser()
            if not interpreter.is_file() or not os.access(interpreter, os.X_OK):
                raise ValueError(
                    "Choose the executable Python from a Vinkulum CAD environment."
                )
            self.stop_playback()
            self.remove_preview()
            output = Path(self.output.text()).expanduser()
            if not output.is_absolute():
                raise ValueError("Choose an absolute directory for saved calculations.")
            output.mkdir(parents=True, exist_ok=True)
            directory = output / str(uuid.uuid4())
            request = bridge.capture(
                self.source,
                directory,
                density=self.density.value(),
                body_id=str(uuid.uuid4()),
                project_id=str(uuid.uuid4()),
                duration=self.duration.value(),
                step=self.step.value(),
                pivot_mm=[field.value() for field in self.pivot],
                axis_world=[field.value() for field in self.axis],
                threads=self.threads.value(),
            )
            prefs = App.ParamGet(PREFERENCES)
            prefs.SetString("EnginePython", str(interpreter))
            prefs.SetString("OutputRoot", str(output))
            job = bridge.Job(self.source, directory, request, interpreter)
            self.job = _active_job = job
            job.completed.connect(
                lambda result: self.completed(directory, request, result)
            )
            job.failed.connect(self.failed)
            job.settled.connect(lambda: self.settled(job))
            self.set_busy(True)
            self.status.setText(
                "Calculating captured geometry. FreeCAD remains available."
            )
            job.start()
        except Exception as error:  # noqa: BLE001 - report native errors at the process/UI boundary
            self.status.setText(str(error))

    def set_busy(self, busy):
        if self.closed:
            return
        for widget in (self.inputs, self.run_button, self.open_button):
            widget.setEnabled(not busy)
        self.cancel_button.setEnabled(busy)

    def cancel(self):
        if self.job is not None:
            self.job.cancel()

    def failed(self, message):
        if not self.closed:
            self.status.setText(message)

    def settled(self, job):
        global _active_job
        if _active_job is job:
            _active_job = None
        if self.job is job:
            self.job = None
        self.set_busy(False)
        job.deleteLater()

    def completed(self, directory, request, result):
        if self.closed:
            return
        self.directory, self.request, self.result = directory, request, result
        self.geometry_dirty = True
        self.slider.setRange(0, len(result["time_s"]) - 1)
        self.slider.setEnabled(True)
        self.play_button.setEnabled(True)
        self.slider.setValue(0)
        self.display_sample(0)
        self.status.setText(
            f"Captured calculation saved in {directory}. Original geometry retained."
        )

    def open_capture(self):
        if self.closed or self.job is not None:
            return
        name = QtWidgets.QFileDialog.getExistingDirectory(
            self.form, "Open calculation folder"
        )
        if name:
            self.load_capture(Path(name))

    def load_capture(self, directory):
        try:
            path = directory / "request.json"
            if path.stat().st_size > 100_000:
                raise ValueError("Captured request exceeds the input budget.")
            request = json.loads(path.read_text())
            result = bridge.admit(self.source, directory, request)
            self.stop_playback()
            self.remove_preview()
            self.completed(directory, request, result)
        except Exception as error:  # noqa: BLE001 - report native errors at the process/UI boundary
            self.status.setText(str(error))

    def remove_preview(self):
        document = App.listDocuments().get(self.document_name)
        if document is not None:
            if (
                self.preview is not None
                and document.getObject(self.preview.Name) is self.preview
            ):
                document.removeObject(self.preview.Name)
            if self.source_visibility is not None and self.source_is_current():
                self.source.Visibility = self.source_visibility
        self.preview = None
        self.source_visibility = None

    def display_sample(self, index):
        if self.closed or self.result is None:
            return
        try:
            if self.geometry_dirty:
                if (
                    bridge.signature(self.source)
                    != self.request["source_geometry_sha256"]
                ):
                    raise ValueError(
                        "Geometry changed. Run again before displaying captured motion."
                    )
                self.geometry_dirty = False
            if self.preview is None:
                self.preview = self.document.addObject(
                    "Part::Feature", "VinkulumPlayback"
                )
                self.preview.Label = "Vinkulum captured motion"
                self.preview.Shape = self.source.Shape.copy()
                self.original_placement = App.Placement(self.preview.Placement)
                self.source_visibility = self.source.Visibility
                self.source.Visibility = False
                self.preview.ViewObject.ShapeColor = (0.92, 0.58, 0.18)
            bridge.set_pose(
                self.preview,
                self.original_placement,
                self.request["properties_si"]["centre_m"],
                self.result["position_m"][index],
                self.result["rotation"][index],
            )
            self.time_label.setText(
                f"Captured time: {self.result['time_s'][index]:.6g} s · native sample {index}"
            )
        except Exception as error:  # noqa: BLE001 - report native errors at the process/UI boundary
            self.stop_playback()
            self.remove_preview()
            self.status.setText(str(error))

    def toggle_playback(self):
        if self.timer.isActive():
            self.stop_playback()
        elif self.result is not None:
            if self.slider.value() == self.slider.maximum():
                self.slider.setValue(0)
            self.play_started = (
                time.monotonic() - self.result["time_s"][self.slider.value()]
            )
            self.play_button.setText("Pause captured motion")
            self.timer.start()

    def tick(self):
        if self.closed:
            return
        import bisect

        index = min(
            len(self.result["time_s"]) - 1,
            max(
                0,
                bisect.bisect_right(
                    self.result["time_s"], time.monotonic() - self.play_started
                )
                - 1,
            ),
        )
        self.slider.setValue(index)
        if index == self.slider.maximum():
            self.stop_playback()

    def stop_playback(self):
        self.timer.stop()
        self.play_button.setText("Play captured motion")

    def dispose(self):
        if self.closed:
            return
        self.stop_playback()
        self.closed = True
        App.removeDocumentObserver(self)
        self.cancel()
        self.remove_preview()

    def reject(self):
        self.dispose()
        Gui.Control.closeDialog()
        return True

    def slotDeletedDocument(self, document):
        if document.Name == self.document_name:
            self.dispose()
            Gui.Control.closeDialog()

    def slotChangedObject(self, obj, property_name):
        if obj is not self.preview and property_name in ("Shape", "Placement", "Group"):
            self.geometry_dirty = True

    def slotRecomputedDocument(self, document):
        if document.Name == self.document_name:
            self.geometry_dirty = True

    def slotStartSaveDocument(self, document, filename):
        if document.Name == self.document_name and not self.closed:
            self.stop_playback()
            self.remove_preview()

    def slotDeletedObject(self, obj):
        if obj is self.preview:
            self.preview = None
            if self.source_visibility is not None and self.source_is_current():
                self.source.Visibility = self.source_visibility
            self.source_visibility = None
        elif obj is self.source and not self.closed:
            self.dispose()
            Gui.Control.closeDialog()


class MotionCommand:
    def GetResources(self):
        return {
            "MenuText": "Motion analysis…",
            "ToolTip": "Calculate the selected solid with Vinkulum",
        }

    def Activated(self):
        try:
            if Gui.Control.activeDialog():
                raise ValueError("Finish or close the current FreeCAD task first.")
            Gui.Control.showDialog(MotionPanel(selected_solid()))
        except Exception as error:  # noqa: BLE001 - report native errors at the process/UI boundary
            QtWidgets.QMessageBox.information(
                Gui.getMainWindow(), "Vinkulum", str(error)
            )

    def IsActive(self):
        return App.ActiveDocument is not None


def open_example():
    path = Path(__file__).parent / "Examples" / "Pendulum.FCStd"
    document = App.openDocument(str(path))
    view = Gui.activeDocument().activeView()
    animated = view.isAnimationEnabled()
    view.setAnimationEnabled(False)
    try:
        view.viewFront()
        # An explicit margin fits immediately. Animated fitting enters a nested
        # event loop that can outlive this viewer if its document is closed.
        view.fitAll(1.15)
    finally:
        view.setAnimationEnabled(animated)
    Gui.Selection.clearSelection()
    Gui.Selection.addSelection(document.getObject("Rod"))


def install():
    global _installed
    if _installed:
        return
    Gui.addCommand("Vinkulum_Motion", MotionCommand())
    menu = QtWidgets.QMenu("Vinkulum", Gui.getMainWindow())
    menu.setObjectName("VinkulumMenu")
    action = QtGui.QAction("Motion analysis…", menu)
    action.triggered.connect(lambda: Gui.runCommand("Vinkulum_Motion"))
    menu.addAction(action)
    menu.addAction("Open pendulum example", open_example)

    def attach_menu(*_):
        menubar = Gui.getMainWindow().menuBar()
        if menu.menuAction() not in menubar.actions():
            menubar.addMenu(menu)

    attach_menu()
    Gui.getMainWindow().workbenchActivated.connect(attach_menu)
    _installed = True
