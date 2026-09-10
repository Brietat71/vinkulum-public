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

from . import analysis as model
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
        self.analysis = model.for_selection(source)
        self.analysis_name = self.analysis.Name
        source = model.validate(self.analysis)
        self.assembly = model.is_assembly(source)
        self.adapter = bridge
        if self.assembly:
            from . import assembly_capture

            self.adapter = assembly_capture
        self.preview_poses = []
        self.preview_static = []
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
            "Native Assembly solids and Revolute joints, gravity −Z. "
            "Joint frames come from the Assembly; all moving solids use this density."
            if self.assembly
            else "One rigid solid, one revolute joint, gravity −Z. "
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
            if not self.assembly:
                fields.addRow(label, row)
        self.threads = QtWidgets.QSpinBox()
        self.threads.setRange(1, 64)
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
        self.reopen_button = QtWidgets.QPushButton("Reopen last calculation")
        self.reopen_button.clicked.connect(self.reopen_capture)
        layout.addWidget(self.reopen_button)
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
        self.bind_inputs()
        App.addDocumentObserver(self)

    def bind_inputs(self):
        self.property_fields = {
            "Density": [self.density],
            "Duration": [self.duration],
            "Step": [self.step],
            "Threads": [self.threads],
            "Pivot": self.pivot,
            "Axis": self.axis,
        }
        if self.assembly:
            self.property_fields.pop("Pivot")
            self.property_fields.pop("Axis")
        self.refresh_inputs()
        for name, controls in self.property_fields.items():
            for index, control in enumerate(controls):
                component = index if name in ("Pivot", "Axis") else None
                control.valueChanged.connect(
                    lambda value, name=name, component=component: model.update(
                        self.analysis, name, value, component
                    )
                )

    def refresh_inputs(self, name=None):
        for property_name, controls in self.property_fields.items():
            if name is not None and name != property_name:
                continue
            value = getattr(self.analysis, property_name)
            values = (
                list(value)
                if property_name in ("Pivot", "Axis")
                else [value.Value if property_name in ("Duration", "Step") else value]
            )
            for control, value in zip(controls, values):
                previous = control.blockSignals(True)
                try:
                    control.setValue(value)
                    control.setToolTip(f"Stored document value: {value!r}")
                finally:
                    control.blockSignals(previous)
        self.reopen_button.setEnabled(
            bool(self.analysis.LastCapture) and self.job is None
        )

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

    def source_exists(self):
        document = App.listDocuments().get(self.document_name)
        return (
            document is not None and document.getObject(self.source_name) is self.source
        )

    def source_is_current(self):
        document = App.listDocuments().get(self.document_name)
        return (
            self.source_exists()
            and document.getObject(self.analysis_name) is self.analysis
            and self.analysis.Source is self.source
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
            request = self.adapter.capture(
                self.source,
                directory,
                **model.inputs(self.analysis),
            )
            prefs = App.ParamGet(PREFERENCES)
            prefs.SetString("EnginePython", str(interpreter))
            prefs.SetString("OutputRoot", str(output))
            job = bridge.Job(
                self.source,
                directory,
                request,
                interpreter,
                adapter=self.adapter if self.assembly else None,
            )
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
        except (
            Exception
        ) as error:  # noqa: BLE001 - report native errors at the process/UI boundary
            self.status.setText(str(error))

    def set_busy(self, busy):
        if self.closed:
            return
        for widget in (self.inputs, self.run_button, self.open_button):
            widget.setEnabled(not busy)
        self.reopen_button.setEnabled(not busy and bool(self.analysis.LastCapture))
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
        if self.closed or not self.source_is_current():
            return
        model.remember_capture(self.analysis, directory)
        self.directory, self.request, self.result = directory, request, result
        self.geometry_dirty = True
        self.slider.setRange(0, len(result["time_s"]) - 1)
        self.slider.setEnabled(True)
        self.play_button.setEnabled(True)
        self.slider.setValue(0)
        self.display_sample(0)
        self.status.setText(
            f"Captured calculation saved in {directory}. "
            "Playback uses captured inputs; document inputs configure the next run."
        )

    def reopen_capture(self):
        if not self.closed and self.job is None and self.analysis.LastCapture:
            self.load_capture(Path(self.analysis.LastCapture))

    def open_capture(self):
        if self.closed or self.job is not None:
            return
        name = QtWidgets.QFileDialog.getExistingDirectory(
            self.form, "Open calculation folder"
        )
        if name:
            self.load_capture(Path(name))

    def load_capture(self, directory):
        if self.closed or self.job is not None:
            return
        try:
            path = directory / "request.json"
            if path.stat().st_size > (2_000_000 if self.assembly else 100_000):
                raise ValueError("Captured request exceeds the input budget.")
            request = json.loads(path.read_text())
            result = self.adapter.admit(self.source, directory, request)
            self.stop_playback()
            self.remove_preview()
            self.completed(directory, request, result)
        except (
            Exception
        ) as error:  # noqa: BLE001 - report native errors at the process/UI boundary
            self.status.setText(str(error))

    def remove_preview(self, deleted=None):
        document = App.listDocuments().get(self.document_name)
        preview, poses, visibility = (
            self.preview,
            self.preview_poses,
            self.source_visibility,
        )
        static = self.preview_static
        self.preview = None
        self.preview_poses = []
        self.preview_static = []
        self.source_visibility = None
        if document is not None:
            for obj in [part for part, _, _ in poses] + static:
                if (
                    obj is not preview
                    and obj is not deleted
                    and document.getObject(obj.Name) is obj
                ):
                    document.removeObject(obj.Name)
            if (
                preview is not None
                and preview is not deleted
                and document.getObject(preview.Name) is preview
            ):
                document.removeObject(preview.Name)
            if visibility is not None and self.source_exists():
                self.source.Visibility = visibility

    def display_sample(self, index):
        if self.closed or self.result is None:
            return
        try:
            if self.geometry_dirty:
                if self.assembly:
                    snapshot, _ = self.adapter.snapshot(self.source)
                    current = bridge.sha(bridge.request_bytes(snapshot))
                    expected = self.request["source_state_sha256"]
                else:
                    current = bridge.signature(self.source)
                    expected = self.request["source_geometry_sha256"]
                if current != expected:
                    raise ValueError(
                        "Geometry changed. Run again before displaying captured motion."
                    )
                self.geometry_dirty = False
            if self.preview is None:
                if self.assembly:
                    self.preview = self.document.addObject(
                        "App::DocumentObjectGroup", "VinkulumPlayback"
                    )
                    _, shapes = self.adapter.snapshot(self.source)
                    ground_name = self.request["source_state"]["grounded_name"]
                    ground = self.document.addObject(
                        "Part::Feature", "VinkulumCapturedBase"
                    )
                    self.preview.addObject(ground)
                    ground.Label = self.document.getObject(ground_name).Label
                    ground.Shape = shapes[ground_name].copy()
                    self.preview_static.append(ground)
                    captures = [
                        (body, shapes[body["name"]]) for body in self.request["bodies"]
                    ]
                else:
                    self.preview = self.document.addObject(
                        "Part::Feature", "VinkulumPlayback"
                    )
                    captures = [(self.request, self.source.Shape.copy())]
                self.preview.Label = "Vinkulum captured motion"
                for body, shape in captures:
                    obj = self.preview
                    if self.assembly:
                        obj = self.document.addObject(
                            "Part::Feature", "VinkulumCapturedBody"
                        )
                        self.preview.addObject(obj)
                        obj.Label = body["label"]
                    obj.Shape = shape.copy()
                    placement = App.Placement(obj.Placement)
                    self.preview_poses.append(
                        (obj, placement, body["properties_si"]["centre_m"])
                    )
                    obj.ViewObject.ShapeColor = (0.92, 0.58, 0.18)
                self.original_placement = self.preview_poses[0][1]
                self.source_visibility = self.source.Visibility
                self.source.Visibility = False
                # Recompute only the independent copies. A document recompute can
                # execute pending source features and alter native joint frames.
                for obj in [
                    part for part, _, _ in self.preview_poses
                ] + self.preview_static:
                    obj.recompute()
                self.preview.recompute()
            for body_index, (obj, placement, centre) in enumerate(self.preview_poses):
                position = self.result["position_m"][index]
                rotation = self.result["rotation"][index]
                if self.assembly:
                    position, rotation = position[body_index], rotation[body_index]
                bridge.set_pose(obj, placement, centre, position, rotation)
            self.time_label.setText(
                f"Captured time: {self.result['time_s'][index]:.6g} s · native sample {index}"
            )
        except (
            Exception
        ) as error:  # noqa: BLE001 - report native errors at the process/UI boundary
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
        if self.closed:
            return
        if obj is self.analysis:
            if property_name == "Source" and obj.Source is not self.source:
                self.reject()
                return
            self.refresh_inputs(property_name)
        is_preview = (
            obj is self.preview
            or obj in self.preview_static
            or any(obj is part for part, _, _ in self.preview_poses)
        )
        if not is_preview and (
            self.assembly or property_name in ("Shape", "Placement", "Group")
        ):
            self.geometry_dirty = True

    def slotRecomputedDocument(self, document):
        if document.Name == self.document_name:
            self.geometry_dirty = True

    def slotStartSaveDocument(self, document, filename):
        if document.Name == self.document_name and not self.closed:
            self.stop_playback()
            self.remove_preview()

    def slotDeletedObject(self, obj):
        if (
            obj is self.preview
            or obj in self.preview_static
            or any(obj is part for part, _, _ in self.preview_poses)
        ):
            self.stop_playback()
            self.remove_preview(deleted=obj)
        elif (obj is self.source or obj is self.analysis) and not self.closed:
            self.dispose()
            Gui.Control.closeDialog()


class MotionCommand:
    def GetResources(self):
        return {
            "MenuText": "Motion analysis…",
            "ToolTip": "Calculate the selected solid or native Assembly with Vinkulum",
        }

    def Activated(self):
        selection = Gui.Selection.getSelection()
        open_analysis(selection[0] if len(selection) == 1 else None)

    def IsActive(self):
        return App.ActiveDocument is not None


def open_analysis(selected):
    try:
        if Gui.Control.activeDialog():
            raise ValueError("Finish or close the current FreeCAD task first.")
        if selected is None:
            raise ValueError(
                "Select one solid, native Assembly or motion analysis first."
            )
        panel = MotionPanel(selected)
        Gui.Control.showDialog(panel)
        return panel
    except Exception as error:  # noqa: BLE001 - native UI boundary
        QtWidgets.QMessageBox.information(Gui.getMainWindow(), "Vinkulum", str(error))
        return None


def _open_example(filename, selected_name):
    path = Path(__file__).parent / "Examples" / filename
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
    Gui.Selection.addSelection(document.getObject(selected_name))


def open_example():
    _open_example("Pendulum.FCStd", "Rod")


def open_assembly_example():
    Gui.activateWorkbench("AssemblyWorkbench")
    _open_example("DoublePendulum.FCStd", "Assembly")


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
    menu.addAction("Open double-pendulum assembly", open_assembly_example)

    def attach_menu(*_):
        menubar = Gui.getMainWindow().menuBar()
        if menu.menuAction() not in menubar.actions():
            menubar.addMenu(menu)

    attach_menu()
    Gui.getMainWindow().workbenchActivated.connect(attach_menu)
    _installed = True
