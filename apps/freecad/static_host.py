"""Native FreeCAD face selection, static inputs and asynchronous calculation.

SPDX-License-Identifier: Apache-2.0
"""

import os
import uuid
from pathlib import Path

import FreeCAD as App
import FreeCADGui as Gui
from PySide6 import QtCore, QtWidgets

from . import host, static_bridge
from . import static_analysis as model
from .analysis import transaction
from .static_job import StaticJob


class BoundaryList(QtWidgets.QListWidget):
    """Open keyboard context menus on the current boundary, independent of the mouse."""

    def keyPressEvent(self, event):
        if event.key() == QtCore.Qt.Key.Key_Menu or (
            event.key() == QtCore.Qt.Key.Key_F10
            and event.modifiers() == QtCore.Qt.KeyboardModifier.ShiftModifier
        ):
            event.accept()
            item = self.currentItem()
            if item is not None:
                self.scrollToItem(item)
                self.customContextMenuRequested.emit(self.visualItemRect(item).center())
            return
        super().keyPressEvent(event)


class StaticPanel:
    def __init__(self, selected):
        self.analysis = model.for_selection(selected)
        self.document = self.analysis.Document
        self.document_name = self.document.Name
        self.source = self.analysis.Source
        self.closed, self.job = False, None
        self.form = QtWidgets.QWidget()
        self.form.setWindowTitle("Vinkulum · linear static analysis")
        layout = QtWidgets.QVBoxLayout(self.form)
        scope = QtWidgets.QLabel(
            "Linear elasticity · fixed supports and pressure.\nSelect faces in the 3D view."
        )
        scope.setWordWrap(True)
        scope.setMinimumHeight(3 * scope.fontMetrics().lineSpacing())
        layout.addWidget(scope)
        self.controls = QtWidgets.QWidget()
        fields = QtWidgets.QFormLayout(self.controls)
        self.fields = {}
        for name, label, scale, minimum, maximum in (
            ("YoungPa", "Young modulus [MPa]", 1e6, 0.000001, 1e9),
            ("Poisson", "Poisson ratio", 1, -0.999999, 0.499999),
            ("MeshSizeMm", "Mesh size [mm]", 1, 0.000001, 1e6),
            ("Threads", "Engine threads", 1, 1, 64),
        ):
            control = host.number(
                getattr(self.analysis, name) / scale,
                minimum,
                maximum,
                0 if name == "Threads" else 6,
            )
            control.valueChanged.connect(
                lambda value, name=name, scale=scale: self.update(
                    name, int(value) if name == "Threads" else value * scale
                )
            )
            self.fields[name] = control
            fields.addRow(label, control)
        self.boundaries = BoundaryList()
        self.boundaries.setFixedHeight(80)
        self.boundaries.setContextMenuPolicy(
            QtCore.Qt.ContextMenuPolicy.CustomContextMenu
        )
        self.boundaries.customContextMenuRequested.connect(self.boundary_menu)
        fields.addRow("Boundary conditions", self.boundaries)
        self.fixed = QtWidgets.QPushButton("Fix selected faces")
        self.fixed.clicked.connect(lambda: self.add("Fixed"))
        fields.addRow(self.fixed)
        self.pressure = host.number(0, -1e9, 1e9)
        fields.addRow("Pressure [MPa] (+ inward)", self.pressure)
        self.add_pressure = QtWidgets.QPushButton("Apply pressure to selected faces")
        self.add_pressure.clicked.connect(lambda: self.add("Pressure"))
        fields.addRow(self.add_pressure)
        self.edit_pressure = QtWidgets.QPushButton("Update selected pressure")
        self.edit_pressure.clicked.connect(self.update_pressure)
        fields.addRow(self.edit_pressure)
        self.remove = QtWidgets.QPushButton("Remove selected condition")
        self.remove.clicked.connect(self.remove_boundary)
        self.boundaries.currentRowChanged.connect(self.select_boundary)
        fields.addRow(self.remove)
        layout.addWidget(self.controls)
        self.settings = QtWidgets.QGroupBox("Engine executables and saved calculations")
        paths = QtWidgets.QFormLayout(self.settings)
        prefs = App.ParamGet(host.PREFERENCES)
        self.paths = {}
        for name, label, fallback in (
            (
                "EnginePython",
                "Engine Python",
                os.environ.get("VINKULUM_FREECAD_PYTHON", ""),
            ),
            ("Gmsh", "OCCT 8 Gmsh", ""),
            ("Calculix", "CalculiX", ""),
            (
                "OutputRoot",
                "Save runs in",
                str(Path(App.getUserAppDataDir()) / "VinkulumRuns"),
            ),
        ):
            field = QtWidgets.QLineEdit(prefs.GetString(name, fallback))
            self.paths[name] = field
            row = QtWidgets.QWidget()
            boxes = QtWidgets.QHBoxLayout(row)
            boxes.setContentsMargins(0, 0, 0, 0)
            boxes.addWidget(field)
            button = QtWidgets.QPushButton("Browse…")
            button.clicked.connect(lambda _=False, name=name: self.browse(name))
            boxes.addWidget(button)
            paths.addRow(label, row)
        self.engine_toggle = QtWidgets.QToolButton()
        self.engine_toggle.setText("Configure engines and output folder")
        self.engine_toggle.setCheckable(True)
        self.engine_toggle.setChecked(True)
        self.engine_toggle.toggled.connect(self.settings.setVisible)
        layout.addWidget(self.engine_toggle)
        layout.addWidget(self.settings)
        self.run_button = QtWidgets.QPushButton("Run static analysis")
        self.run_button.clicked.connect(self.run)
        self.cancel_button = QtWidgets.QPushButton("Cancel calculation")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel)
        layout.addWidget(self.run_button)
        layout.addWidget(self.cancel_button)
        self.status = QtWidgets.QLabel("Ready. Negative pressure pulls outward.")
        self.status.setWordWrap(True)
        self.status.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        layout.addWidget(self.status)
        self.refresh()
        App.addDocumentObserver(self)

    def update(self, name, value):
        with transaction(self.document, "Edit static " + name):
            setattr(self.analysis, name, value)
        model.refresh_result(self.analysis)

    def browse(self, name):
        if name == "OutputRoot":
            value = QtWidgets.QFileDialog.getExistingDirectory(
                self.form, "Save calculations in"
            )
        else:
            value, _ = QtWidgets.QFileDialog.getOpenFileName(
                self.form, "Choose executable"
            )
        if value:
            self.paths[name].setText(value)

    def selected_boundary(self):
        item = self.boundaries.currentItem()
        if item is None:
            return None
        name = item.data(QtCore.Qt.ItemDataRole.UserRole)
        return next((row for row in self.analysis.Boundaries if row.Name == name), None)

    def select_boundary(self, _index=-1):
        row = self.selected_boundary()
        pressure = row is not None and row.Kind == "Pressure"
        self.edit_pressure.setEnabled(pressure)
        self.remove.setEnabled(row is not None)
        if pressure:
            self.pressure.setValue(row.PressurePa / 1e6)

    def refresh(self):
        current = self.boundaries.currentItem()
        name = current.data(QtCore.Qt.ItemDataRole.UserRole) if current else None
        blocker = QtCore.QSignalBlocker(self.boundaries)
        self.boundaries.clear()
        for row in self.analysis.Boundaries:
            item = QtWidgets.QListWidgetItem(
                row.Label
                + (f" · {row.PressurePa / 1e6:g} MPa" if row.Kind == "Pressure" else "")
            )
            item.setData(QtCore.Qt.ItemDataRole.UserRole, row.Name)
            self.boundaries.addItem(item)
            if row.Name == name:
                self.boundaries.setCurrentItem(item)
        del blocker
        self.select_boundary()

    def add(self, kind):
        try:
            refs = model.selected_faces(self.analysis, Gui.Selection.getSelectionEx())
            model.add_boundary(self.analysis, refs, kind, self.pressure.value() * 1e6)
            self.refresh()
            self.status.setText(
                "Boundary added. Select other faces for the next condition."
            )
        except Exception as error:  # noqa: BLE001 - native UI error boundary
            self.status.setText(str(error))

    def boundary_menu(self, point):
        if self.boundaries.itemAt(point) is None:
            return
        self.boundaries.setCurrentItem(self.boundaries.itemAt(point))
        menu = QtWidgets.QMenu(self.boundaries)
        row = self.selected_boundary()
        if row is not None and row.Kind == "Pressure":
            menu.addAction("Edit pressure…", self.focus_pressure)
        menu.addAction("Remove condition", self.remove_boundary)
        menu.exec(self.boundaries.viewport().mapToGlobal(point))

    def focus_pressure(self):
        self.select_boundary()
        self.pressure.setFocus()
        self.pressure.selectAll()

    def update_pressure(self):
        index = self.boundaries.currentRow()
        if (
            not 0 <= index < len(self.analysis.Boundaries)
            or self.analysis.Boundaries[index].Kind != "Pressure"
        ):
            self.status.setText("Select a pressure condition in the list first.")
            return
        with transaction(self.document, "Edit selected pressure"):
            self.analysis.Boundaries[index].PressurePa = self.pressure.value() * 1e6
        self.refresh()
        model.refresh_result(self.analysis)

    def remove_boundary(self):
        index = self.boundaries.currentRow()
        rows = list(self.analysis.Boundaries)
        if 0 <= index < len(rows):
            row = rows.pop(index)
            with transaction(self.document, "Remove static boundary"):
                self.analysis.Boundaries = rows
                self.document.removeObject(row.Name)
            self.refresh()
            model.refresh_result(self.analysis)

    def run(self):
        if self.closed:
            return
        if host._active_job is not None:
            self.status.setText(
                "A Vinkulum calculation is active. Wait or cancel it first."
            )
            return
        try:
            arguments = model.inputs(self.analysis)
            captured_inputs = model.signature(self.analysis)
            paths = {
                name: Path(field.text()).expanduser()
                for name, field in self.paths.items()
            }
            for name, path in paths.items():
                if not path.is_absolute() or (
                    name != "OutputRoot"
                    and (not path.is_file() or not os.access(path, os.X_OK))
                ):
                    raise ValueError(
                        f"Choose an absolute {'directory' if name == 'OutputRoot' else 'executable'} for {name}."
                    )
            directory = paths["OutputRoot"] / str(uuid.uuid4())
            request = static_bridge.capture(self.source, directory, **arguments)
            job = StaticJob(
                self.source,
                directory,
                request,
                paths["EnginePython"],
                paths["Gmsh"],
                paths["Calculix"],
            )
            self.job = host._active_job = job
            for name, path in paths.items():
                App.ParamGet(host.PREFERENCES).SetString(name, str(path))
            job.completed.connect(
                lambda _: self.completed(directory, request, captured_inputs)
            )
            job.failed.connect(self.failed)
            job.settled.connect(lambda: self.settled(job))
            self.engine_toggle.setChecked(False)
            self.set_busy(True)
            self.status.setText(
                "Meshing and calculating the captured solid… FreeCAD remains available."
            )
            job.start()
        except Exception as error:  # noqa: BLE001 - native UI error boundary
            self.status.setText(str(error))

    def completed(self, directory, request, captured_inputs):
        if self.closed:
            return
        try:
            if model.signature(self.analysis) != captured_inputs:
                raise ValueError(
                    "Inputs changed during calculation. Recalculate before importing."
                )
            result, _ = static_bridge.import_result(self.source, directory, request)
            if self.analysis.Result and self.analysis.Result.Mesh:
                self.analysis.Result.Mesh.Visibility = False
            with transaction(self.document, "Record static result"):
                self.analysis.Result = result
                self.analysis.LastCapture = str(directory)
                self.analysis.CapturedInputsSha256 = captured_inputs
            result.Mesh.ViewObject.setNodeColorByScalars(
                result.NodeNumbers, result.DisplacementLengths
            )
            result.Mesh.Visibility = True
            if not self.analysis.SourceHidden:
                self.analysis.OriginalSourceVisibility = self.source.Visibility
                self.analysis.SourceHidden = True
            self.source.Visibility = False
            model.refresh_result(self.analysis)
            self.status.setText(
                "Captured displacement result saved in the document. No general FEM error bound is claimed."
            )
        except Exception as error:  # noqa: BLE001 - native UI error boundary
            self.status.setText(str(error))

    def set_busy(self, busy):
        if not self.closed:
            self.controls.setEnabled(not busy)
            self.settings.setEnabled(not busy)
            self.run_button.setEnabled(not busy)
            self.cancel_button.setEnabled(busy)

    def settled(self, job):
        if host._active_job is job:
            host._active_job = None
        if self.job is job:
            self.job = None
        self.set_busy(False)
        job.deleteLater()

    def failed(self, error):
        if not self.closed:
            self.status.setText(error)

    def cancel(self):
        if self.job:
            self.job.cancel()

    def dispose(self):
        if not self.closed:
            self.closed = True
            App.removeDocumentObserver(self)
            self.cancel()

    def getStandardButtons(self):
        return QtWidgets.QDialogButtonBox.StandardButton.Close.value

    def accept(self):
        return self.reject()

    def reject(self):
        self.dispose()
        Gui.Control.closeDialog()
        return True

    def slotDeletedDocument(self, document):
        if document.Name == self.document_name:
            self.reject()

    def slotDeletedObject(self, obj):
        if obj is self.source or obj is self.analysis:
            self.reject()

    def slotChangedObject(self, obj, name):
        if self.closed:
            return
        if obj is self.analysis and name in self.fields:
            control = self.fields[name]
            scale = 1e6 if name == "YoungPa" else 1
            blocker = QtCore.QSignalBlocker(control)
            control.setValue(getattr(obj, name) / scale)
            del blocker
            model.refresh_result(self.analysis)
        if (
            obj is self.analysis and name == "Boundaries"
        ) or obj in self.analysis.Boundaries:
            self.refresh()
            model.refresh_result(self.analysis)
        if obj is self.analysis and name == "Source":
            self.reject()
        elif obj is self.source and name in ("Shape", "Placement"):
            model.refresh_result(self.analysis)

    def slotRecomputedDocument(self, document):
        if not self.closed and document.Name == self.document_name:
            model.refresh_result(self.analysis)


class StaticCommand:
    def GetResources(self):
        return {
            "MenuText": "Static analysis…",
            "ToolTip": "Select faces and calculate linear elasticity with CalculiX",
        }

    def IsActive(self):
        return App.ActiveDocument is not None

    def Activated(self):
        selected = Gui.Selection.getSelection()
        open_analysis(selected[0] if len(selected) == 1 else None)


def open_analysis(selected):
    try:
        if Gui.Control.activeDialog():
            raise ValueError("Finish or close the current FreeCAD task first.")
        panel = StaticPanel(selected)
        Gui.Control.showDialog(panel)
        return panel
    except Exception as error:  # noqa: BLE001 - native UI error boundary
        QtWidgets.QMessageBox.information(Gui.getMainWindow(), "Vinkulum", str(error))
        return None
