"""Captured-state editing and operator inspection for the optional Pinocchio worker."""

import csv
import json
import os
import uuid
from dataclasses import asdict
from pathlib import Path

import numpy as np
from PySide6.QtCore import QAbstractTableModel, QSettings, Qt, Signal
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QComboBox,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QTableView,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from .articulated import state_vector, tree_links
from .articulated_result import load_operators
from .controls import NumberField
from .document import MAX_PROJECT_BYTES, Project
from .examples3d import double_pendulum
from .model import finite_number, read_json, write_json
from .pinocchio_controller import PinocchioController
from .theme import apply_theme
from .viewport import Viewport


class OperatorTable(QAbstractTableModel):
    def __init__(self, headers, labels, values, parent=None):
        super().__init__(parent)
        self.headers, self.labels = tuple(headers), tuple(labels)
        self.values = np.array(values, dtype=float, copy=True)
        self.values.setflags(write=False)

    def rowCount(self, parent=None):
        return 0 if parent is not None and parent.isValid() else len(self.labels)

    def columnCount(self, parent=None):
        return 0 if parent is not None and parent.isValid() else len(self.headers)

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        if role == Qt.ItemDataRole.TextAlignmentRole and index.column():
            return Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        if role not in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.ToolTipRole):
            return None
        if index.column() == 0:
            return self.labels[index.row()]
        value = float(self.values[index.row(), index.column() - 1])
        return repr(value) if role == Qt.ItemDataRole.ToolTipRole else f"{value:.7g}"

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole:
            if orientation == Qt.Orientation.Horizontal:
                return (
                    self.headers[section] if 0 <= section < len(self.headers) else None
                )
            return str(section + 1)
        return None


class ArticulatedWindow(QMainWindow):
    closed = Signal()

    def __init__(self, project, parent=None, *, capture=None, settings=None):
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("Vinkulum Studio · Articulated operators / Pinocchio")
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.resize(1440, 940)
        self.setMinimumSize(1060, 740)
        self.settings = (
            settings if settings is not None else QSettings("Vinkulum", "Studio")
        )
        self.capture = capture
        self.project = None
        self.result = None
        self.links = ()
        self._loading = False
        self._saved = None
        self.controller = PinocchioController(self)
        self._build()
        apply_theme(self, "dark")
        self.controller.busy_changed.connect(self._busy)
        self.controller.stage_changed.connect(self.status.setText)
        self.controller.problem.connect(self.status.setText)
        self.controller.completed.connect(self._completed)
        self.set_project(project)

    def _build(self):
        toolbar = self.addToolBar("Articulated analysis")
        toolbar.setMovable(False)
        toolbar.addWidget(QLabel("Vinkulum  /  Articulated operators"))
        toolbar.addSeparator()
        self.capture_action = toolbar.addAction(
            "Capture current model", self.capture_model
        )
        self.capture_action.setEnabled(self.capture is not None)
        self.open_action = toolbar.addAction("Open input…", self.open_input)
        self.save_action = toolbar.addAction("Save input…", self.save_input)
        self.open_result_action = toolbar.addAction("Open result…", self.open_result)
        self.example_action = toolbar.addAction("Double pendulum", self.open_example)
        toolbar.addSeparator()
        self.mode = QComboBox()
        self.mode.addItems(("Authored model", "Captured state"))
        self.mode.setAccessibleName("Articulated display source")
        self.mode.currentIndexChanged.connect(self._display)
        toolbar.addWidget(self.mode)
        self.viewport = Viewport(self)
        self.viewport.setAccessibleName("Articulated mechanism at the captured state")
        self.viewport.set_editable(False)
        self.viewport.selected.connect(self._select_body)
        self.tabs = QTabWidget()
        self.input_table = QTableWidget(0, 5)
        self.input_table.setHorizontalHeaderLabels(
            ("Joint · units", "q", "v", "a", "Actuator effort")
        )
        self.input_table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Stretch
        )
        self.input_table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents
        )
        self.input_table.setAccessibleName("Articulated state inputs")
        self.input_table.setToolTip(
            "q: position [rad or m]; v: velocity [rad/s or m/s]; "
            "a: requested acceleration [rad/s² or m/s²]. "
            "Actuator effort is torque [N m] or force [N] along the joint coordinate."
        )
        self.tabs.addTab(self.input_table, "State inputs")
        result_widget = QWidget()
        layout = QVBoxLayout(result_widget)
        layout.setContentsMargins(6, 6, 6, 6)
        row = QHBoxLayout()
        self.channel = QComboBox()
        self.channel.addItems(
            (
                "Joint response",
                "Mass matrix",
                "Intrinsic dτ/dq",
                "Intrinsic dτ/dv",
                "Body kinematics",
                "Body Jacobian",
            )
        )
        self.channel.currentIndexChanged.connect(self._table)
        self.body_choice = QComboBox()
        self.body_choice.currentIndexChanged.connect(self._table)
        self.export_button = QPushButton("Export table…")
        self.export_button.clicked.connect(self.export_table)
        self.export_button.setEnabled(False)
        row.addWidget(self.channel, 1)
        row.addWidget(self.body_choice)
        row.addWidget(self.export_button)
        layout.addLayout(row)
        self.table = QTableView()
        self.table.setAlternatingRowColors(True)
        self.table.setAccessibleName("Captured operator values")
        layout.addWidget(self.table)
        self.table_note = QLabel("Evaluate a state to inspect captured values.")
        self.table_note.setWordWrap(True)
        layout.addWidget(self.table_note)
        self.tabs.addTab(result_widget, "Captured values")
        self.provenance = QPlainTextEdit()
        self.provenance.setReadOnly(True)
        self.tabs.addTab(self.provenance, "Provenance")
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(self.viewport)
        splitter.addWidget(self.tabs)
        splitter.setSizes((540, 300))
        self.setCentralWidget(splitter)

        panel = QWidget()
        controls = QVBoxLayout(panel)
        self.model_summary = QLabel()
        self.model_summary.setWordWrap(True)
        controls.addWidget(self.model_summary)
        self.admission = QLabel()
        self.admission.setWordWrap(True)
        controls.addWidget(self.admission)
        form = QFormLayout()
        self.load_time = NumberField(0.0)
        self.load_time.textEdited.connect(self._pending)
        form.addRow("Load evaluation time [s]", self.load_time)
        controls.addLayout(form)
        self.execution_toggle = QToolButton()
        self.execution_toggle.setText("Execution settings")
        self.execution_toggle.setCheckable(True)
        self.execution_toggle.setToolButtonStyle(
            Qt.ToolButtonStyle.ToolButtonTextBesideIcon
        )
        self.execution_toggle.setArrowType(Qt.ArrowType.RightArrow)
        controls.addWidget(self.execution_toggle)
        self.execution = QWidget()
        execution_form = QFormLayout(self.execution)
        execution_form.setContentsMargins(0, 0, 0, 0)
        self.interpreter = QLineEdit(
            str(
                self.settings.value(
                    "pinocchio/interpreter",
                    os.environ.get("VINKULUM_PINOCCHIO_PYTHON", ""),
                )
            )
        )
        self.interpreter.setPlaceholderText("Pinocchio environment / bin / python")
        self.interpreter.setAccessibleName("Pinocchio environment Python executable")
        choose = QPushButton("…")
        choose.setFixedWidth(28)
        choose.setAccessibleName("Choose Pinocchio Python executable")
        choose.clicked.connect(self.choose_interpreter)
        python_row = QHBoxLayout()
        python_row.addWidget(self.interpreter)
        python_row.addWidget(choose)
        execution_form.addRow("Python", python_row)
        self.timeout = NumberField(30.0)
        execution_form.addRow("Timeout [s]", self.timeout)
        self.output = QLineEdit(
            str(
                self.settings.value(
                    "pinocchio/output",
                    str(Path.home() / "Vinkulum runs" / "articulated"),
                )
            )
        )
        self.output.setAccessibleName("Articulated calculation output folder")
        execution_form.addRow("Save calculations in", self.output)
        controls.addWidget(self.execution)
        self.execution.setVisible(False)
        self.execution_toggle.toggled.connect(self._toggle_execution)
        self.interpreter.textChanged.connect(self._pending)
        self.run_action = QAction("Evaluate state", self)
        self.run_action.setShortcut("Ctrl+Return")
        self.run_action.triggered.connect(self.start)
        self.addAction(self.run_action)
        buttons = QHBoxLayout()
        self.run_button = QPushButton("Evaluate state")
        self.run_button.clicked.connect(self.run_action.trigger)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(lambda: self.controller.cancel())
        self.cancel_button.setEnabled(False)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.cancel_button)
        controls.addLayout(buttons)
        self.pending = QLabel()
        self.pending.setWordWrap(True)
        controls.addWidget(self.pending)
        heading = QLabel("CAPTURED STATE")
        heading.setObjectName("section")
        controls.addWidget(heading)
        self.result_summary = QLabel("No captured result yet.")
        self.result_summary.setWordWrap(True)
        controls.addWidget(self.result_summary)
        scope = QLabel(
            "Single-state rigid-tree operators.\nRevolute and prismatic joints.\nWorld forces and moments.\nNo trajectory integration.\n\nIntrinsic derivatives exclude derivatives of applied world loads."
        )
        scope.setObjectName("muted")
        scope.setWordWrap(True)
        controls.addWidget(scope)
        controls.addStretch()
        dock = QDockWidget("Analysis · SI", self)
        dock.setWidget(panel)
        dock.setMinimumWidth(300)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)
        self.status = QLabel("Ready.")
        self.status.setWordWrap(True)
        self.statusBar().addWidget(self.status, 1)

    def _toggle_execution(self, checked):
        self.execution.setVisible(checked)
        self.execution_toggle.setArrowType(
            Qt.ArrowType.DownArrow if checked else Qt.ArrowType.RightArrow
        )

    def choose_interpreter(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose Python in the Pinocchio environment", self.interpreter.text()
        )
        if path:
            self.interpreter.setText(path)

    def set_project(self, project, state=None):
        self._loading = True
        self.project = project
        self.model_summary.setText(
            f"{project.name}\n{len(project.bodies)} bodies · {len(project.joints)} joints"
        )
        try:
            self.links = tree_links(project)
            values = {
                key: state_vector(
                    (state or {}).get(key), self.links, key, initial=key == "q"
                )
                for key in ("q", "velocity", "acceleration", "effort")
            }
            self.admission.setText(
                f"{len(self.links)} articulated coordinates\nDeclared A→B signs retained."
            )
        except (ValueError, TypeError) as error:
            self.links = ()
            values = {}
            self.admission.setText(str(error))
        self.input_table.setRowCount(len(self.links))
        self.state_fields = {
            key: [] for key in ("q", "velocity", "acceleration", "effort")
        }
        for i, link in enumerate(self.links):
            item = QTableWidgetItem(
                f"{link.joint.name}\n{link.coordinate_unit} · {link.effort_unit}"
            )
            item.setToolTip(
                f"{link.joint.id}\nPositive: B relative to A along joint A local +Z."
            )
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.input_table.setItem(i, 0, item)
            self.input_table.setRowHeight(i, 44)
            for column, key in enumerate(self.state_fields, 1):
                field = NumberField(float(values[key][i]))
                field.setAccessibleName(f"{link.joint.name} {key}")
                field.textEdited.connect(self._pending)
                self.input_table.setCellWidget(i, column, field)
                self.state_fields[key].append(field)
        self.load_time.setText(str((state or {}).get("time_s", 0.0)))
        self._loading = False
        self._saved = self._signature()
        self.mode.setCurrentIndex(0)
        self._display()
        self._pending()

    def edited_state(self):
        if not self.links:
            raise ValueError(
                "The captured model is outside the supported tree contract."
            )
        result = {
            key: state_vector(
                [field.value() for field in fields], self.links, key
            ).tolist()
            for key, fields in self.state_fields.items()
        }
        result["time_s"] = self.load_time.value()
        if (
            not finite_number(result["time_s"])
            or not 0 <= result["time_s"] <= self.project.duration
        ):
            raise ValueError(
                "Load evaluation time must lie within the project duration."
            )
        return result

    def _signature(self):
        try:
            state = self.edited_state()
        except ValueError, TypeError:
            # Preserve invalid edits in the dirty-state comparison too. Invalid
            # fields must not silently become equivalent to an unsupported model.
            return self.project, (
                "invalid",
                self.load_time.text(),
                tuple(
                    tuple(field.text() for field in fields)
                    for fields in self.state_fields.values()
                ),
            )
        return self.project, state

    def _pending(self, *_):
        if self._loading:
            return
        signature = self._signature()
        captured = (
            None if self.result is None else (self.result.project, self.result.state)
        )
        self.pending.setText(
            "Inputs differ from the displayed result."
            if self.result is not None and signature != captured
            else "Each run preserves its model and state."
        )
        ready = (
            bool(self.links)
            and isinstance(signature[1], dict)
            and bool(self.interpreter.text().strip())
            and self.controller.process is None
        )
        self.run_action.setEnabled(ready)
        self.run_button.setEnabled(ready)
        if not self.interpreter.text().strip():
            self.pending.setText(
                self.pending.text()
                + "\nSelect the separate Pinocchio 4.1 environment in Execution settings."
            )

    def _discard_allowed(self):
        signature = self._signature()
        if signature == self._saved or (
            self.result is not None
            and signature == (self.result.project, self.result.state)
        ):
            return True
        return (
            QMessageBox.question(
                self,
                "Unrecorded analysis inputs",
                "Discard the edited analysis inputs?",
                QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel,
            )
            == QMessageBox.StandardButton.Discard
        )

    def capture_model(self):
        if self.capture is not None and self._discard_allowed():
            project = self.capture()
            if project is not None:
                self.set_project(project)

    def open_example(self):
        if self._discard_allowed():
            self.set_project(double_pendulum())

    def open_input(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open articulated analysis input", "", "Articulated input (*.json)"
        )
        if not path:
            return
        try:
            data = read_json(path, MAX_PROJECT_BYTES)
            if (
                not isinstance(data, dict)
                or set(data) != {"format", "schema", "project", "state"}
                or data["format"] != "vinkulum-articulated-input"
                or type(data["schema"]) is not int
                or data["schema"] != 1
            ):
                raise ValueError("Unknown articulated input format.")
            project = Project.from_dict(data["project"])
            links = tree_links(project)
            state = data["state"]
            if not isinstance(state, dict) or set(state) != {
                "q",
                "velocity",
                "acceleration",
                "effort",
                "time_s",
            }:
                raise ValueError("Incomplete articulated state.")
            for key in ("q", "velocity", "acceleration", "effort"):
                state_vector(state[key], links, key)
            if (
                not finite_number(state["time_s"])
                or not 0 <= state["time_s"] <= project.duration
            ):
                raise ValueError("Invalid load time.")
            if self._discard_allowed():
                self.set_project(project, state)
        except (OSError, ValueError, TypeError, KeyError) as error:
            self.status.setText(str(error))

    def save_input(self):
        try:
            state = self.edited_state()
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Save articulated analysis input",
                "articulated-input.json",
                "Articulated input (*.json)",
            )
            if path:
                write_json(
                    path,
                    {
                        "format": "vinkulum-articulated-input",
                        "schema": 1,
                        "project": asdict(self.project),
                        "state": state,
                    },
                )
                self._saved = (self.project, state)
                self.status.setText(f"Analysis inputs saved in {path}")
        except (OSError, ValueError, TypeError) as error:
            self.status.setText(str(error))

    def open_result(self):
        if self.controller.process is not None:
            return
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open captured Pinocchio operators",
            "",
            "Captured result (result.json)",
        )
        if path:
            try:
                result = load_operators(path)
            except (OSError, ValueError, TypeError, KeyError) as error:
                self.status.setText(f"Result rejected: {error}")
                return
            self.controller.last_result = result
            self._completed(result)
            self.status.setText(
                "Captured operators loaded and checked; current analysis inputs preserved."
            )

    def start(self):
        try:
            state = self.edited_state()
            directory = (
                Path(self.output.text()).expanduser()
                / f"operators-{uuid.uuid4().hex[:16]}"
            )
            self.controller.start(
                self.project,
                state,
                directory,
                interpreter=self.interpreter.text().strip(),
                timeout=self.timeout.value(),
            )
        except (OSError, ValueError, TypeError, RuntimeError) as error:
            self.status.setText(str(error))

    def _busy(self, busy):
        for item in (
            self.capture_action,
            self.open_action,
            self.save_action,
            self.open_result_action,
            self.example_action,
            self.input_table,
            self.load_time,
            self.execution,
        ):
            item.setEnabled(not busy)
        self.capture_action.setEnabled(not busy and self.capture is not None)
        self.cancel_button.setEnabled(busy)
        self._pending()

    def _completed(self, result):
        self.result = result
        self.settings.setValue("pinocchio/interpreter", self.interpreter.text())
        self.settings.setValue("pinocchio/output", self.output.text())
        r = result.report
        self.result_summary.setText(
            f"{result.project.name}\nPinocchio {r['engine_version']} · {r['nv']} coordinates\nLoad time: {result.state['time_s']:.7g} s\nKinetic energy: {r['kinetic_energy_J']:.7g} J\nPotential energy: {r['potential_energy_J']:.7g} J\nScientific status: NotAssessed"
        )
        self.body_choice.blockSignals(True)
        self.body_choice.clear()
        for body in result.project.bodies:
            self.body_choice.addItem(body.name, body.id)
        self.body_choice.blockSignals(False)
        self.mode.setCurrentIndex(1)
        self._display()
        self._table()
        self.tabs.setCurrentIndex(1)
        self._pending()
        self.status.setText(
            f"Captured operators and inputs saved in {result.directory}"
        )

    def _display(self, *_):
        if self.mode.currentIndex() == 1 and self.result is not None:
            self.viewport.set_project(self.result.display_project, fit=True)
        elif self.project is not None:
            self.viewport.set_project(self.project, fit=True)

    def _select_body(self, identifier):
        index = self.body_choice.findData(identifier)
        if index >= 0:
            self.body_choice.setCurrentIndex(index)

    def _table(self, *_):
        if self.result is None:
            return
        r = self.result.report
        order = r["coordinate_order"]
        names = [
            f"{i + 1}. {entry['name']} [{entry['coordinate_unit']}; {entry['effort_unit']}]"
            for i, entry in enumerate(order)
        ]
        channel = self.channel.currentIndex()
        self.body_choice.setVisible(channel == 5)
        note = "Rows retain the captured joint identities. q/v/a use each row's rad or m coordinate; efforts use N m or N."
        if channel == 0:
            headers = (
                "Joint · effort unit",
                "q [rad|m]",
                "v [rad/s|m/s]",
                "a requested [rad/s²|m/s²]",
                "Actuator [N m|N]",
                "External [N m|N]",
                "Bias [N m|N]",
                "Inverse [N m|N]",
                "Acceleration [rad/s²|m/s²]",
            )
            values = np.column_stack(
                [
                    *[
                        r["state"][key]
                        for key in ("q", "velocity", "acceleration", "effort")
                    ],
                    *[
                        r[key]
                        for key in (
                            "external_effort",
                            "intrinsic_bias",
                            "inverse_effort",
                            "forward_acceleration",
                        )
                    ],
                ]
            )
        elif channel in (1, 2, 3):
            variable = {1: "acceleration", 2: "q", 3: "velocity"}[channel]
            derivative = {
                1: r["mass_matrix"],
                2: r["intrinsic_inverse_derivatives"]["q"],
                3: r["intrinsic_inverse_derivatives"]["velocity"],
            }[channel]
            suffix = {1: "/s²", 2: "", 3: "/s"}[channel]
            headers = (
                "Joint · effort unit",
                *[
                    f"{i + 1} [{entry['coordinate_unit']}{suffix}]"
                    for i, entry in enumerate(order)
                ],
            )
            values = np.array(derivative)
            note = f"Each entry maps {variable} in its column to effort in its row. Intrinsic RNEA includes gravity and excludes derivatives of applied world loads."
        elif channel == 4:
            names = [body["name"] for body in r["bodies"]]
            headers = (
                "Body",
                "x [m]",
                "y [m]",
                "z [m]",
                "vx [m/s]",
                "vy [m/s]",
                "vz [m/s]",
                "wx [rad/s]",
                "wy [rad/s]",
                "wz [rad/s]",
            )
            values = np.array(
                [body["position_m"] + body["velocity_world"] for body in r["bodies"]]
            )
            note = "Captured centres of mass and velocities, expressed in world axes."
        else:
            body = next(
                (
                    body
                    for body in r["bodies"]
                    if body["id"] == self.body_choice.currentData()
                ),
                r["bodies"][0],
            )
            names = (
                "vx [m/s]",
                "vy [m/s]",
                "vz [m/s]",
                "wx [rad/s]",
                "wy [rad/s]",
                "wz [rad/s]",
            )
            headers = (
                "Body velocity",
                *[
                    f"{i + 1} [{entry['coordinate_unit']}/s]"
                    for i, entry in enumerate(order)
                ],
            )
            values = np.array(body["jacobian"])
            note = f"{body['name']}: J @ joint velocity gives velocity at this body's centre of mass, in world-aligned axes."
        old = self.table.model()
        model = OperatorTable(headers, names, values, self.table)
        self.table.setModel(model)
        if old is not None:
            old.deleteLater()
        self.table.resizeColumnsToContents()
        self.table_note.setText(note)
        self.export_button.setEnabled(True)
        self.provenance.setPlainText(
            json.dumps(
                {
                    key: value
                    for key, value in r.items()
                    if key
                    not in ("mass_matrix", "intrinsic_inverse_derivatives", "bodies")
                },
                indent=2,
            )
        )

    def export_table(self):
        model = self.table.model()
        if model is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export captured operator values", "operators.csv", "CSV (*.csv)"
        )
        if path:
            try:
                with Path(path).open("w", encoding="utf-8", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(model.headers)
                    for name, row in zip(model.labels, model.values):
                        writer.writerow((name, *map(float, row)))
                self.status.setText(f"Captured values exported to {path}")
            except OSError as error:
                self.status.setText(str(error))

    def closeEvent(self, event):
        if not self._discard_allowed():
            event.ignore()
            return
        self.controller.shutdown()
        self.viewport.close()
        event.accept()
        self.closed.emit()
