"""Finite-element study workspace: captured inputs, direct solver, inspectable results."""

import csv
import json
import uuid
from dataclasses import replace
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QAction, QDesktopServices
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QTableView,
    QTextEdit,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from .archive_controller import ArchiveController
from .calculix import StaticStudy, load_static_result, load_study, study_document
from .controls import NumberField
from .finite_elements import POINT_COUNTS
from .model import write_json
from .static_controller import StaticController
from .static_view import StaticTable, StaticViewport
from .theme import apply_theme
from .workspace_link import add_workspace_return


def tension_example():
    """One-cell SI patch example, with lateral supports allowing Poisson strain."""
    return StaticStudy(
        nodes=(
            (0, 0, 0),
            (1, 0, 0),
            (1, 0.1, 0),
            (0, 0.1, 0),
            (0, 0, 0.1),
            (1, 0, 0.1),
            (1, 0.1, 0.1),
            (0, 0.1, 0.1),
        ),
        elements=((1, 2, 3, 4, 5, 6, 7, 8),),
        fixed_dofs=((1, 1), (4, 1), (5, 1), (8, 1), (1, 2), (1, 3), (4, 3)),
        forces=((2, 250, 0, 0), (3, 250, 0, 0), (6, 250, 0, 0), (7, 250, 0, 0)),
        young_pa=210e9,
        poisson=0.3,
    )


class StaticWindow(QMainWindow):
    closed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("Vinkulum Studio · Linear statics / CalculiX")
        self.resize(1380, 880)
        self.setMinimumSize(1040, 720)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.study = None
        self.result = None
        self._loading = False
        self._close_pending = False
        self.controller = StaticController(self)
        self.archive = ArchiveController(self, scheduler=self.controller.scheduler)
        self._build()
        apply_theme(self, "dark")
        self.controller.busy_changed.connect(self._update_busy)
        self.archive.busy_changed.connect(self._update_busy)
        self.archive.stage_changed.connect(self.status.setText)
        self.archive.completed.connect(self._archive_loaded)
        self.archive.problem.connect(self._archive_failed)
        self.archive.cancelled.connect(self.status.setText)
        self.controller.stage_changed.connect(self.status.setText)
        self.controller.problem.connect(self.status.setText)
        self.controller.completed.connect(self._completed)
        self.set_study(tension_example(), "Tension · analytic patch example")

    def _build(self):
        self.setCorner(Qt.Corner.BottomLeftCorner, Qt.DockWidgetArea.LeftDockWidgetArea)
        toolbar = self.addToolBar("Static study")
        add_workspace_return(self, toolbar)
        toolbar.setMovable(False)
        toolbar.addWidget(QLabel("Vinkulum  /  Linear statics"))
        toolbar.addSeparator()
        self.open_action = QAction("Open study…", self)
        self.open_action.triggered.connect(self.open_study)
        toolbar.addAction(self.open_action)
        self.open_result_action = QAction("Open result…", self)
        self.open_result_action.triggered.connect(self.open_result_dialog)
        toolbar.addAction(self.open_result_action)
        self.save_action = QAction("Save study…", self)
        self.save_action.triggered.connect(self.save_study)
        toolbar.addAction(self.save_action)
        self.example_action = QAction("Tension example", self)
        self.example_action.triggered.connect(self.open_example)
        toolbar.addAction(self.example_action)
        toolbar.addSeparator()
        self.mode = QComboBox()
        self.mode.addItems(("Study mesh", "Captured result"))
        self.mode.setAccessibleName("Static display source")
        self.mode.currentIndexChanged.connect(lambda _: self._display(fit=True))
        toolbar.addWidget(self.mode)
        self.viewport = StaticViewport(self)
        self.setCentralWidget(self.viewport)

        panel = QWidget()
        panel.setObjectName("property_fields")
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(6)
        self.title = QLabel()
        self.title.setObjectName("inspector_title")
        self.title.setWordWrap(True)
        layout.addWidget(self.title)
        self.mesh_info = QLabel()
        self.mesh_info.setWordWrap(True)
        self.mesh_info.setObjectName("muted")
        layout.addWidget(self.mesh_info)
        form = QFormLayout()
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapAllRows)
        form.setVerticalSpacing(5)
        self.young = NumberField(210e9)
        self.poisson = NumberField(0.3)
        self.load_factor = NumberField(1)
        self.executable = QLineEdit("ccx")
        self.timeout = NumberField(60)
        self.output = QLineEdit(str(Path.home() / "Vinkulum" / "calculations"))
        self.output.setToolTip(
            "Each run creates a new subdirectory and keeps input, raw output and logs."
        )
        for name, field in (
            ("Young modulus [Pa]", self.young),
            ("Poisson ratio", self.poisson),
            ("Nodal load multiplier", self.load_factor),
        ):
            field.setAccessibleName(name)
            form.addRow(name, field)
        layout.addLayout(form)
        execution = QWidget()
        self.threads = QSpinBox()
        self.threads.setRange(1, self.controller.scheduler.capacity)
        self.threads.setValue(min(4, self.controller.scheduler.capacity))
        advanced = QFormLayout(execution)
        advanced.setContentsMargins(0, 0, 0, 0)
        advanced.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapAllRows)
        for name, field in (
            ("CalculiX executable", self.executable),
            ("CPU threads", self.threads),
            ("Timeout [s]", self.timeout),
            ("Calculation folder", self.output),
        ):
            field.setAccessibleName(name)
            advanced.addRow(name, field)
        execution.hide()
        toggle = QToolButton()
        toggle.setText("Execution settings")
        toggle.setCheckable(True)
        toggle.setArrowType(Qt.ArrowType.RightArrow)
        toggle.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        toggle.toggled.connect(execution.setVisible)
        toggle.toggled.connect(
            lambda checked: toggle.setArrowType(
                Qt.ArrowType.DownArrow if checked else Qt.ArrowType.RightArrow
            )
        )
        layout.addWidget(toggle)
        layout.addWidget(execution)
        for field in (self.young, self.poisson, self.load_factor):
            field.textChanged.connect(self._pending)
        self.pending = QLabel()
        self.pending.setWordWrap(True)
        self.pending.setObjectName("muted")
        layout.addWidget(self.pending)
        buttons = QHBoxLayout()
        self.run_button = QPushButton("Run CalculiX")
        self.run_button.setObjectName("primary")
        self.run_button.clicked.connect(self.start)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self.cancel_work)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.cancel_button)
        layout.addLayout(buttons)
        result_title = QLabel("DISPLAY / CAPTURED RESULT")
        result_title.setObjectName("section")
        layout.addWidget(result_title)
        deformation = QFormLayout()
        self.scale = NumberField(1)
        self.scale.setAccessibleName("Displacement display multiplier")
        self.scale.editingFinished.connect(self._display)
        deformation.addRow("Deformation multiplier", self.scale)
        layout.addLayout(deformation)
        auto_scale = QPushButton("Fit deformation to 10% of extent")
        auto_scale.clicked.connect(self._auto_scale)
        layout.addWidget(auto_scale)
        self.result_summary = QLabel("No captured result.")
        self.result_summary.setWordWrap(True)
        layout.addWidget(self.result_summary)
        self.folder_button = QPushButton("Open captured calculation folder")
        self.folder_button.setEnabled(False)
        self.folder_button.clicked.connect(self._open_folder)
        layout.addWidget(self.folder_button)
        self.limits = QLabel(
            "Linear isotropic elasticity · C3D4 / C3D10 / affine C3D8\nZero supports · nodal loads\nMesh import is independent of the CAD document."
        )
        self.limits.setWordWrap(True)
        self.limits.setObjectName("muted")
        layout.addWidget(self.limits)
        layout.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(panel)
        dock = QDockWidget("Study · SI", self)
        dock.setObjectName("static_study")
        dock.setWidget(scroll)
        dock.setMinimumWidth(290)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)

        output_panel = QWidget()
        out = QVBoxLayout(output_panel)
        out.setContentsMargins(8, 4, 8, 6)
        row = QHBoxLayout()
        self.channels = QComboBox()
        self.channels.addItems(
            (
                "Nodal displacement / reactions",
                "Stress / energy at integration points",
                "Supports / applied nodal forces",
                "Captured provenance",
            )
        )
        self.channels.setAccessibleName("Static result channels")
        self.channels.currentIndexChanged.connect(self._table)
        row.addWidget(self.channels, 1)
        self.export_button = QPushButton("Export table…")
        self.export_button.clicked.connect(self.export_table)
        row.addWidget(self.export_button)
        out.addLayout(row)
        self.table = QTableView()
        self.table.setAccessibleName("Static numeric values with SI units")
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setStretchLastSection(True)
        out.addWidget(self.table)
        self.provenance = QTextEdit()
        self.provenance.setReadOnly(True)
        out.addWidget(self.provenance)
        values = QDockWidget("Values · world frame · unscaled", self)
        values.setObjectName("static_values")
        values.setWidget(output_panel)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, values)
        self.resizeDocks([values], [215], Qt.Orientation.Vertical)
        self.status = QLabel("Ready.")
        self.status.setWordWrap(True)
        self.statusBar().addWidget(self.status, 1)

    def edited_study(self):
        factor = self.load_factor.value()
        if not np.isfinite(factor):
            raise ValueError("Load multiplier must be finite.")
        return replace(
            self.study.scaled_loads(factor),
            young_pa=self.young.value(),
            poisson=self.poisson.value(),
        )

    def _pending(self):
        if self._loading or self.study is None:
            return
        try:
            current = self.edited_study()
            changed = current != self.study
            captured = self.result is not None and current == self.result[0]
            self.pending.setText(
                "These settings are saved in the captured calculation."
                if captured
                else "Edited study settings. Run or save to capture them."
                if changed
                else "Study settings match the loaded document."
            )
            if (
                self.mode.currentIndex() == 1
                and self.result is not None
                and not captured
            ):
                self.pending.setText(
                    "Settings differ from the displayed result. Run to capture a new calculation."
                )
            self._table()
        except (ValueError, TypeError) as error:
            self.pending.setText(str(error))

    def _discard_allowed(self):
        if self.study is None:
            return True
        try:
            current = self.edited_study()
            if current == self.study or (
                self.result is not None and current == self.result[0]
            ):
                return True
        except ValueError, TypeError:
            pass
        return (
            QMessageBox.question(
                self,
                "Unsaved study settings",
                "Discard changes to the study settings? Captured calculations remain on disk.",
                QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel,
            )
            == QMessageBox.StandardButton.Discard
        )

    def set_study(self, study, name="Static study"):
        if self.controller.process is not None:
            raise RuntimeError(
                "Stop the active calculation before replacing its study."
            )
        self._loading = True
        self.study = study
        self.title.setText(name)
        self.limits.setText(
            "Linear isotropic elasticity · C3D4 / C3D10 / affine C3D8\nZero supports · nodal loads\n"
            + (
                "Captured CAD mesh and face conditions included."
                if study.mesh_binding is not None
                else "Mesh import is independent of the CAD document."
            )
        )
        self.mesh_info.setText(
            f"{len(study.nodes):,} nodes · {len(study.elements):,} {study.element_type} elements\n{len(study.fixed_dofs)} constrained DOFs · {len(study.forces)} loaded nodes"
        )
        self.young.setText(repr(study.young_pa))
        self.poisson.setText(repr(study.poisson))
        self.load_factor.setText("1.0")
        self._loading = False
        self.mode.setCurrentIndex(0)
        self._pending()
        self._display(fit=True)

    def open_study(self):
        if self.controller.process is not None or self.archive.busy:
            return
        if not self._discard_allowed():
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open finite-element study", "", "Static study (*.ccx.json *.json)"
        )
        if path:
            self.archive.start(load_study, path, context=("study", Path(path).stem))

    def open_example(self):
        if self._discard_allowed():
            self.set_study(tension_example(), "Tension · analytic patch example")

    def open_result_dialog(self):
        if self.controller.process is not None or self.archive.busy:
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open archived CalculiX result", "", "Captured result (result.json)"
        )
        if not path:
            return
        self.archive.start(load_static_result, path, context=("result", None))

    def _archive_failed(self, error):
        self.status.setText(f"Archive rejected: {error}")

    def _archive_loaded(self, result, context):
        if self._close_pending:
            return
        if context[0] == "study":
            self.set_study(result, context[1])
            return
        self.controller.last_result = result
        self._present_result(result, archived=True)

    def cancel_work(self):
        self.controller.cancel()
        self.archive.cancel()

    def _update_busy(self, *_args):
        self._busy(self.controller.process is not None or self.archive.busy)

    def save_study(self):
        try:
            study = self.edited_study()
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Save finite-element study",
                "study.ccx.json",
                "Static study (*.ccx.json)",
            )
            if path:
                write_json(
                    path,
                    study_document(study),
                )
                self.set_study(study, Path(path).stem)
                self.status.setText(
                    "Study saved. Previous captured results retain their original inputs."
                )
        except (OSError, ValueError, TypeError) as error:
            self.status.setText(str(error))

    def start(self):
        if self.archive.busy:
            return
        try:
            study = self.edited_study()
            root = Path(self.output.text()).expanduser() / (
                "static-" + uuid.uuid4().hex[:16]
            )
            self.controller.start(
                study,
                root,
                executable=self.executable.text(),
                timeout=self.timeout.value(),
                threads=self.threads.value(),
            )
        except (OSError, ValueError, TypeError, RuntimeError) as error:
            self.status.setText(str(error))

    def _busy(self, busy):
        for item in (
            self.open_action,
            self.open_result_action,
            self.save_action,
            self.example_action,
            self.young,
            self.poisson,
            self.load_factor,
            self.executable,
            self.threads,
            self.timeout,
            self.output,
            self.run_button,
        ):
            item.setEnabled(not busy)
        self.cancel_button.setEnabled(busy)
        if not busy and self._close_pending:
            QTimer.singleShot(0, self.close)

    def _completed(self, result):
        if not self._close_pending:
            self._present_result(result)

    def _present_result(self, result, *, archived=False):
        self.result = result
        study, report, root = result
        transport = report.get("input_transport")
        transport_note = (
            f"\nInput conversion: max coordinate change {transport['coordinates_m']['max_absolute_change']:.3g} m; solver geometry rechecked."
            if transport
            else ""
        )
        self.folder_button.setEnabled(True)
        self.result_summary.setText(
            f"CalculiX {report['engine_version']} · captured input\n{len(study.nodes)} nodes · {len(study.elements)} elements\nE = {study.young_pa:.7g} Pa · ν = {study.poisson:.7g}\nElastic energy: {report['strain_energy_J']:.7g} J{transport_note}\nScientific status: NotAssessed\nNo discretisation-error certificate."
        )
        self.mode.setCurrentIndex(1)
        self._display(fit=True)
        self._pending()
        self.status.setText(
            f"Archived result loaded. Study, raw output and stored values checked: {root}"
            if archived
            else f"Completed. Result and diagnostics saved in {root}"
        )

    def _display(self, *args, fit=False):
        if self.study is None:
            return
        try:
            scale = self.scale.value()
            if not np.isfinite(scale) or not 0 <= scale <= 1e12:
                raise ValueError(
                    "Deformation multiplier must be finite and between 0 and 10¹²."
                )
            if self.mode.currentIndex() == 1:
                if self.result is None:
                    self.status.setText("No captured result yet. Run CalculiX first.")
                    self.mode.setCurrentIndex(0)
                    return
                study, report, _ = self.result
                if report.get("schema") in (3, 4):
                    study = study.solver_study
                self.viewport.set_study(study, report, scale=scale, fit=fit)
                self.status.setText(
                    f"Captured displacement |u| in metres · displayed deformation ×{scale:.7g}. Table values are unscaled."
                )
            else:
                self.viewport.set_study(self.study, fit=fit)
            self._table()
        except (ValueError, TypeError) as error:
            self.status.setText(str(error))

    def _auto_scale(self):
        if self.result is None:
            return
        study, report, _ = self.result
        maximum = np.linalg.norm(report["displacements"], axis=1).max()
        factor = (
            min(1e12, 0.1 * np.ptp(study.nodes, axis=0).max() / maximum)
            if maximum > 0
            else 1
        )
        self.scale.setText(repr(float(factor)))
        self.mode.setCurrentIndex(1)
        self._display(fit=True)

    def _table(self):
        if self.study is None:
            return
        channel = self.channels.currentIndex()
        captured = self.result if self.mode.currentIndex() == 1 else None
        if captured:
            study, report, _ = captured
            if report.get("schema") in (3, 4):
                study = study.solver_study
        else:
            try:
                study = self.edited_study()
            except ValueError, TypeError:
                study = self.study
            report = None
        headers, rows = (), ()
        if channel == 0 and report:
            headers = (
                "Node",
                "ux [m]",
                "uy [m]",
                "uz [m]",
                "Rx [N]",
                "Ry [N]",
                "Rz [N]",
            )
            rows = np.column_stack(
                (
                    np.arange(1, len(study.nodes) + 1),
                    report["displacements"],
                    report["reactions"],
                )
            )
        elif channel == 1 and report:
            headers = (
                "Element",
                "IP",
                "sxx [Pa]",
                "syy [Pa]",
                "szz [Pa]",
                "sxy [Pa]",
                "sxz [Pa]",
                "syz [Pa]",
                "ENER [J/m³]",
            )
            rows = np.column_stack(
                (
                    np.repeat(
                        np.arange(1, len(study.elements) + 1),
                        POINT_COUNTS[study.element_type],
                    ),
                    np.tile(
                        np.arange(1, POINT_COUNTS[study.element_type] + 1),
                        len(study.elements),
                    ),
                    report["stress"],
                    report["energy_density"],
                )
            )
        elif channel == 2:
            headers = (
                "Node",
                "Fixed X",
                "Fixed Y",
                "Fixed Z",
                "Fx [N]",
                "Fy [N]",
                "Fz [N]",
            )
            rows = np.zeros((len(study.nodes), 7))
            rows[:, 0] = np.arange(1, len(study.nodes) + 1)
            for node, axis in study.fixed_dofs:
                rows[node - 1, axis] = 1
            for node, x, y, z in study.forces:
                rows[node - 1, 4:] = x, y, z
        self.table_model = StaticTable(headers, rows, self.table)
        old = self.table.model()
        self.table.setModel(self.table_model)
        if old is not None:
            old.deleteLater()
        self.table.resizeColumnsToContents()
        self.table.setVisible(channel != 3)
        self.provenance.setVisible(channel == 3)
        self.export_button.setEnabled(bool(headers))
        metadata = (
            {
                key: value
                for key, value in report.items()
                if key
                not in (
                    "displacements",
                    "external_forces",
                    "reactions",
                    "stress",
                    "energy_density",
                )
            }
            if report
            else {
                "status": "No captured result selected",
                "study": "Use Run CalculiX, then select Captured result.",
            }
        )
        self.provenance.setPlainText(json.dumps(metadata, indent=2))

    def export_table(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Export unscaled values", "static-values.csv", "CSV (*.csv)"
        )
        if path:
            try:
                with Path(path).open("w", encoding="utf-8", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(self.table_model.headers)
                    writer.writerows(self.table_model.rows.tolist())
            except OSError as error:
                self.status.setText(str(error))

    def _open_folder(self):
        if self.result:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.result[2])))

    def closeEvent(self, event):
        if not self._close_pending and not self._discard_allowed():
            event.ignore()
            return
        if self.controller.process is not None or self.archive.busy:
            self._close_pending = True
            self.cancel_work()
            event.ignore()
            return
        self.archive.shutdown()
        self.controller.shutdown()
        self.viewport.shutdown()
        self.closed.emit()
        event.accept()
