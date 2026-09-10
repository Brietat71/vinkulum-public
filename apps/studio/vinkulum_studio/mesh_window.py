"""CAD-to-statics workspace: captured mesh, explicit face conditions, SI study."""

import os
import uuid
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, QThread, QTimer, Signal
from PySide6.QtGui import QAction, QColor, QUndoCommand, QUndoStack
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QToolButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .boundary_display import boundary_frames, effective_values
from .calculix import load_study, study_document
from .controls import NumberField, VectorField, line_icon
from .mesh_binding import (
    MeshBinding,
    MeshCondition,
    mesh_measurements,
    surface_integrals,
)
from .mesh_view import SURFACE_COLORS, MeshViewport
from .meshing import HxtMeshRequest, MeshRequest, load_mesh
from .meshing_controller import MeshingController
from .model import finite_number, write_json
from .theme import apply_theme
from .workspace_link import add_workspace_return

KIND_NAMES = {
    "support": "Support",
    "pressure": "Pressure",
    "total_force": "Total force",
}


class ConditionDialog(QDialog):
    def __init__(self, kind, surfaces, parent=None, existing=None):
        super().__init__(parent)
        self.kind, self.surfaces, self.existing = (
            kind,
            tuple(sorted(surfaces)),
            existing,
        )
        self.condition = None
        self.setWindowTitle(
            ("Edit " if existing else "Add ") + KIND_NAMES[kind].lower()
        )
        self.setMinimumWidth(380)
        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.name = QLineEdit(existing.name if existing else KIND_NAMES[kind])
        self.name.setAccessibleName("Boundary condition name")
        form.addRow("Name", self.name)
        faces = QLabel(", ".join(str(i) for i in self.surfaces))
        faces.setWordWrap(True)
        form.addRow("Selected mesh faces", faces)
        self.axes = [QCheckBox(axis) for axis in ("X", "Y", "Z")]
        if kind == "support":
            row = QHBoxLayout()
            for i, field in enumerate(self.axes, 1):
                field.setChecked(existing is None or i in existing.axes)
                field.setAccessibleName("Fix world " + field.text())
                row.addWidget(field)
            form.addRow("Fixed world directions", row)
            hint = "Zero displacement in the selected directions, at every node of these faces."
        elif kind == "pressure":
            self.pressure = NumberField(existing.values[0] if existing else 1e5)
            self.pressure.setAccessibleName("Pressure in pascals")
            form.addRow("Pressure [Pa]", self.pressure)
            hint = "Positive pressure acts into the solid, along each local surface normal."
        else:
            self.force = VectorField(
                existing.values if existing else (0.0, 0.0, -1000.0), ("X", "Y", "Z")
            )
            self.force.setAccessibleName("Total force in world coordinates, newtons")
            form.addRow("Total force [N]", self.force)
            hint = "World components of one total force, distributed over the union of the selected faces."
        layout.addLayout(form)
        label = QLabel(hint)
        label.setWordWrap(True)
        label.setObjectName("muted")
        layout.addWidget(label)
        self.error = QLabel()
        self.error.setWordWrap(True)
        layout.addWidget(self.error)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def accept(self):
        try:
            self.condition = MeshCondition(
                self.existing.id if self.existing else str(uuid.uuid4()),
                self.name.text().strip(),
                self.kind,
                self.surfaces,
                tuple(i for i, field in enumerate(self.axes, 1) if field.isChecked())
                if self.kind == "support"
                else (),
                ()
                if self.kind == "support"
                else (self.pressure.value(),)
                if self.kind == "pressure"
                else self.force.values(),
            )
        except (ValueError, TypeError) as error:
            self.error.setText(str(error))
            return
        super().accept()


class ConditionEdit(QUndoCommand):
    def __init__(self, window, conditions, label):
        super().__init__(label)
        self.window, self.before, self.after = (
            window,
            window.conditions,
            tuple(conditions),
        )

    def redo(self):
        self.window._set_conditions(self.after)

    def undo(self):
        self.window._set_conditions(self.before)


class StudyTask(QThread):
    """Bounded archive/study checks off the GUI thread; publish only when idle."""

    def __init__(self, operation, parent):
        super().__init__(parent)
        self.operation, self.value, self.error = operation, None, None

    def run(self):
        try:
            self.value = self.operation()
        except (OSError, ValueError, TypeError, RuntimeError, KeyError) as error:
            self.error = str(error)


def surface_areas(solid):
    points = np.array(solid.mesh.nodes)
    return {
        surface.id: sum(
            sum(surface_integrals(tuple(map(tuple, points[np.array(row) - 1])))[0])
            for row in surface.triangles
        )
        for surface in solid.surfaces
    }


def mesh_presentation(solid, root=None):
    """Prepare geometric display data with the mesh's background checks."""
    return (
        solid,
        root,
        surface_areas(solid),
        mesh_measurements(solid),
        boundary_frames(solid),
    )


class MeshWindow(QMainWindow):
    closed = Signal()

    def __init__(self, body, parent=None, *, capture=None, settings=None):
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle("Vinkulum Studio · CAD to linear statics")
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.resize(1440, 900)
        self.setMinimumSize(1100, 740)
        self.capture, self.settings = capture, settings
        self.source, self.solid, self.mesh_root = None, None, None
        self.conditions, self.selected_faces = (), set()
        self._saved_state = None
        self._task, self._static_window = None, None
        self._close_pending = False
        self._face_items = {}
        self.controller = MeshingController(self)
        self.undo = QUndoStack(self)
        self._build()
        apply_theme(self, "dark")
        self.controller.busy_changed.connect(self._update_busy)
        self.controller.stage_changed.connect(self.status.setText)
        self.controller.failed.connect(
            lambda kind, message: self.status.setText(
                f"{kind.replace('_', ' ').capitalize()}: {message}"
            )
        )
        self.controller.completed.connect(self._meshed)
        self.set_source(body)

    @property
    def busy(self):
        return self.controller.busy or self._task is not None

    def _dock(self, title, name, panel, side, width):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(panel)
        dock = QDockWidget(title, self)
        dock.setObjectName(name)
        dock.setWidget(scroll)
        dock.setMinimumWidth(width)
        self.addDockWidget(side, dock)
        return dock

    @staticmethod
    def _section(layout, text):
        label = QLabel(text)
        label.setObjectName("section")
        layout.addWidget(label)

    def _build(self):
        bar = self.addToolBar("CAD study")
        add_workspace_return(self, bar)
        bar.setMovable(False)
        brand = QLabel("Vinkulum  /  CAD → Linear statics")
        brand.setObjectName("brand")
        bar.addWidget(brand)
        bar.addSeparator()
        self.capture_action = QAction("Capture selected solid", self)
        self.capture_action.triggered.connect(self.capture_selected)
        self.capture_action.setEnabled(self.capture is not None)
        bar.addAction(self.capture_action)
        self.open_mesh_action = QAction("Open mesh…", self)
        self.open_mesh_action.triggered.connect(self.open_mesh)
        bar.addAction(self.open_mesh_action)
        self.open_study_action = QAction("Open CAD study…", self)
        self.open_study_action.triggered.connect(self.open_study)
        bar.addAction(self.open_study_action)
        bar.addSeparator()
        self.undo_action = self.undo.createUndoAction(self, "Undo condition change")
        self.undo_action.setShortcut("Ctrl+Z")
        self.undo_action.setIcon(line_icon("undo"))
        self.redo_action = self.undo.createRedoAction(self, "Redo condition change")
        self.redo_action.setShortcut("Ctrl+Shift+Z")
        self.redo_action.setIcon(line_icon("redo"))
        bar.addAction(self.undo_action)
        bar.addAction(self.redo_action)

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        viewbar = QHBoxLayout()
        viewbar.setContentsMargins(10, 4, 8, 4)
        self.view_title = QLabel("Captured CAD solid")
        viewbar.addWidget(self.view_title, 1)
        self.edges = QCheckBox("Mesh edges")
        self.edges.setChecked(True)
        viewbar.addWidget(self.edges)
        self.symbols = QCheckBox("Symbols")
        self.symbols.setChecked(True)
        self.symbols.setAccessibleName("Show boundary direction symbols")
        self.symbols.setToolTip(
            "Sampled directions: open arrows for pressure, solid arrows for total force, bars for constrained world axes. Symbol lengths do not encode force magnitudes."
        )
        viewbar.addWidget(self.symbols)
        views = QComboBox()
        views.setAccessibleName("Mesh camera direction")
        for label, direction in (
            ("Isometric", "iso"),
            ("Front", "front"),
            ("Side", "side"),
            ("Top", "top"),
        ):
            views.addItem(label, direction)
        viewbar.addWidget(views)
        fit = QToolButton()
        fit.setIcon(line_icon("fit"))
        fit.setToolTip("Fit captured solid")
        fit.setAccessibleName("Fit captured solid")
        viewbar.addWidget(fit)
        layout.addLayout(viewbar)
        self.viewport = MeshViewport(self)
        self.viewport.surface_picked.connect(self._picked)
        self.viewport.surface_hovered.connect(self._hovered)
        self.edges.toggled.connect(self.viewport.set_edges)
        self.symbols.toggled.connect(self.viewport.glyph_layer.set_visible)
        fit.clicked.connect(self.viewport.fit_scene)
        views.currentIndexChanged.connect(
            lambda _: self.viewport.camera(views.currentData())
        )
        layout.addWidget(self.viewport, 1)
        legend = QHBoxLayout()
        legend.setSpacing(12)
        legend.setContentsMargins(10, 5, 10, 5)
        for kind, text in (
            ("selected", "Selected"),
            ("support", "Support"),
            ("pressure", "Pressure"),
            ("total_force", "Force"),
            ("multiple", "Combined"),
        ):
            color = "#" + "".join(f"{v:02x}" for v in SURFACE_COLORS[kind])
            label = QLabel(f'<span style="color:{color}">●</span> {text}')
            legend.addWidget(label)
        legend.addStretch(1)
        layout.addLayout(legend)
        self.symbol_hint = QLabel(
            "Sampled directions · open arrows: pressure · solid arrows: force · bars: fixed axes"
        )
        self.symbol_hint.setObjectName("muted")
        self.symbol_hint.setWordWrap(True)
        self.symbol_hint.setContentsMargins(10, 0, 10, 2)
        layout.addWidget(self.symbol_hint)
        self.face_hint = QLabel(
            "  Click a face · Ctrl/Shift+click to toggle · Drag to orbit · Wheel to zoom"
        )
        self.face_hint.setObjectName("muted")
        layout.addWidget(self.face_hint)
        self.setCentralWidget(central)

        setup = QWidget()
        setup.setObjectName("property_fields")
        left = QVBoxLayout(setup)
        left.setContentsMargins(12, 8, 12, 10)
        left.setSpacing(6)
        self.source_title = QLabel()
        self.source_title.setObjectName("inspector_title")
        self.source_title.setWordWrap(True)
        left.addWidget(self.source_title)
        self.source_info = QLabel()
        self.source_info.setWordWrap(True)
        self.source_info.setObjectName("muted")
        left.addWidget(self.source_info)
        self._section(left, "01  VOLUME MESH")
        form = QFormLayout()
        self.size = NumberField(40.0)
        self.size.setAccessibleName("Target mesh size in millimetres")
        self.order = QComboBox()
        self.order.addItem("Quadratic · C3D10", 2)
        self.order.addItem("Linear · C3D4", 1)
        self.order.setAccessibleName("Tetrahedron order")
        form.addRow("Target size [mm]", self.size)
        form.addRow("Elements", self.order)
        left.addLayout(form)
        self.executable = QLineEdit(
            str(self.settings.value("meshing/executable", "gmsh"))
            if self.settings
            else os.environ.get("VINKULUM_GMSH_EXECUTABLE", "gmsh")
        )
        self.executable.setAccessibleName("Gmsh executable built with OCCT 8")
        self.output = QLineEdit(str(Path.home() / "Vinkulum" / "meshes"))
        self.output.setAccessibleName("Mesh capture directory")
        self.timeout = NumberField(120.0)
        self.timeout.setAccessibleName("Mesh timeout in seconds")
        execution = QWidget()
        advanced = QFormLayout(execution)
        advanced.setContentsMargins(0, 0, 0, 0)
        advanced.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapAllRows)
        engine_row = QHBoxLayout()
        engine_row.addWidget(self.executable, 1)
        browse = QPushButton("…")
        browse.setMaximumWidth(30)
        browse.setAccessibleName("Choose Gmsh executable")
        browse.clicked.connect(self.choose_engine)
        engine_row.addWidget(browse)
        advanced.addRow("Gmsh · OCCT 8 minimum", engine_row)
        self.algorithm = QComboBox()
        self.algorithm.addItem("HXT · parallel 3D", "hxt")
        self.algorithm.addItem("Delaunay · legacy 3D", "delaunay")
        self.algorithm.setAccessibleName("Volume meshing algorithm")
        self.threads = QSpinBox()
        self.threads.setRange(1, self.controller.scheduler.capacity)
        self.threads.setValue(min(4, self.controller.scheduler.capacity))
        self.threads.setAccessibleName("Meshing CPU threads")
        advanced.addRow("3D algorithm", self.algorithm)
        advanced.addRow("CPU threads", self.threads)
        advanced.addRow("Timeout [s]", self.timeout)
        advanced.addRow("Mesh capture directory", self.output)
        execution.hide()
        toggle = QToolButton()
        toggle.setText("Execution settings")
        toggle.setCheckable(True)
        toggle.setArrowType(Qt.ArrowType.RightArrow)
        toggle.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        toggle.toggled.connect(execution.setVisible)
        toggle.toggled.connect(
            lambda value: toggle.setArrowType(
                Qt.ArrowType.DownArrow if value else Qt.ArrowType.RightArrow
            )
        )
        left.addWidget(toggle)
        left.addWidget(execution)
        self.execution_panel = execution
        buttons = QHBoxLayout()
        self.generate = QPushButton("Generate mesh")
        self.generate.clicked.connect(self.start)
        self.cancel = QPushButton("Cancel")
        self.cancel.clicked.connect(self.cancel_work)
        buttons.addWidget(self.generate, 1)
        buttons.addWidget(self.cancel)
        left.addLayout(buttons)
        self.mesh_info = QLabel("Generate a mesh to select boundary faces.")
        self.mesh_info.setWordWrap(True)
        self.mesh_info.setObjectName("muted")
        left.addWidget(self.mesh_info)
        self._section(left, "02  MATERIAL · SI")
        material = QFormLayout()
        self.young, self.poisson = NumberField(210e9), NumberField(0.3)
        self.load_factor = NumberField(1.0)
        self.young.setAccessibleName("Young modulus in pascals")
        self.poisson.setAccessibleName("Poisson ratio")
        self.load_factor.setAccessibleName("Boundary load multiplier")
        material.addRow("Young modulus [Pa]", self.young)
        material.addRow("Poisson ratio", self.poisson)
        material.addRow("Load multiplier", self.load_factor)
        left.addLayout(material)
        note = QLabel(
            "Linear isotropic elasticity. Boundary conditions use world directions."
        )
        note.setObjectName("muted")
        note.setWordWrap(True)
        left.addWidget(note)
        self._section(left, "03  STATIC STUDY")
        self.study_info = QLabel("Add supports and loads on the right.")
        self.study_info.setWordWrap(True)
        left.addWidget(self.study_info)
        self.open_calculix = QPushButton("Open in CalculiX workspace")
        self.open_calculix.setObjectName("primary")
        self.open_calculix.clicked.connect(lambda: self.prepare_study("open"))
        self.save = QPushButton("Save CAD study…")
        self.save.clicked.connect(lambda: self.prepare_study("save"))
        left.addWidget(self.open_calculix)
        left.addWidget(self.save)
        limit = QLabel(
            "Local element checks do not certify solution accuracy. Compare mesh refinement before interpreting stress concentrations."
        )
        limit.setObjectName("muted")
        limit.setWordWrap(True)
        left.addWidget(limit)
        left.addStretch(1)
        self._dock(
            "Study setup",
            "cad_mesh_setup",
            setup,
            Qt.DockWidgetArea.LeftDockWidgetArea,
            300,
        )

        panel = QWidget()
        panel.setObjectName("property_fields")
        right = QVBoxLayout(panel)
        right.setContentsMargins(10, 8, 10, 10)
        right.setSpacing(6)
        self.selection_info = QLabel("No mesh faces yet")
        self.selection_info.setObjectName("inspector_title")
        right.addWidget(self.selection_info)
        self.filter = QLineEdit()
        self.filter.setPlaceholderText("Filter faces…")
        self.filter.setAccessibleName("Filter captured mesh faces")
        self.filter.textChanged.connect(self._filter_faces)
        right.addWidget(self.filter)
        self.faces = QTreeWidget()
        self.faces.setHeaderLabels(("Face", "Area [mm²]"))
        self.faces.setRootIsDecorated(False)
        self.faces.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.faces.setAccessibleName("Captured mesh boundary faces")
        self.faces.setMinimumHeight(180)
        self.faces.setMaximumHeight(260)
        self.faces.setColumnWidth(0, 145)
        self.faces.itemSelectionChanged.connect(self._tree_selection)
        right.addWidget(self.faces)
        row = QHBoxLayout()
        self.select_all = QPushButton("Select visible")
        self.select_all.clicked.connect(
            lambda: self.set_selected_faces(
                {i for i, item in self._face_items.items() if not item.isHidden()}
            )
        )
        self.clear_selection = QPushButton("Clear")
        self.clear_selection.clicked.connect(lambda: self.set_selected_faces(()))
        row.addWidget(self.select_all)
        row.addWidget(self.clear_selection)
        right.addLayout(row)
        self._section(right, "BOUNDARY CONDITIONS")
        self.condition_list = QListWidget()
        self.condition_list.setAccessibleName("Applied supports and loads")
        self.condition_list.setMinimumHeight(100)
        self.condition_list.setMaximumHeight(190)
        self.condition_list.currentRowChanged.connect(self._condition_selected)
        self.condition_list.itemDoubleClicked.connect(lambda _: self.edit_condition())
        right.addWidget(self.condition_list)
        row = QHBoxLayout()
        self.add_buttons = {}
        for kind, text in (
            ("support", "Support…"),
            ("pressure", "Pressure…"),
            ("total_force", "Force…"),
        ):
            button = QPushButton(text)
            button.clicked.connect(lambda _, kind=kind: self.add_condition(kind))
            self.add_buttons[kind] = button
            row.addWidget(button)
        right.addLayout(row)
        row = QHBoxLayout()
        self.edit = QPushButton("Edit…")
        self.edit.clicked.connect(self.edit_condition)
        self.remove = QPushButton("Remove")
        self.remove.clicked.connect(self.remove_condition)
        row.addWidget(self.edit)
        row.addWidget(self.remove)
        right.addLayout(row)
        explanation = QLabel(
            "Select faces in the view or list, then add a condition. Edit uses the current face selection. A new mesh starts with no conditions."
        )
        explanation.setWordWrap(True)
        explanation.setObjectName("muted")
        right.addWidget(explanation)
        right.addStretch(1)
        self._dock(
            "Faces & conditions",
            "cad_mesh_faces",
            panel,
            Qt.DockWidgetArea.RightDockWidgetArea,
            305,
        )
        self.status = QLabel("Ready.")
        self.status.setWordWrap(True)
        self.statusBar().addWidget(self.status, 1)
        self.undo.indexChanged.connect(self._update_busy)
        self.load_factor.valueChanged.connect(self._refresh_condition_display)

    def _state(self):
        return (
            self.solid,
            self.conditions,
            self.young.value(),
            self.poisson.value(),
            self.load_factor.value(),
        )

    def _discard_allowed(self):
        try:
            if self._state() == self._saved_state or (
                not self.conditions and self.solid is None
            ):
                return True
        except ValueError, TypeError:
            pass
        return (
            QMessageBox.question(
                self,
                "Unsaved CAD study",
                "Discard the current material or boundary-condition edits? Saved mesh captures and calculations remain on disk.",
                QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel,
            )
            == QMessageBox.StandardButton.Discard
        )

    def set_source(self, body):
        if self.busy:
            raise RuntimeError("Stop the current operation before changing the source.")
        size = min(1e6, max(0.001, max(body.dimensions) * 1000 / 3))
        MeshRequest(body, size)  # Validate before changing any current state.
        self.source, self.solid, self.mesh_root = body, None, None
        self.source_title.setText(body.name)
        self.source_info.setText(
            f"Captured OCCT {body.cad.occt_version} solid\nCAD volume: {body.cad.volume_m3:.7g} m³\nChanges in the model require a new capture."
        )
        self.size.setText(repr(size))
        self.viewport.set_source(body)
        self.view_title.setText("Captured CAD solid")
        self.mesh_info.setText("Generate a mesh to select boundary faces.")
        self._face_items.clear()
        self.faces.clear()
        self.set_selected_faces(())
        self._set_conditions(())
        self.load_factor.setText("1.0")
        self.undo.clear()
        self._saved_state = self._state()
        self._update_busy()

    def capture_selected(self):
        if self.busy or self.capture is None or not self._discard_allowed():
            return
        try:
            body = self.capture()
            if body is not None:
                self.set_source(body)
        except (ValueError, TypeError, RuntimeError) as error:
            self.status.setText(str(error))

    def choose_engine(self):
        path, _ = QFileDialog.getOpenFileName(self, "Choose Gmsh built with OCCT 8")
        if path:
            self.executable.setText(path)

    def start(self):
        if self.busy:
            return
        if (
            self.conditions
            and QMessageBox.question(
                self,
                "Regenerate mesh",
                "A successful new mesh will clear every boundary condition. Existing face numbers will not be reused. Continue?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel,
                QMessageBox.StandardButton.Cancel,
            )
            != QMessageBox.StandardButton.Yes
        ):
            return
        try:
            request_type = (
                HxtMeshRequest if self.algorithm.currentData() == "hxt" else MeshRequest
            )
            request = request_type(
                self.source,
                self.size.value(),
                self.order.currentData(),
                self.threads.value(),
            )
            directory = Path(self.output.text()).expanduser() / str(uuid.uuid4())
            self.controller.start(
                request,
                directory,
                executable=self.executable.text(),
                timeout=self.timeout.value(),
            )
            if self.settings:
                self.settings.setValue("meshing/executable", self.executable.text())
        except (OSError, ValueError, TypeError, RuntimeError) as error:
            self.status.setText(str(error))

    def _meshed(self, result):
        if self._close_pending:
            return
        solid, root = result
        self._background(
            lambda: mesh_presentation(solid, root),
            self._present_mesh,
            "Measuring captured boundary surfaces…",
        )

    def _background(self, operation, callback, label):
        if self.busy:
            raise RuntimeError("An operation is already running.")
        worker = StudyTask(operation, self)
        self._task = worker
        self.status.setText(label)
        self._update_busy()
        worker.finished.connect(lambda: self._task_finished(worker, callback))
        worker.start()

    def _task_finished(self, worker, callback):
        if worker is not self._task:
            return
        self._task = None
        if not worker.isInterruptionRequested() and not self._close_pending:
            if worker.error:
                self.status.setText(worker.error)
            else:
                try:
                    callback(worker.value)
                except (OSError, ValueError, TypeError, RuntimeError) as error:
                    self.status.setText(str(error))
        else:
            self.status.setText("Operation cancelled. Previous study preserved.")
        worker.deleteLater()
        self._update_busy()

    def _present_mesh(self, capture):
        solid, root, areas, measurements, frames = capture
        self.source, self.solid, self.mesh_root = solid.request.body, solid, root
        self.source_title.setText(self.source.name)
        self.source_info.setText(
            f"Captured OCCT {solid.occ_version} solid\nCAD volume: {self.source.cad.volume_m3:.7g} m³\nChanges in the model require a new capture."
        )
        self.size.setText(repr(solid.request.size_mm))
        self.order.setCurrentIndex(self.order.findData(solid.request.order))
        self.algorithm.setCurrentIndex(
            self.algorithm.findData(
                "hxt" if isinstance(solid.request, HxtMeshRequest) else "delaunay"
            )
        )
        self.threads.setValue(
            min(solid.request.threads, self.controller.scheduler.capacity)
        )
        self.viewport.set_solid(solid, frames=frames)
        self._areas = areas
        self.faces.blockSignals(True)
        self.faces.clear()
        self._face_items = {}
        for surface in solid.surfaces:
            item = QTreeWidgetItem((surface.name, f"{areas[surface.id] * 1e6:.6g}"))
            item.setData(0, Qt.ItemDataRole.UserRole, surface.id)
            item.setTextAlignment(
                1, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            )
            item.setToolTip(
                0,
                f"Captured face {surface.id} · {len(surface.triangles)} boundary triangles",
            )
            item.setToolTip(1, f"Quadrature area: {areas[surface.id]:.17g} m²")
            self.faces.addTopLevelItem(item)
            self._face_items[surface.id] = item
        self.faces.blockSignals(False)
        self.filter.clear()
        self.set_selected_faces(())
        self._set_conditions(())
        self.undo.clear()
        self._saved_state = self._state()
        mesh = solid.mesh
        self.view_title.setText(
            f"{mesh.element_type} · {len(solid.surfaces)} captured faces"
        )
        self.mesh_info.setText(
            f"{len(mesh.nodes):,} nodes · {len(mesh.elements):,} elements\nGmsh {solid.gmsh_version} · OCCT {solid.occ_version}\nMesh size: {solid.request.size_mm:.5g} mm\nVolume difference from CAD: {measurements['relative_volume_difference']:.3%}\nVolume agreement does not bound stress error."
        )
        self.status.setText("Mesh captured. Select faces and add supports or loads.")
        self._update_busy()

    def _picked(self, identifier, additive):
        if self.busy:
            return
        selected = set(self.selected_faces) if additive else set()
        if identifier is not None:
            selected.symmetric_difference_update((identifier,))
        self.set_selected_faces(selected)

    def _hovered(self, identifier):
        if identifier is not None and identifier in self._face_items and not self.busy:
            self.face_hint.setText(
                f"  {self._face_items[identifier].text(0)} · area {self._areas[identifier] * 1e6:.6g} mm² · click to select"
            )
        else:
            self.face_hint.setText(
                "  Click a face · Ctrl/Shift+click to toggle · Drag to orbit · Wheel to zoom"
            )

    def set_selected_faces(self, identifiers):
        self.selected_faces = set(identifiers) & set(self._face_items)
        self.faces.blockSignals(True)
        for identifier, item in self._face_items.items():
            item.setSelected(identifier in self.selected_faces)
        self.faces.blockSignals(False)
        self.viewport.set_selection(self.selected_faces)
        self.selection_info.setText(
            f"{len(self.selected_faces)} selected / {len(self._face_items)} faces"
        )
        self._update_busy()

    def _tree_selection(self):
        selected = {
            item.data(0, Qt.ItemDataRole.UserRole)
            for item in self.faces.selectedItems()
        }
        selected |= {i for i in self.selected_faces if self._face_items[i].isHidden()}
        self.set_selected_faces(selected)

    def _filter_faces(self, query):
        self.faces.blockSignals(True)
        for item in self._face_items.values():
            item.setHidden(query.casefold() not in item.text(0).casefold())
        self.faces.blockSignals(False)

    def _set_conditions(self, conditions):
        if self.solid is not None:
            MeshBinding(self.solid, tuple(conditions))
        self.conditions = tuple(conditions)
        current = self.condition_list.currentItem()
        selected_id = current.data(Qt.ItemDataRole.UserRole) if current else None
        self.condition_list.blockSignals(True)
        self.condition_list.clear()
        for condition in self.conditions:
            item = QListWidgetItem()
            item.setData(Qt.ItemDataRole.UserRole, condition.id)
            self.condition_list.addItem(item)
            if condition.id == selected_id:
                self.condition_list.setCurrentItem(item)
        self.condition_list.blockSignals(False)
        self._refresh_condition_display()
        supports = sum(c.kind == "support" for c in self.conditions)
        self.study_info.setText(
            f"{supports} support(s) · {len(self.conditions) - supports} load(s)\nStudy and transported geometry are checked before opening the solver workspace."
        )
        self._update_busy()

    def _refresh_condition_display(self):
        try:
            factor = self.load_factor.value()
        except ValueError, TypeError:
            factor = None
        invalid = self.viewport.set_conditions(self.conditions, factor)
        for index, condition in enumerate(self.conditions):
            values = effective_values(condition, factor)
            detail = (
                "Fix " + " / ".join("XYZ"[a - 1] for a in condition.axes)
                if condition.kind == "support"
                else "Invalid multiplier or load range"
                if values is None
                else f"Pressure {values[0]:.6g} Pa"
                if condition.kind == "pressure"
                else "Total force (" + ", ".join(f"{v:.5g}" for v in values) + ") N"
            )
            item = self.condition_list.item(index)
            item.setText(
                f"{condition.name} · {detail}\n{len(condition.surfaces)} mesh face(s)"
            )
            item.setForeground(
                QColor(
                    *(
                        SURFACE_COLORS[condition.kind]
                        if values is not None
                        else (246, 143, 143)
                    )
                )
            )
            item.setToolTip(
                f"{KIND_NAMES[condition.kind]} · world frame\nFaces: "
                + ", ".join(map(str, condition.surfaces))
                + (
                    "\nZero displacement on the selected world axes."
                    if condition.kind == "support"
                    else f"\nStored values: {condition.values!r}\nLoad multiplier: {factor!r}\nApplied values: {values!r}\nSymbol size is independent of magnitude; arrows are not nodal forces."
                )
            )
        self._multiplier_valid = finite_number(factor)
        self.symbol_hint.setText(
            "Correct the load multiplier or load values to display their directions."
            if invalid or not self._multiplier_valid
            else "Sampled directions · open arrows: pressure · solid arrows: force · bars: fixed axes"
        )
        self._invalid_symbol_loads = invalid
        self._update_busy()

    def _condition_selected(self, index):
        if 0 <= index < len(self.conditions):
            self.set_selected_faces(self.conditions[index].surfaces)
        self._update_busy()

    def add_condition(self, kind):
        if self.busy or not self.selected_faces:
            return
        if len(self.conditions) >= 100:
            self.status.setText("This study supports at most 100 boundary conditions.")
            return
        dialog = ConditionDialog(kind, self.selected_faces, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.undo.push(
                ConditionEdit(
                    self,
                    (*self.conditions, dialog.condition),
                    "Add " + KIND_NAMES[kind].lower(),
                )
            )

    def edit_condition(self):
        index = self.condition_list.currentRow()
        if (
            self.busy
            or not 0 <= index < len(self.conditions)
            or not self.selected_faces
        ):
            return
        condition = self.conditions[index]
        dialog = ConditionDialog(condition.kind, self.selected_faces, self, condition)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            changed = (
                self.conditions[:index]
                + (dialog.condition,)
                + self.conditions[index + 1 :]
            )
            if changed != self.conditions:
                self.undo.push(ConditionEdit(self, changed, "Edit boundary condition"))

    def remove_condition(self):
        index = self.condition_list.currentRow()
        if not self.busy and 0 <= index < len(self.conditions):
            self.undo.push(
                ConditionEdit(
                    self,
                    self.conditions[:index] + self.conditions[index + 1 :],
                    "Remove boundary condition",
                )
            )

    def prepare_study(self, action):
        if self.busy or self.solid is None:
            return
        if action == "open" and self._static_window is not None:
            if self._static_window.controller.process is not None:
                self.status.setText(
                    "Wait for the active CalculiX calculation before replacing its study."
                )
                return
            if not self._static_window._discard_allowed():
                return
        try:
            binding = MeshBinding(self.solid, self.conditions, self.load_factor.value())
            young, poisson = self.young.value(), self.poisson.value()
        except (ValueError, TypeError) as error:
            self.status.setText(str(error))
            return

        def validate():
            study = binding.study(young, poisson)
            _ = study.input_transport
            return study

        self._background(
            validate,
            lambda study: self._study_ready(study, action),
            "Checking supports, material and actual solver geometry…",
        )

    def _study_ready(self, study, action):
        if action == "save":
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Save CAD static study",
                "cad-study.ccx.json",
                "CAD static study (*.ccx.json)",
            )
            if path:
                write_json(path, study_document(study))
                self._saved_state = self._state()
                self.status.setText(f"Captured CAD study saved: {path}")
            return
        from .static_window import StaticWindow

        if self._static_window is None:
            self._static_window = StaticWindow(self)
            self._static_window.closed.connect(
                lambda: setattr(self, "_static_window", None)
            )
        self._static_window.set_study(study, self.source.name + " · captured CAD mesh")
        self._static_window.show()
        self._static_window.raise_()
        self._static_window.activateWindow()
        self.status.setText(
            "Validated study opened in the CalculiX workspace. Save it to preserve the boundary conditions."
        )

    def open_mesh(self):
        if self.busy or not self._discard_allowed():
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open captured mesh result", "", "Captured mesh (result.json)"
        )
        if path:

            def read():
                solid, root = load_mesh(path)
                return mesh_presentation(solid, root)

            self._background(
                read, self._present_mesh, "Rechecking captured CAD mesh and raw files…"
            )

    def open_study(self):
        if self.busy or not self._discard_allowed():
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open CAD static study", "", "CAD static study (*.ccx.json *.json)"
        )
        if path:

            def read():
                study = load_study(path)
                if study.mesh_binding is None:
                    raise ValueError(
                        "This study has no CAD mesh capture. Open it in the CalculiX workspace."
                    )
                return study, mesh_presentation(study.mesh_binding.solid)

            def present(value):
                study, presentation = value
                binding = study.mesh_binding
                self._present_mesh(presentation)
                self.load_factor.setText(repr(binding.load_factor))
                self._set_conditions(binding.conditions)
                self.young.setText(repr(study.young_pa))
                self.poisson.setText(repr(study.poisson))
                self._saved_state = self._state()
                self.status.setText(f"Captured CAD study opened: {path}")

            self._background(read, present, "Checking captured CAD study…")

    def cancel_work(self):
        self.controller.cancel()
        if self._task is not None:
            self._task.requestInterruption()
        self.status.setText("Cancelling; waiting for the current checks to stop…")

    def _update_busy(self, *_):
        busy = self.busy
        self.capture_action.setEnabled(not busy and self.capture is not None)
        for item in (
            self.open_mesh_action,
            self.open_study_action,
            self.generate,
            self.size,
            self.order,
            self.execution_panel,
            self.young,
            self.poisson,
            self.load_factor,
            self.faces,
            self.filter,
            self.condition_list,
        ):
            item.setEnabled(not busy)
        self.cancel.setEnabled(busy)
        self.viewport.selection_enabled = not busy
        self.select_all.setEnabled(not busy and self.solid is not None)
        self.clear_selection.setEnabled(not busy and bool(self.selected_faces))
        for button in self.add_buttons.values():
            button.setEnabled(not busy and bool(self.selected_faces))
        selected = 0 <= self.condition_list.currentRow() < len(self.conditions)
        self.edit.setEnabled(not busy and selected and bool(self.selected_faces))
        self.remove.setEnabled(not busy and selected)
        ready = (
            not busy
            and self.solid is not None
            and any(c.kind == "support" for c in self.conditions)
            and not getattr(self, "_invalid_symbol_loads", ())
            and getattr(self, "_multiplier_valid", True)
        )
        self.save.setEnabled(ready)
        self.open_calculix.setEnabled(ready)
        self.undo_action.setEnabled(not busy and self.undo.canUndo())
        self.redo_action.setEnabled(not busy and self.undo.canRedo())
        if self._close_pending and not busy:
            QTimer.singleShot(0, self.close)

    def closeEvent(self, event):
        if self.busy:
            if not self._close_pending and not self._discard_allowed():
                event.ignore()
                return
            self._close_pending = True
            self.cancel_work()
            event.ignore()
            return
        if not self._close_pending and not self._discard_allowed():
            event.ignore()
            return
        if self._static_window is not None and not self._static_window.close():
            self._close_pending = False
            event.ignore()
            return
        self.controller.shutdown()
        self.viewport.shutdown()
        self.closed.emit()
        event.accept()
