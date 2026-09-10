"""Edit a solid's feature graph and preview a complete transaction before applying."""

from dataclasses import asdict, replace

from PySide6.QtCore import Qt, QTimer, Slot
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .cad_controller import CadController
from .cad_history import PARAMETERS, recipe_for_body
from .controls import NumberField, VectorField
from .editor import euler_matrix, matrix_euler
from .model import finite_number
from .theme import apply_theme
from .viewport import Viewport


class CadHistoryDialog(QDialog):
    def __init__(self, project, identifier, parent=None):
        super().__init__(parent)
        self.original = next((b for b in project.bodies if b.id == identifier), None)
        if self.original is None:
            raise ValueError("Select a body to edit its CAD features.")
        self.project = project
        self.recipe = recipe_for_body(self.original)
        self.preview_body = self.result_body = None
        self.preview_project = self.result_project = None
        self._close_pending = None
        self._feature_id = None
        self._loading = False
        self._preview_signature = None
        self._preview_ok = False
        self.controller = CadController(self, project=project)
        self.controller.completed.connect(self._completed)
        self.controller.failed.connect(self._failed)
        self.controller.busy_changed.connect(self._busy)
        self.setWindowTitle("CAD features · " + self.original.name)
        self.resize(1240, 840)
        self.setMinimumSize(980, 700)
        self._build()
        apply_theme(self, getattr(parent, "theme_name", "dark"))
        self._populate()
        self._display(fit=True)

    def _build(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        bar = QHBoxLayout()
        title = QLabel("CAD  /  " + self.original.name)
        title.setObjectName("section")
        bar.addWidget(title, 1)
        self.show_assembly = QCheckBox("Show assembly")
        self.show_assembly.toggled.connect(self._display)
        bar.addWidget(self.show_assembly)
        self.mode = QComboBox()
        self.mode.addItems(("Original part", "Last preview"))
        self.mode.model().item(1).setEnabled(False)
        self.mode.setAccessibleName("CAD display source")
        self.mode.currentIndexChanged.connect(self._display)
        bar.addWidget(self.mode)
        layout.addLayout(bar)
        panel = QWidget()
        controls = QVBoxLayout(panel)
        controls.setContentsMargins(0, 0, 4, 0)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(("Feature", "Operation"))
        self.tree.setAccessibleName("CAD feature dependencies")
        self.tree.setMinimumHeight(150)
        self.tree.currentItemChanged.connect(self._selected)
        controls.addWidget(self.tree, 1)
        fields = QWidget()
        fields.setObjectName("property_fields")
        self.form = QFormLayout(fields)
        self.form.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setWidget(fields)
        controls.addWidget(scroll, 2)
        self.hint = QLabel()
        self.hint.setWordWrap(True)
        self.hint.setObjectName("muted")
        controls.addWidget(self.hint)
        density_form = QFormLayout()
        if self.original.cad:
            density = self.original.mass / self.original.cad.volume_m3
        else:
            import math

            d = self.original.dimensions
            volume = (
                math.prod(d)
                if self.original.shape == "box"
                else 4 / 3 * math.pi * d[0] ** 3
                if self.original.shape == "sphere"
                else math.pi * d[0] ** 2 * d[1]
            )
            density = self.original.mass / volume
        self.density = NumberField(density)
        self.density.setAccessibleName("CAD homogeneous density in kg per cubic metre")
        self.density.textEdited.connect(self._pending)
        density_form.addRow("Density [kg/m³]", self.density)
        controls.addLayout(density_form)
        mass_note = QLabel(
            "Regeneration recomputes homogeneous mass and inertia from the solid. Attachments keep their world positions."
        )
        mass_note.setObjectName("muted")
        mass_note.setWordWrap(True)
        controls.addWidget(mass_note)
        self.summary = QLabel("No preview computed yet.")
        self.summary.setWordWrap(True)
        controls.addWidget(self.summary)
        self.viewport = Viewport(self)
        self.viewport.set_editable(False)
        self.viewport.setAccessibleName("Original or previewed CAD solid")
        splitter = QSplitter()
        splitter.addWidget(panel)
        splitter.addWidget(self.viewport)
        splitter.setSizes((360, 860))
        splitter.setStretchFactor(1, 1)
        panel.setMinimumWidth(310)
        layout.addWidget(splitter, 1)
        self.status = QLabel(
            "Edit an input feature, then preview the regenerated solid."
        )
        self.status.setWordWrap(True)
        self.controller.stage_changed.connect(self.status.setText)
        layout.addWidget(self.status)
        buttons = QHBoxLayout()
        self.preview_button = QPushButton("Preview")
        self.preview_button.setObjectName("primary")
        self.preview_button.clicked.connect(self.start_preview)
        self.stop_button = QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(lambda: self.controller.cancel())
        self.apply_button = QPushButton("Apply to model")
        self.apply_button.setEnabled(False)
        self.apply_button.clicked.connect(self.accept)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        for button in (
            self.preview_button,
            self.stop_button,
            self.apply_button,
            cancel,
        ):
            button.setAutoDefault(False)
        buttons.addWidget(self.preview_button)
        buttons.addWidget(self.stop_button)
        buttons.addStretch()
        buttons.addWidget(self.apply_button)
        buttons.addWidget(cancel)
        layout.addLayout(buttons)
        self.preview_action = QAction("Preview CAD", self)
        self.preview_action.setShortcut("Ctrl+Return")
        self.preview_action.triggered.connect(self.start_preview)
        self.addAction(self.preview_action)

    def _populate(self):
        self._loading = True
        self.items = {}
        for feature in self.recipe.ordered():
            item = QTreeWidgetItem((feature.name, feature.kind.replace("_", " ")))
            item.setData(0, Qt.ItemDataRole.UserRole, feature.id)
            item.setToolTip(
                0,
                f"Stable identity: {feature.id}\nInput identities: {', '.join(feature.inputs) or 'None'}",
            )
            if feature.id == self.recipe.root:
                font = item.font(0)
                font.setBold(True)
                item.setFont(0, font)
                item.setToolTip(1, "Current output solid")
            self.tree.addTopLevelItem(item)
            self.items[feature.id] = item
        self.tree.resizeColumnToContents(0)
        self._loading = False
        self.tree.setCurrentItem(self.tree.topLevelItem(0))

    def _selected(self, item, previous):
        if self._loading or item is None:
            return
        if self._feature_id is not None:
            try:
                self.recipe = self.edited_recipe()
            except (ValueError, TypeError) as error:
                self.tree.blockSignals(True)
                self.tree.setCurrentItem(previous)
                self.tree.blockSignals(False)
                self.status.setText(str(error))
                return
        self._feature_id = item.data(0, Qt.ItemDataRole.UserRole)
        self._fields()
        self._pending()

    def _fields(self):
        self._loading = True
        while self.form.rowCount():
            self.form.removeRow(0)
        feature = next(f for f in self.recipe.features if f.id == self._feature_id)
        self.name = QLineEdit(feature.name)
        self.name.setAccessibleName("CAD feature name")
        self.form.addRow("Name", self.name)
        self.name.textEdited.connect(self._pending)
        self.dimensions = []
        for label, value in zip(PARAMETERS[feature.kind], feature.dimensions_mm):
            field = NumberField(value)
            field.setAccessibleName(label + " in millimetres")
            field.textEdited.connect(self._pending)
            self.form.addRow(label + " [mm]", field)
            self.dimensions.append(field)
        if feature.kind == "extrude_sketch":
            self.profile_button = QPushButton("Edit sketch…")
            self.profile_button.clicked.connect(self.edit_profile)
            self.form.addRow("XY profile", self.profile_button)
        self.position = self.rotation = None
        if not feature.inputs:
            self.position = VectorField(feature.position_mm, ("X", "Y", "Z"))
            self.position.setAccessibleName(
                "Feature position in design frame, millimetres"
            )
            self.form.addRow("Position [mm]", self.position)
            self._original_angles = tuple(matrix_euler(feature.orientation))
            self._original_rotation = feature.orientation
            self.rotation = VectorField(self._original_angles, ("X", "Y", "Z"))
            self.rotation.setAccessibleName("Feature orientation Rz Ry Rx in degrees")
            self.form.addRow("Rz Ry Rx [°]", self.rotation)
            self.position.textEdited.connect(self._pending)
            self.rotation.textEdited.connect(self._pending)
        if feature.kind == "snapshot":
            message = "Captured solid input. Its source geometry is frozen; placement remains editable."
            if feature.source_body_id:
                message += (
                    " Changing the original tool body does not update this input."
                )
        elif feature.kind == "fillet":
            message = "Round every edge of the input solid. Failure preserves the previous preview and model."
        elif feature.inputs:
            names = {f.id: f.name for f in self.recipe.features}
            message = (
                "Inputs: "
                + ", ".join(names[key] for key in feature.inputs)
                + ". Select an input feature to edit its geometry."
            )
        else:
            message = "Dimensions and placement are in the part's design frame. Moving the centre of mass does not move that frame."
        self.hint.setText(message)
        self._loading = False

    def edited_recipe(self):
        if self._feature_id is None:
            return self.recipe
        feature = next(f for f in self.recipe.features if f.id == self._feature_id)
        changes = {
            "name": self.name.text(),
            "dimensions_mm": tuple(f.value() for f in self.dimensions),
        }
        if self.position is not None:
            changes["position_mm"] = self.position.values()
            angles = self.rotation.values()
            changes["orientation"] = (
                self._original_rotation
                if angles == self._original_angles
                else euler_matrix(angles)
            )
        return self.recipe.replace_feature(replace(feature, **changes))

    def edit_profile(self):
        from .sketch_dialog import SketchDialog

        try:
            recipe = self.edited_recipe()
            feature = next(f for f in recipe.features if f.id == self._feature_id)
            dialog = SketchDialog(feature.profile, self)
            if dialog.exec() == QDialog.DialogCode.Accepted:
                self.recipe = recipe.replace_feature(
                    replace(feature, profile=dialog.profile)
                )
                self._pending()
        except (ValueError, TypeError) as error:
            self.status.setText(str(error))

    def _signature(self):
        recipe = self.edited_recipe()
        density = self.density.value()
        if not finite_number(density) or not 0 < density <= 1e6:
            raise ValueError("Density must be positive and at most 10⁶ kg/m³.")
        return recipe, density

    def _pending(self, *_):
        if self._loading:
            return
        try:
            signature = self._signature()
            if self._feature_id is not None:
                self.items[self._feature_id].setText(0, self.name.text())
        except ValueError, TypeError:
            signature = None
        idle = self.controller.process is None
        self.preview_button.setEnabled(idle and signature is not None)
        self.preview_action.setEnabled(idle and signature is not None)
        self.apply_button.setEnabled(
            idle and self._preview_ok and signature == self._preview_signature
        )
        if (
            self.preview_body is not None
            and signature != self._preview_signature
            and idle
        ):
            self.status.setText(
                "Parameters differ from the last preview. Preview again before applying."
            )

    def start_preview(self):
        if self.controller.process is not None:
            return
        try:
            recipe, density = self._signature()
            self.recipe = recipe
            self._running_signature = (recipe, density)
            self._preview_ok = False
            self.status.setText("Regenerating CAD features…")
            self.controller.start(
                {
                    "operation": "regenerate",
                    "a": asdict(self.original),
                    "recipe": asdict(recipe),
                    "density": density,
                }
            )
        except (OSError, ValueError, TypeError, RuntimeError) as error:
            self.status.setText(str(error))
            self._pending()

    def _busy(self, busy):
        if not busy and self._close_pending is not None:
            QTimer.singleShot(0, self._finish_close)
        self.tree.setEnabled(not busy)
        self.form.parentWidget().setEnabled(not busy)
        self.density.setEnabled(not busy)
        self.stop_button.setEnabled(busy)
        self._pending()

    def _completed(self, result):
        if self._close_pending is not None:
            return
        admission = self.controller.last_admission
        body = admission.body
        self.preview_body = body
        self.preview_project = admission.project
        self._preview_signature = self._running_signature
        self._preview_ok = True
        self.mode.model().item(1).setEnabled(True)
        self.mode.setCurrentIndex(1)
        self._display()
        self.summary.setText(
            f"Preview · {body.cad.volume_m3:.8g} m³\nMass: {body.mass:.8g} kg\n{len(body.cad.recipe.features)} features · OCCT {body.cad.occt_version}"
        )
        self.summary.setToolTip(
            f"Exact stored volume: {body.cad.volume_m3!r} m³\nExact stored mass: {body.mass!r} kg"
        )
        self._pending()
        self.status.setText(
            "Preview ready. Apply to model commits one undoable change."
        )

    def _failed(self, kind, message):
        self._preview_ok = False
        self._pending()
        self.status.setText(message)
        self.status.setToolTip(self.controller.last_log)

    def _display(self, *_args, fit=False):
        body = (
            self.preview_body
            if self.mode.currentIndex() == 1 and self.preview_body is not None
            else self.original
        )
        project = (
            self.preview_project
            if self.mode.currentIndex() == 1 and self.preview_project is not None
            else self.project
        )
        if not self.show_assembly.isChecked():
            project = replace(project, bodies=(body,), joints=(), loads=())
        self.viewport.set_project(project, fit=fit)

    def accept(self):
        self._pending()
        if self.apply_button.isEnabled():
            self.result_body = self.preview_body
            self.result_project = self.preview_project
            super().accept()

    @Slot()
    def _finish_close(self):
        if self._close_pending is not None:
            self.done(self._close_pending)

    def done(self, result):
        if self.controller.process is not None:
            self._close_pending = QDialog.DialogCode.Rejected
            self.controller.cancel()
            return
        self.controller.shutdown()
        self.viewport.close()
        super().done(result)
