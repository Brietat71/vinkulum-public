"""Rigid-mechanism editor: document transactions, 3D authoring and captured runs."""

from bisect import bisect_right
from dataclasses import asdict, replace
from functools import wraps
import json
import math
import time

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
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
    QSlider,
    QTextEdit,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from . import __version__
from .controller import Controller
from .dialogs import JointDialog, LawDialog, numbers, numeric_text
from .document import (
    Body,
    History,
    Joint,
    Law,
    Load,
    Project,
    import_g0,
    load_project,
    new_id,
    pendulum,
    save_project,
)
from .series_view import SeriesView
from .viewport import Viewport
from .examples3d import EXAMPLES


def report_edit_errors(method):
    @wraps(method)
    def checked(self, *args, **kwargs):
        try:
            return method(self, *args, **kwargs)
        except (ValueError, TypeError) as exc:
            self.status.setText(str(exc))
            return False

    return checked


def euler_matrix(angles):
    x, y, z = np.radians(angles)
    cx, sx, cy, sy, cz, sz = (
        math.cos(x),
        math.sin(x),
        math.cos(y),
        math.sin(y),
        math.cos(z),
        math.sin(z),
    )
    R = (
        np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
        @ np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
        @ np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    )
    return tuple(float(v) for v in R.flat)


def matrix_euler(values):
    R = np.array(values).reshape(3, 3)
    pitch = math.atan2(-R[2, 0], math.hypot(R[0, 0], R[1, 0]))
    if abs(math.cos(pitch)) > 1e-10:
        roll, yaw = math.atan2(R[2, 1], R[2, 2]), math.atan2(R[1, 0], R[0, 0])
    else:
        roll, yaw = 0.0, math.atan2(-R[0, 1], R[1, 1])
    return tuple(math.degrees(v) for v in (roll, pitch, yaw))


class EditorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"Vinkulum Studio {__version__} · Mécanismes 3D")
        self.resize(1440, 950)
        self.history = History(pendulum())
        self.controller = Controller(self)
        self.result = None
        self.selection = None
        self.fields = {}
        self.original_fields = {}
        self.dirty_fields = False
        self._refreshing = False
        self._path = None
        self._saved = self.history.current
        self._rendered_settings = None
        self.timer = QTimer(self)
        self.timer.setInterval(33)
        self.timer.timeout.connect(self._tick)
        self._series = []
        self._make_layout()
        self.controller.busy_changed.connect(self._busy)
        self.controller.completed.connect(self._completed)
        self.controller.problem.connect(self._problem)
        self.viewport.selected.connect(self.select_object)
        self.viewport.pose_committed.connect(self.move_body)
        self._refresh(fit=True)

    @property
    def project(self):
        return self.history.current

    @property
    def display_project(self):
        return (
            self.result.project
            if self.result is not None and self.mode.currentIndex() == 1
            else self.project
        )

    def _make_layout(self):
        self._design_actions = []
        toolbar = self.addToolBar("Conception")
        toolbar.setMovable(False)
        for label, slot, shortcut in (
            ("Nouveau", self.new_project, "Ctrl+N"),
            ("Ouvrir", self.open_dialog, "Ctrl+O"),
            ("Enregistrer", self.save_dialog, "Ctrl+S"),
            ("Importer G0", self.import_dialog, None),
            ("Annuler", self.undo, "Ctrl+Z"),
            ("Rétablir", self.redo, "Ctrl+Shift+Z"),
        ):
            action = QAction(label, self)
            action.triggered.connect(slot)
            if shortcut:
                action.setShortcut(QKeySequence(shortcut))
            toolbar.addAction(action)
            if label in {"Annuler", "Rétablir"}:
                self._design_actions.append(action)
        toolbar.addSeparator()
        examples = QComboBox()
        examples.addItem("Exemples…")
        examples.addItems(EXAMPLES)
        examples.activated.connect(
            lambda index: self.load_example(examples.itemText(index)) if index else None
        )
        toolbar.addWidget(examples)
        for shape, label in (
            ("box", "Boîte"),
            ("cylinder", "Cylindre"),
            ("sphere", "Sphère"),
        ):
            action = QAction(label, self)
            action.triggered.connect(lambda checked=False, s=shape: self.add_body(s))
            action.setShortcut(
                {"box": "Alt+B", "cylinder": "Alt+C", "sphere": "Alt+S"}[shape]
            )
            toolbar.addAction(action)
            self._design_actions.append(action)
        for label, slot in (
            ("Liaison", self.add_joint_dialog),
            ("Charge", self.add_load),
            ("Dupliquer", self.duplicate),
            ("Supprimer", self.delete),
        ):
            action = QAction(label, self)
            action.triggered.connect(slot)
            action.setShortcut(
                {
                    "Liaison": "Alt+L",
                    "Charge": "Alt+F",
                    "Dupliquer": "Ctrl+D",
                    "Supprimer": "Ctrl+Delete",
                }[label]
            )
            toolbar.addAction(action)
            self._design_actions.append(action)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabel("MÉCANISME")
        self.tree.currentItemChanged.connect(
            lambda current, previous: (
                self.select_object(current.data(0, Qt.ItemDataRole.UserRole))
                if current
                else None
            )
        )
        dock = QDockWidget("Objets", self)
        dock.setWidget(self.tree)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)
        dock.setMinimumWidth(210)
        self.properties = QWidget()
        self.form = QFormLayout(self.properties)
        self.form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow
        )
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.properties)
        dock = QDockWidget("Propriétés · unités SI", self)
        dock.setWidget(scroll)
        dock.setMinimumWidth(360)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        central = QWidget()
        layout = QVBoxLayout(central)
        self.setCentralWidget(central)
        top = QHBoxLayout()
        layout.addLayout(top)
        self.mode = QComboBox()
        self.mode.addItems(["Conception", "Résultat du calcul"])
        self.mode.model().item(1).setEnabled(False)
        self.mode.currentIndexChanged.connect(self._mode_changed)
        top.addWidget(self.mode)
        for direction, label in (
            ("iso", "Isométrique"),
            ("front", "Face"),
            ("side", "Côté"),
            ("top", "Dessus"),
        ):
            button = QPushButton(label)
            button.clicked.connect(
                lambda checked=False, d=direction: self.viewport.camera(d)
            )
            top.addWidget(button)
        button = QPushButton("Cadrer sélection")
        button.clicked.connect(lambda: self.viewport.camera("iso", selection=True))
        top.addWidget(button)
        top.addStretch()
        self.viewport = Viewport()
        layout.addWidget(self.viewport, 1)
        self.result_label = QLabel("Aucun résultat calculé.")
        self.result_label.setWordWrap(True)
        layout.addWidget(self.result_label)
        runbar = QHBoxLayout()
        layout.addLayout(runbar)
        self.duration = QLineEdit("2.0")
        self.step = QLineEdit("0.005")
        self.duration.setMaximumWidth(85)
        self.step.setMaximumWidth(90)
        runbar.addWidget(QLabel("Durée [s]"))
        runbar.addWidget(self.duration)
        runbar.addWidget(QLabel("Pas [s]"))
        runbar.addWidget(self.step)
        self.run_button = QPushButton("Calculer la conception")
        self.run_button.clicked.connect(self.run)
        self.run_button.setShortcut("Ctrl+Return")
        runbar.addWidget(self.run_button)
        self.stop_button = QPushButton("Arrêter")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop)
        runbar.addWidget(self.stop_button)
        self.export_button = QPushButton("Exporter CSV")
        self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self.export_dialog)
        runbar.addWidget(self.export_button)
        self.provenance_button = QPushButton("Provenance")
        self.provenance_button.setEnabled(False)
        self.provenance_button.clicked.connect(self.provenance)
        runbar.addWidget(self.provenance_button)
        self.series_combo = QComboBox()
        self.series_combo.currentIndexChanged.connect(self._select_series)
        layout.addWidget(self.series_combo)
        self.curve = SeriesView()
        self.curve.setMaximumHeight(220)
        layout.addWidget(self.curve)
        playback = QHBoxLayout()
        layout.addLayout(playback)
        self.play_button = QPushButton("Lecture")
        self.play_button.clicked.connect(self.toggle_play)
        self.play_button.setEnabled(False)
        playback.addWidget(self.play_button)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self._frame)
        self.slider.sliderPressed.connect(self.pause)
        self.slider.actionTriggered.connect(lambda _: self.pause())
        playback.addWidget(self.slider, 1)
        self.time_label = QLabel("t = —")
        playback.addWidget(self.time_label)
        self.status = QLabel("Prêt. Sélectionnez un corps pour le modifier.")
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.TextFormat.PlainText)
        layout.addWidget(self.status)
        self.diagnostics = QTreeWidget()
        self.diagnostics.setHeaderLabels(["Objet", "Diagnostic"])
        self.diagnostics.itemActivated.connect(
            lambda item, column: self.select_object(
                item.data(0, Qt.ItemDataRole.UserRole)
            )
        )
        self.diagnostics.setMaximumHeight(130)
        layout.addWidget(self.diagnostics)
        note = QLabel(
            "Conception : déplacer ne résout pas les liaisons. Résultats : précision non certifiée. Flèches orange = forces, violettes = moments ; longueur symbolique, directions mondiales."
        )
        note.setWordWrap(True)
        note.setStyleSheet("color:#657487;font-size:11px")
        layout.addWidget(note)

    def object(self, identifier=None, display=False):
        identifier = self.selection if identifier is None else identifier
        project = self.display_project if display else self.project
        return next(
            (
                o
                for o in (*project.bodies, *project.joints, *project.loads)
                if o.id == identifier
            ),
            None,
        )

    def _refresh(self, fit=False):
        self._refreshing = True
        try:
            self.tree.blockSignals(True)
            self.tree.clear()
            shown = self.display_project
            for label, objects in (
                ("Corps", shown.bodies),
                ("Liaisons", shown.joints),
                ("Charges", shown.loads),
            ):
                group = QTreeWidgetItem([label])
                self.tree.addTopLevelItem(group)
                for obj in objects:
                    item = QTreeWidgetItem([obj.name])
                    item.setData(0, Qt.ItemDataRole.UserRole, obj.id)
                    group.addChild(item)
                    if obj.id == self.selection:
                        self.tree.setCurrentItem(item)
                group.setExpanded(True)
            self.tree.blockSignals(False)
            self._properties()
            settings = (self.project.duration, self.project.step)
            if settings != self._rendered_settings:
                self.duration.setText(repr(self.project.duration))
                self.step.setText(repr(self.project.step))
                self._rendered_settings = settings
            self.properties.setEnabled(self.mode.currentIndex() == 0)
            for action in self._design_actions:
                action.setEnabled(self.mode.currentIndex() == 0)
            self.diagnostics.clear()
            for diagnostic in self.project.diagnostics():
                obj = self.object(diagnostic.object_id)
                item = QTreeWidgetItem(
                    [obj.name if obj else self.project.name, diagnostic.message]
                )
                item.setData(0, Qt.ItemDataRole.UserRole, diagnostic.object_id)
                self.diagnostics.addTopLevelItem(item)
            self._show_scene(fit)
            self._busy(self.controller.process is not None)
            self._caption()
        finally:
            self._refreshing = False

    def _text(self, key, label, value):
        text = (
            value
            if isinstance(value, str)
            else numeric_text(value)
            if isinstance(value, tuple)
            else repr(value)
        )
        field = QLineEdit(text)
        field.setAccessibleName(label)
        field.textEdited.connect(self._pending)
        self.form.addRow(label, field)
        self.fields[key] = field
        self.original_fields[key] = (text, value)
        return field

    def _combo(self, key, label, options, current):
        combo = QComboBox()
        for text, data in options:
            combo.addItem(text, data)
        combo.setCurrentIndex(combo.findData(current))
        combo.currentIndexChanged.connect(self._pending)
        self.form.addRow(label, combo)
        self.fields[key] = combo
        return combo

    def _pending(self, *args):
        if not self._refreshing:
            self.dirty_fields = True
            self.status.setText(
                "Champs modifiés, non encore appliqués au document. Cliquez Appliquer ou lancez le calcul pour les valider."
            )

    def _properties(self):
        while self.form.rowCount():
            self.form.removeRow(0)
        self.fields = {}
        self.original_fields = {}
        self.dirty_fields = False
        obj = self.object(display=True)
        if obj is None:
            self._text("gravity", "Gravité X,Y,Z [m/s²]", self.display_project.gravity)
        else:
            self._text("name", "Nom", obj.name)
            self.form.addRow(QLabel(f"Identité : {obj.id[:8]}"))
            if isinstance(obj, Body):
                self.form.addRow(
                    QLabel(
                        {
                            "box": "Boîte",
                            "sphere": "Sphère",
                            "cylinder": "Cylindre · axe local Z",
                        }[obj.shape]
                    )
                )
                self._text(
                    "dimensions",
                    {
                        "box": "Dimensions X,Y,Z [m]",
                        "sphere": "Rayon [m]",
                        "cylinder": "Rayon, hauteur [m]",
                    }[obj.shape],
                    obj.dimensions,
                )
                self._text("mass", "Masse [kg]", obj.mass)
                self._text("position", "Position X,Y,Z [m]", obj.position)
                self._text(
                    "orientation",
                    "Angles X,Y,Z [°] · Rz Ry Rx",
                    matrix_euler(obj.orientation),
                )
                self._combo(
                    "inertia_mode",
                    "Inertie",
                    [("Homogène", "homogeneous"), ("Explicite", "explicit")],
                    obj.inertia_mode,
                )
                self._text(
                    "explicit_inertia",
                    "Tenseur 3×3, par lignes [kg·m²]",
                    obj.explicit_inertia,
                )
                diagonal = np.diag(np.array(obj.inertia()).reshape(3, 3))
                label = QLabel(
                    "Diagonale utilisée : " + numeric_text(tuple(diagonal)) + " kg·m²"
                )
                label.setWordWrap(True)
                self.form.addRow(label)
            elif isinstance(obj, Joint):
                self.form.addRow(QLabel(obj.kind.capitalize()))
                options = [("Bâti", None)] + [
                    (b.name, b.id) for b in self.display_project.bodies
                ]
                for key, label in (("a", "Corps A"), ("b", "Corps B")):
                    opts = list(options)
                    value = getattr(obj, key)
                    if value is not None and not any(v == value for _, v in opts):
                        opts.append(("Corps manquant", value))
                    self._combo(key, label, opts, value)
                self._text("pa", "Ancrage A local X,Y,Z [m]", obj.pa)
                self._text("pb", "Ancrage B local X,Y,Z [m]", obj.pb)
                self._text("ra", "Repère A · angles X,Y,Z [°]", matrix_euler(obj.ra))
                self._text("rb", "Repère B · angles X,Y,Z [°]", matrix_euler(obj.rb))
                if obj.kind in {"pivot", "glissiere"}:
                    button = QPushButton(
                        "Configurer le mouvement imposé"
                        if obj.motion is None
                        else f"Mouvement : {obj.motion.kind}"
                    )
                    button.clicked.connect(self.edit_motion)
                    self.form.addRow(button)
                    button = QPushButton("Libérer le mouvement")
                    button.clicked.connect(lambda: self.set_motion(None))
                    self.form.addRow(button)
            else:
                self._combo(
                    "body",
                    "Corps",
                    [(b.name, b.id) for b in self.display_project.bodies]
                    + (
                        []
                        if any(b.id == obj.body for b in self.display_project.bodies)
                        else [("Corps manquant", obj.body)]
                    ),
                    obj.body,
                )
                self._text("point", "Point local X,Y,Z [m]", obj.point)
                for field, unit, label in (
                    ("force", "N", "Force"),
                    ("moment", "N·m", "Moment"),
                ):
                    for index, axis in enumerate("XYZ"):
                        button = QPushButton(
                            f"{label} {axis} [{unit}] : {getattr(obj, field)[index].kind}"
                        )
                        button.clicked.connect(
                            lambda checked=False, f=field, i=index, u=unit: (
                                self.edit_load_law(f, i, u)
                            )
                        )
                        self.form.addRow(button)
        button = QPushButton("Appliquer les propriétés")
        button.clicked.connect(self.apply_properties)
        self.form.addRow(button)

    def _value(self, key, count=None):
        field = self.fields[key]
        if isinstance(field, QComboBox):
            return field.currentData()
        text = field.text()
        if text == self.original_fields[key][0]:
            return self.original_fields[key][1]
        if key == "name":
            return text.strip()
        values = numbers(text, count)
        return values[0] if count == 1 and key == "mass" else values

    def apply_properties(self):
        if self.mode.currentIndex() == 1:
            return True
        if not self.dirty_fields:
            return True
        try:
            obj = self.object()
            if obj is None:
                project = replace(self.project, gravity=self._value("gravity", 3))
            else:
                values = {"name": self._value("name")}
                if isinstance(obj, Body):
                    values.update(
                        dimensions=self._value("dimensions", len(obj.dimensions)),
                        mass=self._value("mass", 1),
                        position=self._value("position", 3),
                        inertia_mode=self._value("inertia_mode"),
                        explicit_inertia=self._value("explicit_inertia", 9),
                    )
                    if (
                        self.fields["orientation"].text()
                        != self.original_fields["orientation"][0]
                    ):
                        values["orientation"] = euler_matrix(
                            self._value("orientation", 3)
                        )
                elif isinstance(obj, Joint):
                    values.update(
                        a=self._value("a"),
                        b=self._value("b"),
                        pa=self._value("pa", 3),
                        pb=self._value("pb", 3),
                    )
                    for key in ("ra", "rb"):
                        if self.fields[key].text() != self.original_fields[key][0]:
                            values[key] = euler_matrix(self._value(key, 3))
                else:
                    values.update(
                        body=self._value("body"), point=self._value("point", 3)
                    )
                project = self.project.replace_object(replace(obj, **values))
            self.history.commit(project)
            self._refresh()
            self.status.setText("Propriétés appliquées au document.")
            return True
        except (ValueError, TypeError) as exc:
            self.status.setText(str(exc))
            return False

    def select_object(self, identifier):
        if self._refreshing or identifier == self.selection:
            return
        if not self.apply_properties():
            self.tree.blockSignals(True)
            current = None
            for group_index in range(self.tree.topLevelItemCount()):
                group = self.tree.topLevelItem(group_index)
                for index in range(group.childCount()):
                    item = group.child(index)
                    if item.data(0, Qt.ItemDataRole.UserRole) == self.selection:
                        current = item
            self.tree.setCurrentItem(current)
            self.tree.blockSignals(False)
            self.viewport.select(self.selection)
            return
        self.selection = identifier
        self._refresh()

    def _commit(self, project, selection=None, fit=False):
        self.history.commit(project)
        if selection is not None:
            self.selection = selection
        self.mode.blockSignals(True)
        self.mode.setCurrentIndex(0)
        self.mode.blockSignals(False)
        self._refresh(fit)

    def add_body(self, shape):
        if not self.apply_properties():
            return
        dimensions = {"box": (0.2, 0.2, 0.2), "sphere": (0.1,), "cylinder": (0.1, 0.4)}[
            shape
        ]
        body = Body(
            new_id(), f"Corps {len(self.project.bodies) + 1}", shape, dimensions
        )
        try:
            self._commit(self.project.replace_object(body), body.id, True)
        except ValueError as exc:
            self.status.setText(str(exc))

    @report_edit_errors
    def add_joint_dialog(self):
        if not self.apply_properties() or not self.project.bodies:
            return
        dialog = JointDialog(self.project, self.selection, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self._commit(self.project.replace_object(dialog.joint), dialog.joint.id)

    @report_edit_errors
    def add_load(self):
        if not self.apply_properties():
            return
        body = self.object()
        if not isinstance(body, Body):
            self.status.setText("Sélectionnez le corps auquel appliquer une charge.")
            return
        load = Load(new_id(), f"Charge {len(self.project.loads) + 1}", body.id)
        self._commit(self.project.replace_object(load), load.id)

    def edit_motion(self):
        if not self.apply_properties():
            return
        joint = self.object()
        if not isinstance(joint, Joint):
            return
        dialog = LawDialog(
            joint.motion or Law(),
            "rad" if joint.kind == "pivot" else "m",
            self.project.duration,
            self,
        )
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.set_motion(dialog.law)

    @report_edit_errors
    def set_motion(self, law):
        if not self.apply_properties():
            return
        joint = self.object()
        if isinstance(joint, Joint):
            self._commit(self.project.replace_object(replace(joint, motion=law)))

    @report_edit_errors
    def edit_load_law(self, field, index, unit):
        if not self.apply_properties():
            return
        load = self.object()
        if not isinstance(load, Load):
            return
        dialog = LawDialog(
            getattr(load, field)[index], unit, self.project.duration, self
        )
        if dialog.exec() == QDialog.DialogCode.Accepted:
            laws = list(getattr(load, field))
            laws[index] = dialog.law
            self._commit(
                self.project.replace_object(replace(load, **{field: tuple(laws)}))
            )

    def move_body(self, identifier, position, orientation):
        if not self.apply_properties():
            self._show_scene()
            return
        body = self.object(identifier)
        if not isinstance(body, Body):
            self._show_scene()
            return
        try:
            self._commit(
                self.project.replace_object(
                    replace(
                        body, position=tuple(position), orientation=tuple(orientation)
                    )
                ),
                identifier,
            )
        except ValueError as exc:
            self.status.setText(str(exc))
            self._show_scene()

    @report_edit_errors
    def duplicate(self):
        if not self.apply_properties():
            return
        obj = self.object()
        if obj:
            duplicate = replace(obj, id=new_id(), name=(obj.name + " copie")[:128])
            self._commit(self.project.replace_object(duplicate), duplicate.id)

    @report_edit_errors
    def delete(self):
        if self.object():
            project = self.project.remove(self.selection)
            self.selection = None
            self._commit(project)

    def undo(self):
        self.history.undo()
        self._refresh()

    def redo(self):
        self.history.redo()
        self._refresh()

    def run(self):
        duration, step = self.duration.text(), self.step.text()
        if not self.apply_properties():
            return
        try:
            project = replace(
                self.project, duration=numbers(duration, 1)[0], step=numbers(step, 1)[0]
            )
            self.history.commit(project)
            self._refresh()
            self.controller.start(self.project)
            self.status.setText(
                "Calcul en cours sur un instantané. La conception reste éditable."
            )
        except (ValueError, RuntimeError, OSError) as exc:
            self._problem(str(exc))

    def stop(self):
        self.controller.cancel()
        self.stop_button.setEnabled(False)
        self.status.setText("Arrêt demandé…")

    def _busy(self, busy):
        self.run_button.setEnabled(not busy and not self.project.diagnostics())
        self.stop_button.setEnabled(busy)

    def _completed(self, result):
        self.pause()
        self.result = result
        self.mode.model().item(1).setEnabled(True)
        self.slider.blockSignals(True)
        self.slider.setRange(0, len(result.time) - 1)
        self.slider.setValue(0)
        self.slider.blockSignals(False)
        self._series = result.series()
        self.series_combo.clear()
        self.series_combo.addItems([label for label, unit, values in self._series])
        for widget in (
            self.play_button,
            self.slider,
            self.export_button,
            self.provenance_button,
        ):
            widget.setEnabled(True)
        if (
            not self.dirty_fields
            and not self.viewport.dragging
            and result.project.id == self.project.id
        ):
            self.mode.setCurrentIndex(1)
        if self.mode.currentIndex() == 1:
            self._refresh()
        self._frame(0)
        self._caption()
        self.status.setText(
            f"Calcul terminé · {len(result.time)} échantillons · erreur de trajectoire non évaluée."
        )

    def _problem(self, message):
        self.status.setText(
            message + (" Le résultat précédent est conservé." if self.result else "")
        )
        self._caption(previous=True)

    def _caption(self, previous=False):
        if self.result:
            p = self.result.project
            state = (
                "Résultat précédent conservé"
                if previous or p != self.project
                else "Résultat calculé"
            )
            self.result_label.setText(
                f"{state} · {self.result.run_id[:8]} · {p.name}, révision {p.revision} · {len(p.bodies)} corps · durée {p.duration:g} s, pas {p.step:g} s"
            )

    def _show_scene(self, fit=False):
        result_mode = self.mode.currentIndex() == 1 and self.result is not None
        self.viewport.set_project(
            self.result.project if result_mode else self.project, fit=fit
        )
        self.viewport.set_editable(not result_mode)
        self.viewport.select(self.selection)
        if result_mode:
            self.viewport.set_poses(
                self.result.poses(self.slider.value()),
                self.result.time[self.slider.value()],
            )

    def _mode_changed(self, index):
        if self._refreshing:
            return
        if index == 1 and self.dirty_fields:
            self.mode.blockSignals(True)
            self.mode.setCurrentIndex(0)
            self.mode.blockSignals(False)
            if not self.apply_properties():
                return
            self.mode.blockSignals(True)
            self.mode.setCurrentIndex(1)
            self.mode.blockSignals(False)
        self._refresh()

    def _select_series(self, index):
        if self.result is not None and 0 <= index < len(self._series):
            label, unit, values = self._series[index]
            self.curve.set_series(self.result.time, values, label, unit)

    def _frame(self, index):
        if self.result is None:
            return
        if self.mode.currentIndex() == 1:
            self.viewport.set_poses(self.result.poses(index), self.result.time[index])
        self.curve.index = index
        self.curve.update()
        self.time_label.setText(f"t = {self.result.time[index]:.5f} s")

    def toggle_play(self):
        if self.timer.isActive():
            self.pause()
            return
        if self.result is None:
            return
        if self.slider.value() == self.slider.maximum():
            self.slider.setValue(0)
        self.mode.setCurrentIndex(1)
        self._play_origin = time.monotonic() - float(
            self.result.time[self.slider.value()]
        )
        self.timer.start()
        self.play_button.setText("Pause")

    def pause(self):
        self.timer.stop()
        self.play_button.setText("Lecture")

    def _tick(self):
        elapsed = time.monotonic() - self._play_origin
        index = min(
            len(self.result.time) - 1,
            max(0, bisect_right(self.result.time, elapsed) - 1),
        )
        self.slider.setValue(index)
        if elapsed >= self.result.time[-1]:
            self.pause()

    def save(self, path):
        duration, step = self.duration.text(), self.step.text()
        if not self.apply_properties():
            return False
        project = replace(
            self.project, duration=numbers(duration, 1)[0], step=numbers(step, 1)[0]
        )
        save_project(path, project)
        self.history.commit(project)
        self._saved = self.project
        self._path = path
        self.status.setText("Projet enregistré.")
        return True

    def load(self, path, legacy=False):
        project = import_g0(path) if legacy else load_project(path)
        self.pause()
        self.history = History(project)
        self.selection = None
        self._path = None if legacy else path
        self._saved = project
        self.mode.setCurrentIndex(0)
        self._refresh(fit=True)

    def _discard_allowed(self):
        settings_changed = (self.duration.text(), self.step.text()) != (
            repr(self.project.duration),
            repr(self.project.step),
        )
        if (
            self.project == self._saved
            and not self.dirty_fields
            and not settings_changed
        ):
            return True
        answer = QMessageBox.question(
            self,
            "Projet modifié",
            "Abandonner les modifications non enregistrées ?",
            QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel,
        )
        return answer == QMessageBox.StandardButton.Discard

    def new_project(self):
        if self._discard_allowed():
            self.pause()
            self.history = History(Project(new_id()))
            self._saved = self.project
            self._path = None
            self.selection = None
            self.mode.setCurrentIndex(0)
            self._refresh(fit=True)

    def load_example(self, name):
        if name in EXAMPLES and self._discard_allowed():
            self.pause()
            self.history = History(EXAMPLES[name]())
            self._saved = self.project
            self._path = None
            self.selection = None
            self.mode.setCurrentIndex(0)
            self._refresh(fit=True)

    def _file(self, save, caption, pattern, action):
        picker = QFileDialog.getSaveFileName if save else QFileDialog.getOpenFileName
        path, _ = picker(self, caption, str(self._path or ""), pattern)
        if path:
            try:
                action(path)
            except (OSError, ValueError, TypeError, KeyError) as exc:
                self.status.setText(str(exc))

    def save_dialog(self):
        self._file(True, "Enregistrer le projet", "Projet (*.json)", self.save)

    def open_dialog(self):
        if self._discard_allowed():
            self._file(False, "Ouvrir le projet", "Projet (*.json)", self.load)

    def import_dialog(self):
        if self._discard_allowed():
            self._file(
                False,
                "Importer un pendule G0",
                "Paramètres G0 (*.json)",
                lambda path: self.load(path, True),
            )

    def export_dialog(self):
        if self.result:
            self._file(
                True,
                "Exporter le résultat affiché",
                "CSV (*.csv)",
                self.result.export_csv,
            )

    def provenance(self):
        if self.result is None:
            return
        dialog = QDialog(self)
        dialog.setWindowTitle("Provenance du résultat")
        dialog.resize(720, 600)
        layout = QVBoxLayout(dialog)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setPlainText(
            json.dumps(
                {
                    "run_id": self.result.run_id,
                    "project": asdict(self.result.project),
                    "manifest": json.loads(self.result.manifest_json),
                },
                indent=2,
                ensure_ascii=False,
            )
        )
        layout.addWidget(text)
        button = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        button.rejected.connect(dialog.reject)
        layout.addWidget(button)
        dialog.exec()

    def closeEvent(self, event):
        if not self._discard_allowed():
            event.ignore()
            return
        self.pause()
        self.controller.shutdown()
        self.viewport.shutdown()
        event.accept()
