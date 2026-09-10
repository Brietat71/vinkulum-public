"""Rigid-mechanism editor: document transactions, 3D authoring and captured runs."""

import json
import math
import time
from bisect import bisect_right
from dataclasses import asdict, replace
from functools import wraps

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QTreeWidgetItem,
    QVBoxLayout,
)

from . import __version__
from .controller import Controller
from .controls import NumberField, VectorField, line_icon
from .dialogs import JointDialog, LawDialog, numbers
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
from .examples3d import EXAMPLES
from .labels import JOINT_LABELS, LAW_LABELS
from .run_archive import RunArchive, series_keys
from .workspace import Workspace


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


class EditorWindow(Workspace, QMainWindow):
    def __init__(self, settings=None):
        super().__init__()
        self.setWindowTitle(f"Vinkulum Studio {__version__} · 3D mechanisms")
        self.resize(1440, 950)
        self.settings = settings
        self.history = History(pendulum())
        self.controller = Controller(self)
        self.result = None
        self.run_archive = RunArchive()
        self.selection = None
        self.fields = {}
        self.original_fields = {}
        self.dirty_fields = False
        self._refreshing = False
        self._path = None
        self._saved = self.history.current
        self._rendered_settings = None
        self._static_window = None
        self._mesh_window = None
        self._articulated_window = None
        self.timer = QTimer(self)
        self.timer.setInterval(33)
        self.timer.timeout.connect(self._tick)
        self._series = []
        self._make_layout()
        self.controller.busy_changed.connect(self._busy)
        self.controller.completed.connect(self._completed)
        self.controller.problem.connect(self._problem)
        self.controller.stage_changed.connect(self.status.setText)
        self.viewport.selected.connect(self.select_object)
        self.viewport.pose_committed.connect(self.move_body)
        self._refresh(fit=True)
        QTimer.singleShot(0, self.viewport.fit_scene)

    @property
    def project(self):
        return self.history.current

    def open_static_study(self):
        from .static_window import StaticWindow

        if self._static_window is None:
            self._static_window = StaticWindow(self)
            self._static_window.closed.connect(
                lambda: setattr(self, "_static_window", None)
            )
        self._static_window.show()
        self._static_window.raise_()
        self._static_window.activateWindow()

    def open_articulated_study(self):
        from .articulated_window import ArticulatedWindow

        if self._articulated_window is None:
            if not self.apply_properties():
                return

            def capture():
                return self.project if self.apply_properties() else None

            self._articulated_window = ArticulatedWindow(
                self.project, self, capture=capture, settings=self.settings
            )
            self._articulated_window.closed.connect(
                lambda: setattr(self, "_articulated_window", None)
            )
        self._articulated_window.show()
        self._articulated_window.raise_()
        self._articulated_window.activateWindow()

    def open_cad_study(self):
        from .mesh_window import MeshWindow

        def capture():
            if not self.apply_properties():
                return None
            body = self.object()
            if not isinstance(body, Body) or body.cad is None:
                raise ValueError("Select one CAD solid to create a static study.")
            return body

        try:
            if self._mesh_window is None:
                body = capture()
                if body is None:
                    return
                self._mesh_window = MeshWindow(
                    body, self, capture=capture, settings=self.settings
                )
                self._mesh_window.closed.connect(
                    lambda: setattr(self, "_mesh_window", None)
                )
            self._mesh_window.show()
            self._mesh_window.raise_()
            self._mesh_window.activateWindow()
        except (ValueError, TypeError, RuntimeError) as error:
            self._problem(str(error))

    def open_cad(self, operation=None):
        if self.mode.currentIndex() != 0 or self.controller.process is not None:
            return
        if not self.apply_properties():
            return
        from importlib.metadata import PackageNotFoundError, version

        try:
            if int(version("cadquery-ocp-novtk").split(".")[0]) < 8:
                raise ValueError("OCCT 8 or later is required.")
            version("build123d")
        except PackageNotFoundError, ValueError:
            self.status.setText(
                "CAD module unavailable: install the Studio CAD environment (OCCT 8 or later)."
            )
            return
        from .cad_dialog import CadDialog
        from .model import atomic_text

        dialog = CadDialog(self.project, self.selection, self)
        if isinstance(operation, str):
            dialog.operation.setCurrentIndex(dialog.operation.findData(operation))
        captured = self.project
        if dialog.exec() != QDialog.DialogCode.Accepted or not dialog.result_data:
            return
        try:
            if "step" in dialog.result_data:
                atomic_text(dialog.output_path, dialog.result_data["step"])
                self.status.setText(
                    "Part exported to STEP with its world pose and millimetre units."
                )
            else:
                if self.project is not captured:
                    raise ValueError(
                        "The document changed during the CAD operation; result not applied."
                    )
                body = dialog.result_body
                if body is None or dialog.result_project is None:
                    raise ValueError("No validated CAD document returned.")
                self._commit(dialog.result_project, fit=True)
                self.select_object(body.id)
                self.status.setText(
                    "CAD solid added · mass and inertia computed from BREP in SI units."
                )
        except (ValueError, OSError) as error:
            self.status.setText(str(error))

    def open_cad_history(self):
        if self.mode.currentIndex() != 0 or self.controller.process is not None:
            return
        if not self.apply_properties():
            return
        from .cad_history_dialog import CadHistoryDialog

        try:
            captured = self.project
            dialog = CadHistoryDialog(captured, self.selection, self)
            if (
                dialog.exec() != QDialog.DialogCode.Accepted
                or dialog.result_body is None
            ):
                return
            if self.project is not captured:
                raise ValueError(
                    "The document changed during CAD editing; preview not applied."
                )
            if dialog.result_project is None:
                raise ValueError("No validated CAD document returned.")
            self._commit(dialog.result_project, fit=True)
            self.select_object(dialog.result_body.id)
            self.status.setText(
                "CAD features regenerated · mass, inertia and attachments updated."
            )
        except (ValueError, TypeError, OSError) as error:
            self.status.setText(str(error))

    @property
    def display_project(self):
        return (
            self.result.project
            if self.result is not None and self.mode.currentIndex() == 1
            else self.project
        )

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
                ("Bodies", shown.bodies),
                ("Joints", shown.joints),
                ("Loads", shown.loads),
            ):
                group = QTreeWidgetItem([label])
                self.tree.addTopLevelItem(group)
                for obj in objects:
                    kind = (
                        {
                            "box": "Box",
                            "sphere": "Sphere",
                            "cylinder": "Cylinder",
                            "cad": "CAD solid",
                        }[obj.shape]
                        if isinstance(obj, Body)
                        else JOINT_LABELS[obj.kind]
                        if isinstance(obj, Joint)
                        else "Load"
                    )
                    item = QTreeWidgetItem([obj.name, kind])
                    item.setToolTip(0, f"{obj.name} · {obj.id}")
                    item.setIcon(
                        0,
                        line_icon(
                            obj.shape
                            if isinstance(obj, Body)
                            else "joint"
                            if isinstance(obj, Joint)
                            else "load",
                            self.palette().text().color().name(),
                        ),
                    )
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
            self._filter_tree()
            self._workspace_state()
        finally:
            self._refreshing = False

    def _text(self, key, label, value):
        text = (
            value
            if isinstance(value, str)
            else ", ".join(format(v, ".12g") for v in value)
            if isinstance(value, tuple)
            else repr(value)
        )
        if isinstance(value, tuple):
            axes = (
                tuple("XYZ")
                if len(value) == 3
                else ("R", "H")
                if len(value) == 2
                else ("R",)
                if len(value) == 1
                else tuple(str(i + 1) for i in range(len(value)))
            )
            field = VectorField(value, axes)
            text = field.text()
        elif isinstance(value, (int, float)):
            field = NumberField(value)
            text = field.text()
        else:
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
            self.apply_button.setEnabled(True)
            self._workspace_state()
            self.status.setText(
                "Properties changed but not applied. Click Apply or run a calculation to validate them."
            )

    def _section(self, title):
        label = QLabel(title)
        label.setObjectName("section")
        self.form.addRow(label)

    def _properties(self):
        while self.form.rowCount():
            self.form.removeRow(0)
        self.fields = {}
        self.original_fields = {}
        self.dirty_fields = False
        obj = self.object(display=True)
        if obj is None:
            self._text("gravity", "Gravity X,Y,Z [m/s²]", self.display_project.gravity)
        else:
            self._text("name", "Name", obj.name)
            self.fields["name"].setToolTip(f"Stable identity: {obj.id}")
            if isinstance(obj, Body):
                self._section("Geometry & mass")
                self._text(
                    "dimensions",
                    {
                        "box": "Dimensions [m]",
                        "sphere": "Radius [m]",
                        "cylinder": "Radius, height [m]",
                        "cad": "BREP extents [m]",
                    }[obj.shape],
                    obj.dimensions,
                )
                if obj.cad:
                    self.fields["dimensions"].setEnabled(False)
                    summary = QLabel(
                        f"Volume: {obj.cad.volume_m3:.6g} m³ · Frame at centre of mass"
                    )
                    summary.setToolTip(
                        f"Exact BREP geometry · OCCT {obj.cad.occt_version}\n"
                        f"build123d {obj.cad.build123d_version}\nVolume exact : {obj.cad.volume_m3!r} m³"
                    )
                    summary.setWordWrap(True)
                    summary.setObjectName("muted")
                    self.form.addRow(summary)
                    edit_features = QPushButton("Edit CAD features…")
                    edit_features.clicked.connect(self.open_cad_history)
                    edit_features.setEnabled(
                        self.mode.currentIndex() == 0
                        and self.controller.process is None
                    )
                    self.form.addRow(edit_features)
                self._text("mass", "Mass [kg]", obj.mass)
                self._section("Transform · world")
                self._text("position", "Position [m]", obj.position)
                self._text(
                    "orientation",
                    "Orientation [°] · Rz Ry Rx",
                    matrix_euler(obj.orientation),
                )
                self._section("Inertia · local frame")
                self._combo(
                    "inertia_mode",
                    "Inertia",
                    [("Homogeneous", "homogeneous"), ("Explicit", "explicit")],
                    obj.inertia_mode,
                )
                self._text(
                    "explicit_inertia",
                    "3×3 tensor, row-major [kg·m²]",
                    obj.explicit_inertia,
                )
                diagonal = np.diag(np.array(obj.inertia()).reshape(3, 3))
                label = QLabel(
                    "Active diagonal: "
                    + ", ".join(format(v, ".6g") for v in diagonal)
                    + " kg·m²"
                )
                label.setWordWrap(True)
                label.setObjectName("muted")
                self.form.addRow(label)
                self.form.setRowVisible(
                    self.fields["explicit_inertia"], obj.inertia_mode == "explicit"
                )
                self.fields["inertia_mode"].currentIndexChanged.connect(
                    lambda _: self.form.setRowVisible(
                        self.fields["explicit_inertia"],
                        self.fields["inertia_mode"].currentData() == "explicit",
                    )
                )
            elif isinstance(obj, Joint):
                self.form.addRow(QLabel(JOINT_LABELS[obj.kind]))
                options = [("Ground", None)] + [
                    (b.name, b.id) for b in self.display_project.bodies
                ]
                for key, label in (("a", "Body A"), ("b", "Body B")):
                    opts = list(options)
                    value = getattr(obj, key)
                    if value is not None and not any(v == value for _, v in opts):
                        opts.append(("Missing body", value))
                    self._combo(key, label, opts, value)
                self._text("pa", "Local anchor A X,Y,Z [m]", obj.pa)
                self._text("pb", "Local anchor B X,Y,Z [m]", obj.pb)
                self._text("ra", "Frame A · angles X,Y,Z [°]", matrix_euler(obj.ra))
                self._text("rb", "Frame B · angles X,Y,Z [°]", matrix_euler(obj.rb))
                if obj.kind in {"pivot", "glissiere"}:
                    button = QPushButton(
                        "Configure prescribed motion"
                        if obj.motion is None
                        else f"Motion: {LAW_LABELS[obj.motion.kind]}"
                    )
                    button.clicked.connect(self.edit_motion)
                    self.form.addRow(button)
                    button = QPushButton("Release motion")
                    button.clicked.connect(lambda: self.set_motion(None))
                    self.form.addRow(button)
            else:
                self._combo(
                    "body",
                    "Bodies",
                    [(b.name, b.id) for b in self.display_project.bodies]
                    + (
                        []
                        if any(b.id == obj.body for b in self.display_project.bodies)
                        else [("Missing body", obj.body)]
                    ),
                    obj.body,
                )
                self._text("point", "Local point X,Y,Z [m]", obj.point)
                for field, unit, label in (
                    ("force", "N", "Force"),
                    ("moment", "N·m", "Moment"),
                ):
                    for index, axis in enumerate("XYZ"):
                        button = QPushButton(
                            f"{label} {axis} [{unit}]: {LAW_LABELS[getattr(obj, field)[index].kind]}"
                        )
                        button.clicked.connect(
                            lambda checked=False, f=field, i=index, u=unit: (
                                self.edit_load_law(f, i, u)
                            )
                        )
                        self.form.addRow(button)
        button = QPushButton("Apply properties")
        button.setObjectName("primary")
        button.clicked.connect(self.apply_properties)
        self.apply_button = button
        self.apply_button.setEnabled(False)
        self.form.addRow(button)

    def _value(self, key, count=None):
        field = self.fields[key]
        if isinstance(field, QComboBox):
            return field.currentData()
        if isinstance(field, VectorField):
            return field.values()
        if isinstance(field, NumberField):
            return field.value()
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
                        self.fields["orientation"].values()
                        != self.original_fields["orientation"][1]
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
                        if self.fields[key].values() != self.original_fields[key][1]:
                            values[key] = euler_matrix(self._value(key, 3))
                else:
                    values.update(
                        body=self._value("body"), point=self._value("point", 3)
                    )
                project = self.project.replace_object(replace(obj, **values))
            applied_from_button = self.apply_button.hasFocus()
            self.history.commit(project)
            self._refresh()
            if applied_from_button:
                # The focused button was destroyed with the old form. Its new
                # replacement is disabled until the next edit; resume at the
                # first property instead of dropping keyboard focus entirely.
                field = next(iter(self.fields.values()))
                if isinstance(field, VectorField):
                    field = field.components[0]
                field.setFocus(Qt.FocusReason.TabFocusReason)
            self.status.setText("Properties applied to the document.")
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
        body = Body(new_id(), f"Body {len(self.project.bodies) + 1}", shape, dimensions)
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
            self.status.setText("Select the body to apply a load to.")
            return
        load = Load(new_id(), f"Load {len(self.project.loads) + 1}", body.id)
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
            duplicate = replace(obj, id=new_id(), name=(obj.name + " copy")[:128])
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
            self.status.setText("Running on a snapshot. The design remains editable.")
            self.controller.start(self.project)
        except (ValueError, RuntimeError, OSError) as exc:
            self._problem(str(exc))

    def stop(self):
        self.stop_button.setEnabled(False)
        self.status.setText("Stopping…")
        self.controller.cancel()

    def _busy(self, busy):
        allowed = not busy and not self.project.diagnostics()
        self.run_button.setEnabled(allowed)
        self.stop_button.setEnabled(busy)
        self.commands["run"].setEnabled(allowed)
        self.commands["stop"].setEnabled(busy)
        self.run_state.setText(
            "Running · captured snapshot" if busy else "Ready to run"
        )

    def _completed(self, result):
        removed = self.run_archive.add(result)
        self._update_run_choices(result.run_id)
        self._display_result(result)
        if removed:
            self.status.setText(
                self.status.text()
                + " The oldest results have left the session history."
            )

    def _display_result(self, result):
        self.pause()
        self.result = result
        self.mode.model().item(1).setEnabled(True)
        self.slider.blockSignals(True)
        self.slider.setRange(0, len(result.time) - 1)
        self.slider.setValue(0)
        self.slider.blockSignals(False)
        self.time_input.blockSignals(True)
        self.time_input.setRange(float(result.time[0]), float(result.time[-1]))
        self.time_input.setSingleStep(float(result.project.step))
        self.time_input.blockSignals(False)
        self.time_input.setEnabled(True)
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
        if self.mode.currentIndex() == 1:
            self.docks["results"].show()
            self.docks["results"].raise_()
        self.status.setText(
            f"Calculation completed · {len(result.time)} samples · trajectory error not assessed."
        )

    def _problem(self, message):
        self.status.setText(
            message + (" The previous result is preserved." if self.result else "")
        )
        self._caption(previous=True)

    def _caption(self, previous=False):
        if self.result:
            p = self.result.project
            state = (
                "Previous result preserved"
                if previous or p != self.project
                else "Computed result"
            )
            self.result_label.setText(
                f"{state} · {self.result.run_id[:8]} · {p.name}, revision {p.revision} · {len(p.bodies)} bodies · duration {p.duration:g} s, step {p.step:g} s"
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
            kind = series_keys(self.result.project)[index][1]
            label += " · world" if kind in {"position", "velocity"} else " · joint"
            self.curve.set_series(self.result.time, values, label, unit)
            self.sample_model.set_series(self.result.time, values, label, unit)
            self.update_comparison()

    def _frame(self, index):
        if self.result is None:
            return
        if self.mode.currentIndex() == 1:
            self.viewport.set_poses(self.result.poses(index), self.result.time[index])
        self.time_input.blockSignals(True)
        self.time_input.setValue(float(self.result.time[index]))
        self.time_input.blockSignals(False)
        selection = self.sample_table.selectionModel()
        selection.blockSignals(True)
        self.sample_table.selectRow(index)
        selection.blockSignals(False)
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
        if self.mode.currentIndex() != 1:
            return
        self._play_origin = (
            time.monotonic()
            - float(self.result.time[self.slider.value()]) / self.speed.currentData()
        )
        self.timer.start()
        self.play_button.setText("Pause")

    def pause(self):
        self.timer.stop()
        self.play_button.setText("Play")

    def _tick(self):
        elapsed = (time.monotonic() - self._play_origin) * self.speed.currentData()
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
        self._workspace_state()
        self.status.setText("Project saved.")
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
            "Modified project",
            "Save changes before continuing?",
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Save,
        )
        if answer == QMessageBox.StandardButton.Save:
            return bool(self.save_dialog())
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
                result = action(path)
                return result is not False
            except (OSError, ValueError, TypeError, KeyError) as exc:
                self.status.setText(str(exc))
        return False

    def save_dialog(self):
        return self._file(True, "Save project", "Project (*.json)", self.save)

    def open_dialog(self):
        if self._discard_allowed():
            self._file(False, "Open project", "Project (*.json)", self.load)

    def import_dialog(self):
        if self._discard_allowed():
            self._file(
                False,
                "Import G0 pendulum",
                "G0 parameters (*.json)",
                lambda path: self.load(path, True),
            )

    def export_dialog(self):
        if self.result:
            self._file(
                True,
                "Export displayed result",
                "CSV (*.csv)",
                self.result.export_csv,
            )

    def provenance(self):
        if self.result is None:
            return
        dialog = QDialog(self)
        dialog.setWindowTitle("Result provenance")
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
        if self._static_window is not None and not self._static_window.close():
            event.ignore()
            return
        if self._mesh_window is not None and not self._mesh_window.close():
            event.ignore()
            return
        if (
            self._articulated_window is not None
            and not self._articulated_window.close()
        ):
            event.ignore()
            return
        self.pause()
        self.save_workspace()
        self.controller.shutdown()
        self.viewport.shutdown()
        event.accept()
