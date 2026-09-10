"""Transactional planar sketch editing with persistent entities and dimensional intent."""

from dataclasses import replace
from uuid import uuid4

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QUndoCommand, QUndoStack
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSplitter,
    QToolBar,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .controls import NumberField
from .sketch import (
    KINDS,
    Sketch,
    SketchConstraint,
    SketchEdge,
    SketchPoint,
    bracket_profile,
    dimension_literal,
    settled,
    solve,
)
from .sketch_view import SketchCanvas
from .theme import apply_theme


class SketchEdit(QUndoCommand):
    def __init__(self, editor, before, after, text, selection):
        super().__init__(text)
        self.editor, self.before, self.after = editor, before, after
        self.old_selection, self.selection = editor.selection, selection

    def redo(self):
        self.editor._show(self.after, self.selection)

    def undo(self):
        self.editor._show(self.before, self.old_selection)


class SketchDialog(QDialog):
    def __init__(self, sketch=None, parent=None):
        super().__init__(parent)
        self.profile = sketch if sketch is not None else bracket_profile()
        try:
            self.profile = settled(self.profile)
        except ValueError:
            # Imported drafts outside the solved-coordinate domain stay repairable.
            pass
        self.selection = ()
        self._pending_fields = self._loading = False
        self._poly_first = None
        self.undo = QUndoStack(self)
        self.setWindowTitle("Sketch · " + self.profile.name)
        self.resize(1220, 840)
        self.setMinimumSize(1040, 700)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        toolbar = QToolBar("Sketch tools")
        toolbar.setMovable(False)
        toolbar.addWidget(QLabel("  SKETCH  /  XY · mm  "))
        self.select_action = toolbar.addAction("Select · S")
        self.line_action = toolbar.addAction("Polyline · L")
        self.close_action = toolbar.addAction("Close profile")
        toolbar.addSeparator()
        self.empty_action = toolbar.addAction("New empty")
        self.bracket_action = toolbar.addAction("Bracket example")
        toolbar.addSeparator()
        self.undo_action = self.undo.createUndoAction(self, "Undo")
        self.redo_action = self.undo.createRedoAction(self, "Redo")
        self.undo_action.setShortcut("Ctrl+Z")
        self.redo_action.setShortcut("Ctrl+Shift+Z")
        toolbar.addAction(self.undo_action)
        toolbar.addAction(self.redo_action)
        layout.addWidget(toolbar)
        left = QWidget()
        controls = QVBoxLayout(left)
        controls.setContentsMargins(0, 0, 4, 0)
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(("Entity", "Definition"))
        self.tree.setAccessibleName("Sketch points, segments and dimensions")
        self.tree.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.tree.itemSelectionChanged.connect(self._tree_selected)
        controls.addWidget(self.tree, 3)
        tools = QHBoxLayout()
        self.constraint_kind = QComboBox()
        for key, name in KINDS.items():
            self.constraint_kind.addItem(name, key)
        self.constraint_kind.setAccessibleName("Constraint type")
        tools.addWidget(self.constraint_kind, 1)
        self.add_constraint = QPushButton("Constrain")
        self.add_constraint.clicked.connect(self.constrain)
        tools.addWidget(self.add_constraint)
        controls.addLayout(tools)
        help_text = QLabel(
            "Select a segment for H/V or a signed X/Y distance; select one point to fix it, or two points for a distance. Ctrl/Shift adds to the selection."
        )
        help_text.setWordWrap(True)
        help_text.setObjectName("muted")
        controls.addWidget(help_text)
        form_widget = QWidget()
        self.form = QFormLayout(form_widget)
        self.form.setContentsMargins(0, 5, 0, 0)
        controls.addWidget(form_widget)
        field_buttons = QHBoxLayout()
        self.update_button = QPushButton("Update entity")
        self.update_button.clicked.connect(self.apply_fields)
        self.revert_button = QPushButton("Revert fields")
        self.revert_button.clicked.connect(self._fields)
        field_buttons.addWidget(self.update_button)
        field_buttons.addWidget(self.revert_button)
        controls.addLayout(field_buttons)
        self.remove_button = QPushButton("Remove selected")
        self.remove_button.clicked.connect(self.remove_selected)
        controls.addWidget(self.remove_button)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        controls.addWidget(self.summary)
        self.canvas = SketchCanvas()
        self.canvas.picked.connect(self._picked)
        self.canvas.add_point.connect(self.append_point)
        self.canvas.moved.connect(
            lambda sketch: self._commit(sketch, "Move sketch point", self.selection)
        )
        self.canvas.diagnostic.connect(self._message)
        self.canvas.dimension_activated.connect(self._activate_dimension)
        self.canvas.tool_changed.connect(self._tool_changed)
        right = QWidget()
        view = QVBoxLayout(right)
        view.setContentsMargins(0, 0, 0, 0)
        viewbar = QHBoxLayout()
        self.snap = QComboBox()
        self.snap.setAccessibleName("Sketch point snap spacing in millimetres")
        for text, value in (
            ("Snap off", 0),
            ("Snap 1 mm", 1),
            ("Snap 5 mm", 5),
            ("Snap 10 mm", 10),
        ):
            self.snap.addItem(text, value)
        self.snap.setCurrentIndex(2)
        self.snap.currentIndexChanged.connect(
            lambda _: setattr(self.canvas, "snap_mm", self.snap.currentData())
        )
        viewbar.addWidget(self.snap)
        viewbar.addStretch(1)
        fit = QPushButton("Fit · F")
        fit.clicked.connect(self.canvas.fit)
        viewbar.addWidget(fit)
        view.addLayout(viewbar)
        view.addWidget(self.canvas, 1)
        hint = QLabel(
            "Drag free points · Middle drag: pan · Wheel: zoom · Double-click a dimension to edit"
        )
        hint.setObjectName("muted")
        hint.setWordWrap(True)
        view.addWidget(hint)
        splitter = QSplitter()
        splitter.addWidget(left)
        splitter.addWidget(right)
        left.setMinimumWidth(330)
        splitter.setSizes((350, 850))
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter, 1)
        self.status = QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        buttons = QHBoxLayout()
        self.accept_button = QPushButton("Use profile")
        self.accept_button.setObjectName("primary")
        self.accept_button.clicked.connect(self.accept)
        cancel = QPushButton("Cancel")
        cancel.clicked.connect(self.reject)
        buttons.addStretch()
        buttons.addWidget(self.accept_button)
        buttons.addWidget(cancel)
        layout.addLayout(buttons)
        for button in self.findChildren(QPushButton):
            button.setAutoDefault(False)
        self.select_action.triggered.connect(lambda: self.canvas.set_tool("select"))
        self.line_action.triggered.connect(lambda: self.canvas.set_tool("polyline"))
        self.close_action.triggered.connect(self.close_profile)
        self.empty_action.triggered.connect(
            lambda: self.preset(Sketch(str(uuid4()), "Profile"))
        )
        self.bracket_action.triggered.connect(lambda: self.preset(bracket_profile()))
        apply_theme(self, getattr(parent, "theme_name", "dark"))
        self._show(self.profile, ())
        self._tool_changed("select")

    def _message(self, text):
        self.status.setText(text)

    def _tool_changed(self, tool):
        self._poly_first = None
        self.close_action.setEnabled(tool == "polyline")
        self._message(
            "Click each vertex, then click the first point or choose Close profile. Escape returns to selection."
            if tool == "polyline"
            else "Select an entity or a dimension. Decimal dimensions are in millimetres."
        )

    def preset(self, sketch):
        if not self._flush():
            return
        self._commit(sketch, "Replace sketch profile", ())
        self.canvas.set_tool("select")
        self.canvas.fit()

    def _commit(self, sketch, text, selection, *, drag=None):
        try:
            sketch = settled(sketch, drag=drag)
            if sketch != self.profile:
                self.undo.push(
                    SketchEdit(self, self.profile, sketch, text, tuple(selection))
                )
            else:
                self._show(sketch, selection)
            return True
        except (ValueError, TypeError) as error:
            self._message(str(error))
            return False

    def _show(self, sketch, selection):
        self.profile = sketch
        self._loading = True
        self.tree.clear()
        self.items = {}
        point_names = {p.id: p.name for p in sketch.points}
        try:
            result = solve(sketch)
            problem = None
        except ValueError as error:
            result, problem = None, str(error)
        self._solution = result
        for title, key, rows in (
            ("Dimensions", "constraint", sketch.constraints),
            ("Segments", "edge", sketch.edges),
            ("Points", "point", sketch.points),
        ):
            group = QTreeWidgetItem((title, ""))
            group.setFlags(group.flags() & ~Qt.ItemFlag.ItemIsSelectable)
            self.tree.addTopLevelItem(group)
            for row in rows:
                detail = (
                    KINDS[row.kind]
                    + (
                        " · " + ", ".join(row.values_mm) + " mm"
                        if row.values_mm
                        else ""
                    )
                    if key == "constraint"
                    else f"{point_names[row.start]} → {point_names[row.end]}"
                    if key == "edge"
                    else ", ".join(f"{v:.6g}" for v in row.xy_mm) + " mm"
                )
                item = QTreeWidgetItem((row.name, detail))
                item.setData(0, Qt.ItemDataRole.UserRole, (key, row.id))
                item.setToolTip(0, row.id)
                item.setToolTip(1, detail)
                if result and row.id in result.conflicts:
                    item.setText(1, "Conflicting · " + detail)
                    item.setForeground(1, QColor("#ef9292"))
                group.addChild(item)
                self.items[(key, row.id)] = item
            group.setExpanded(True)
        self.selection = tuple(
            tuple(key) for key in selection if tuple(key) in self.items
        )
        self._sync_selection()
        self.tree.resizeColumnToContents(0)
        self._loading = False
        self.canvas.conflicts = result.conflicts if result else ()
        self.canvas.set_sketch(sketch, self.selection)
        self._fields()
        valid = False
        if result and result.conflicts:
            names = [c.name for c in sketch.constraints if c.id in result.conflicts]
            self.summary.setText("Conflicting dimensions · " + ", ".join(names))
            self._message(
                "Edit or remove a highlighted dimension. The conflicting cycle cannot be satisfied exactly."
            )
        elif result:
            dof = result.degrees_of_freedom
            text = (
                "Empty sketch"
                if not sketch.points
                else "Fully dimensioned"
                if dof == 0
                else "Underconstrained"
            )
            self.summary.setText(
                f"{text} · {dof} remaining degrees of freedom\n{result.redundant_equations} redundant equations"
            )
            try:
                sketch.polygon()
                valid = True
                self._message(
                    "Closed profile ready for extrusion. Free coordinates retain their current positions."
                    if dof
                    else "Closed, fully dimensioned profile ready for extrusion."
                )
            except ValueError as error:
                self._message(str(error))
        else:
            self.summary.setText("Dimension evaluation unavailable")
            self._message(problem)
        self.accept_button.setEnabled(valid)

    def _sync_selection(self):
        self.tree.blockSignals(True)
        self.tree.clearSelection()
        for key in self.selection:
            self.items[key].setSelected(True)
        self.tree.blockSignals(False)

    def select(self, selection):
        if not self._flush():
            self._sync_selection()
            return False
        self.selection = tuple(
            tuple(key) for key in selection if tuple(key) in self.items
        )
        self._sync_selection()
        self.canvas.set_sketch(self.profile, self.selection)
        self._fields()
        return True

    def _tree_selected(self):
        if not self._loading:
            self.select(
                tuple(
                    item.data(0, Qt.ItemDataRole.UserRole)
                    for item in self.tree.selectedItems()
                )
            )

    def _picked(self, key, additive):
        selection = list(self.selection) if additive else []
        if key in selection:
            selection.remove(key)
        elif key:
            selection.append(key)
        self.select(selection)

    def _fields(self):
        self._loading = True
        while self.form.rowCount():
            self.form.removeRow(0)
        self.name_field = None
        self.value_fields = []
        if len(self.selection) == 1:
            kind, key = self.selection[0]
            rows = getattr(
                self.profile,
                {"point": "points", "edge": "edges", "constraint": "constraints"}[kind],
            )
            row = next(r for r in rows if r.id == key)
            self.name_field = QLineEdit(row.name)
            self.name_field.setAccessibleName("Sketch entity name")
            self.name_field.textEdited.connect(self._pending)
            self.form.addRow("Name", self.name_field)
            if kind == "point":
                result = self._solution
                for axis, value in enumerate(row.xy_mm):
                    field = NumberField(value)
                    field.setAccessibleName(f"Point {'XY'[axis]} in millimetres")
                    field.valueChanged.connect(self._pending)
                    field.setEnabled(
                        bool(result and result.free_axes.get(key, (False, False))[axis])
                    )
                    self.value_fields.append(field)
                    self.form.addRow(f"{'XY'[axis]} [mm]", field)
            elif kind == "constraint":
                for axis, value in enumerate(row.values_mm):
                    field = QLineEdit(value)
                    field.setAccessibleName(
                        f"{row.name} decimal dimension {axis + 1} in millimetres"
                    )
                    field.textEdited.connect(self._pending)
                    self.value_fields.append(field)
                    self.form.addRow(
                        f"{'XY'[axis] if row.kind == 'fixed' else 'Value'} [mm]", field
                    )
                note = (
                    "Signed coordinate of the second point minus the first."
                    if row.kind.startswith("distance_")
                    else KINDS[row.kind]
                )
                self.form.addRow(QLabel(note))
        self._loading = False
        self._pending_fields = False
        self.update_button.setEnabled(False)
        self.revert_button.setEnabled(False)
        self.remove_button.setEnabled(bool(self.selection))

    def _pending(self, *_):
        if not self._loading:
            self._pending_fields = True
            self.update_button.setEnabled(True)
            self.revert_button.setEnabled(True)

    def _flush(self):
        return self.apply_fields() if self._pending_fields else True

    def apply_fields(self):
        if len(self.selection) != 1 or self.name_field is None:
            return True
        kind, key = self.selection[0]
        attribute = {"point": "points", "edge": "edges", "constraint": "constraints"}[
            kind
        ]
        rows = getattr(self.profile, attribute)
        row = next(r for r in rows if r.id == key)
        try:
            changes = {"name": self.name_field.text()}
            if kind == "point":
                changes["xy_mm"] = tuple(field.value() for field in self.value_fields)
            elif kind == "constraint":
                changes["values_mm"] = tuple(
                    field.text().strip() for field in self.value_fields
                )
            updated = replace(row, **changes)
            candidate = replace(
                self.profile,
                **{attribute: tuple(updated if r.id == key else r for r in rows)},
            )
            drag = (key, updated.xy_mm) if kind == "point" else None
            return self._commit(
                candidate, "Edit sketch entity", self.selection, drag=drag
            )
        except (ValueError, TypeError) as error:
            self._message(str(error))
            return False

    def _activate_dimension(self, key):
        if self.select((("constraint", key),)) and self.value_fields:
            self.value_fields[0].setFocus()
            self.value_fields[0].selectAll()

    def constrain(self):
        if not self._flush():
            return
        kind = self.constraint_kind.currentData()
        keys = [key for entity, key in self.selection if entity == "point"]
        edges = [key for entity, key in self.selection if entity == "edge"]
        try:
            values = ()
            if kind in ("horizontal", "vertical"):
                if len(edges) != 1 or len(self.selection) != 1:
                    raise ValueError(
                        "Select one segment for a horizontal or vertical constraint."
                    )
                entities = tuple(edges)
            else:
                if (
                    kind.startswith("distance_")
                    and len(edges) == 1
                    and len(self.selection) == 1
                ):
                    edge = next(e for e in self.profile.edges if e.id == edges[0])
                    keys = [edge.start, edge.end]
                required = 1 if kind == "fixed" else 2
                if len(keys) != required:
                    raise ValueError(
                        "Select one point to fix, or two points / one segment for a distance."
                    )
                entities = tuple(keys)
                solution = solve(self.profile)
                if solution.conflicts:
                    raise ValueError(
                        "Resolve conflicting dimensions before adding another constraint."
                    )
                points = solution.exact_coordinates
                values = (
                    tuple(dimension_literal(v) for v in points[keys[0]])
                    if kind == "fixed"
                    else (
                        dimension_literal(
                            points[keys[1]][int(kind == "distance_y")]
                            - points[keys[0]][int(kind == "distance_y")]
                        ),
                    )
                )
            row = SketchConstraint(str(uuid4()), KINDS[kind], kind, entities, values)
            self._commit(
                replace(self.profile, constraints=(*self.profile.constraints, row)),
                "Add sketch constraint",
                (("constraint", row.id),),
            )
            self._activate_dimension(row.id)
        except (ValueError, TypeError) as error:
            self._message(str(error))

    def remove_selected(self):
        if not self._flush():
            return
        selected = set(self.selection)
        try:
            candidate = replace(
                self.profile,
                points=tuple(
                    p for p in self.profile.points if ("point", p.id) not in selected
                ),
                edges=tuple(
                    e for e in self.profile.edges if ("edge", e.id) not in selected
                ),
                constraints=tuple(
                    c
                    for c in self.profile.constraints
                    if ("constraint", c.id) not in selected
                ),
            )
            self._commit(candidate, "Remove sketch entities", ())
        except ValueError as error:
            self._message(
                str(error)
                + " Select dependent constraints and segments too, or remove them first."
            )

    def append_point(self, point, existing):
        if not self._flush():
            return
        try:
            last = self.canvas.last_point
            if existing == last and last is not None:
                return
            points, edges = self.profile.points, self.profile.edges
            if existing is None:
                row = SketchPoint(
                    str(uuid4()), f"Point {len(points) + 1}", tuple(point)
                )
                points, existing = (*points, row), row.id
            if last:
                edges = (
                    *edges,
                    SketchEdge(
                        str(uuid4()), f"Segment {len(edges) + 1}", last, existing
                    ),
                )
            if self._commit(
                replace(self.profile, points=points, edges=edges),
                "Draw sketch segment",
                (("point", existing),),
            ):
                if self._poly_first == existing and last:
                    self.canvas.set_tool("select")
                else:
                    self._poly_first = self._poly_first or existing
                    self.canvas.last_point = existing
                    self.canvas.update()
        except (ValueError, TypeError) as error:
            self._message(str(error))

    def close_profile(self):
        if self._poly_first and self.canvas.last_point != self._poly_first:
            self.append_point(None, self._poly_first)

    def accept(self):
        if not self._flush():
            return
        try:
            self.profile.polygon()
        except ValueError as error:
            self._message(str(error))
            return
        super().accept()
