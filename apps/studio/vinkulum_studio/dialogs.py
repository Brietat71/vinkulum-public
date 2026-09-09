"""Typed authoring dialogs. Inputs are numbers and tables, never executable code."""

import numpy as np
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)

from .document import Law, joint_at
from .labels import JOINT_LABELS, LAW_LABELS
from .model import finite_number
from .series_view import SeriesView


def numbers(text, count=None):
    try:
        values = tuple(float(v.strip()) for v in text.replace(";", ",").split(","))
    except ValueError as exc:
        raise ValueError("Enter numbers separated by commas.") from exc
    if (
        not all(finite_number(v) for v in values)
        or count is not None
        and len(values) != count
    ):
        raise ValueError(f"{count or 'Expected'} finite values required.")
    return values


def numeric_text(values):
    return ", ".join(repr(float(v)) for v in values)


class LawDialog(QDialog):
    def __init__(self, law=Law(), unit="", duration=2.0, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Time law [{unit}]")
        self.resize(560, 480)
        self.law = law
        self.duration = duration
        self.unit = unit
        layout = QVBoxLayout(self)
        self.kind = QComboBox()
        for value, label in LAW_LABELS.items():
            self.kind.addItem(label, value)
        layout.addWidget(self.kind)
        self.hint = QLabel()
        self.hint.setWordWrap(True)
        layout.addWidget(self.hint)
        self.input = QPlainTextEdit()
        self.input.setMaximumHeight(150)
        layout.addWidget(self.input)
        self.error = QLabel()
        self.error.setWordWrap(True)
        layout.addWidget(self.error)
        self.plot = SeriesView()
        layout.addWidget(self.plot)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        preview = QPushButton("Update preview")
        preview.clicked.connect(self.preview)
        layout.addWidget(preview)
        layout.addWidget(buttons)
        self.kind.currentIndexChanged.connect(
            lambda _: self._type_changed(self.kind.currentData())
        )
        self.kind.setCurrentIndex(self.kind.findData(law.kind))
        self._type_changed(law.kind)
        self.input.setPlainText(
            "\n".join(
                numeric_text(law.values[i : i + 2])
                for i in range(0, len(law.values), 2)
            )
            if law.kind == "table"
            else numeric_text(law.values)
        )
        self.preview()

    def _type_changed(self, kind):
        self.hint.setText(
            {
                "constante": f"Value [{self.unit}]",
                "lineaire": f"Initial value [{self.unit}], slope [{self.unit}/s]",
                "table": f"One time [s], value [{self.unit}] pair per line. Linear interpolation; endpoint values are held outside the table.",
            }[kind]
        )
        self.input.setPlainText(
            {
                "constante": "0",
                "lineaire": "0, 1",
                "table": f"0, 0\n{self.duration:g}, 1",
            }[kind]
        )

    def read(self):
        text = self.input.toPlainText()
        if self.kind.currentData() == "table":
            values = tuple(
                v for row in text.splitlines() if row.strip() for v in numbers(row, 2)
            )
        else:
            values = numbers(text, 1 if self.kind.currentData() == "constante" else 2)
        return Law(self.kind.currentData(), values)

    def preview(self):
        try:
            law = self.read()
            times = np.linspace(0, self.duration, 201)
            values = np.array([law.value(t) for t in times])
            if not np.isfinite(values).all():
                raise ValueError(
                    "The law exceeds the finite-number range over this duration."
                )
            self.plot.set_series(times, values, "Prescribed value", self.unit)
            self.error.setText("")
            return law
        except ValueError as exc:
            self.error.setText(str(exc))
            return None

    def _accept(self):
        law = self.preview()
        if law is not None:
            self.law = law
            self.accept()


class JointDialog(QDialog):
    def __init__(self, project, selected=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Create joint")
        self.project = project
        self.joint = None
        layout = QVBoxLayout(self)
        form = QFormLayout()
        layout.addLayout(form)
        self.kind = QComboBox()
        for value, label in JOINT_LABELS.items():
            self.kind.addItem(label, value)
        form.addRow("Type", self.kind)
        self.a = QComboBox()
        self.b = QComboBox()
        for combo in (self.a, self.b):
            combo.addItem("Ground", None)
            for body in project.bodies:
                combo.addItem(body.name, body.id)
        self.b.setCurrentIndex(max(1, self.b.findData(selected)))
        form.addRow("Body A", self.a)
        form.addRow("Body B", self.b)
        self.point = QLineEdit("0, 0, 0")
        self.axis = QLineEdit("0, 0, 1")
        form.addRow("World point X,Y,Z [m]", self.point)
        form.addRow("World axis X,Y,Z", self.axis)
        note = QLabel(
            "Both local anchors are defined at the chosen point. Their placement remains explicit after any body modification."
        )
        note.setWordWrap(True)
        layout.addWidget(note)
        self.error = QLabel()
        self.error.setWordWrap(True)
        layout.addWidget(self.error)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _accept(self):
        try:
            self.joint = joint_at(
                self.project,
                self.kind.currentData(),
                self.a.currentData(),
                self.b.currentData(),
                numbers(self.point.text(), 3),
                numbers(self.axis.text(), 3),
            )
            self.accept()
        except (ValueError, StopIteration) as exc:
            self.error.setText(str(exc))
