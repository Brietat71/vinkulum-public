"""Typed authoring dialogs. Inputs are numbers and tables, never executable code."""

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
import numpy as np

from .document import Law, joint_at
from .model import finite_number
from .series_view import SeriesView


def numbers(text, count=None):
    try:
        values = tuple(float(v.strip()) for v in text.replace(";", ",").split(","))
    except ValueError as exc:
        raise ValueError("Saisir des nombres séparés par des virgules.") from exc
    if (
        not all(finite_number(v) for v in values)
        or count is not None
        and len(values) != count
    ):
        raise ValueError(f"{count or 'Des'} valeurs finies requises.")
    return values


def numeric_text(values):
    return ", ".join(repr(float(v)) for v in values)


class LawDialog(QDialog):
    def __init__(self, law=Law(), unit="", duration=2.0, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Loi temporelle [{unit}]")
        self.resize(560, 480)
        self.law = law
        self.duration = duration
        self.unit = unit
        layout = QVBoxLayout(self)
        self.kind = QComboBox()
        self.kind.addItems(["constante", "lineaire", "table"])
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
        preview = QPushButton("Actualiser l'aperçu")
        preview.clicked.connect(self.preview)
        layout.addWidget(preview)
        layout.addWidget(buttons)
        self.kind.currentTextChanged.connect(self._type_changed)
        self.kind.setCurrentText(law.kind)
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
                "constante": f"Valeur [{self.unit}]",
                "lineaire": f"Valeur initiale [{self.unit}], pente [{self.unit}/s]",
                "table": f"Une paire temps [s], valeur [{self.unit}] par ligne. Interpolation linéaire, valeurs maintenues hors table.",
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
        if self.kind.currentText() == "table":
            values = tuple(
                v for row in text.splitlines() if row.strip() for v in numbers(row, 2)
            )
        else:
            values = numbers(text, 1 if self.kind.currentText() == "constante" else 2)
        return Law(self.kind.currentText(), values)

    def preview(self):
        try:
            law = self.read()
            times = np.linspace(0, self.duration, 201)
            values = np.array([law.value(t) for t in times])
            if not np.isfinite(values).all():
                raise ValueError(
                    "La loi dépasse le domaine des nombres finis sur cette durée."
                )
            self.plot.set_series(times, values, "Consigne", self.unit)
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
        self.setWindowTitle("Créer une liaison")
        self.project = project
        self.joint = None
        layout = QVBoxLayout(self)
        form = QFormLayout()
        layout.addLayout(form)
        self.kind = QComboBox()
        self.kind.addItems(["pivot", "rotule", "glissiere", "encastrement"])
        form.addRow("Type", self.kind)
        self.a = QComboBox()
        self.b = QComboBox()
        for combo in (self.a, self.b):
            combo.addItem("Bâti", None)
            for body in project.bodies:
                combo.addItem(body.name, body.id)
        self.b.setCurrentIndex(max(1, self.b.findData(selected)))
        form.addRow("Corps A", self.a)
        form.addRow("Corps B", self.b)
        self.point = QLineEdit("0, 0, 0")
        self.axis = QLineEdit("0, 0, 1")
        form.addRow("Point monde X,Y,Z [m]", self.point)
        form.addRow("Axe monde X,Y,Z", self.axis)
        note = QLabel(
            "Les deux ancrages locaux sont définis au point choisi. Leur placement reste explicite après toute modification des corps."
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
                self.kind.currentText(),
                self.a.currentData(),
                self.b.currentData(),
                numbers(self.point.text(), 3),
                numbers(self.axis.text(), 3),
            )
            self.accept()
        except (ValueError, StopIteration) as exc:
            self.error.setText(str(exc))
