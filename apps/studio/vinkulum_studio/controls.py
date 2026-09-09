"""Compact vector inputs and original scalable line icons for the desktop UI."""

from functools import lru_cache

from PySide6.QtCore import QByteArray, Qt, Signal
from PySide6.QtGui import QIcon, QPainter, QPixmap
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import QGridLayout, QHBoxLayout, QLabel, QLineEdit, QWidget


class VectorField(QWidget):
    textEdited = Signal(str)

    def __init__(self, values, labels, parent=None):
        super().__init__(parent)
        layout = QGridLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setHorizontalSpacing(6)
        layout.setVerticalSpacing(6)
        self.components = []
        self.original_values = tuple(values)
        self.original_texts = []
        for i, (value, label) in enumerate(zip(values, labels)):
            cell = QWidget()
            cell.setObjectName("vector_cell")
            row = QHBoxLayout(cell)
            row.setContentsMargins(7, 0, 0, 0)
            row.setSpacing(3)
            axis = QLabel(label)
            axis.setObjectName("axis_label")
            row.addWidget(axis)
            field = QLineEdit(format(value, ".6g"))
            self.original_texts.append(field.text())
            field.setCursorPosition(0)
            field.setObjectName("component")
            field.setMinimumWidth(20)
            field.setAccessibleName(label)
            field.setToolTip(f"{label} · {value:.17g}")
            field.textEdited.connect(lambda _: self.textEdited.emit(self.text()))
            row.addWidget(field, 1)
            self.components.append(field)
            layout.addWidget(cell, i // 3, i % 3)
        self.setFocusProxy(self.components[0])

    def text(self):
        return ", ".join(field.text() for field in self.components)

    def values(self):
        return tuple(
            original if field.text() == text else float(field.text())
            for field, text, original in zip(
                self.components, self.original_texts, self.original_values
            )
        )

    def setText(self, text):
        values = text.split(",")
        if len(values) != len(self.components):
            raise ValueError("Nombre de composantes incorrect.")
        for field, value in zip(self.components, values):
            field.setText(value.strip())

    def setAccessibleName(self, name):
        super().setAccessibleName(name)
        for field in self.components:
            field.setAccessibleName(name + " · " + field.accessibleName())


PATHS = {
    "open": '<path d="M3 7h6l2 3h10l-3 10H3z M3 7V4h7l2 3h7v3"/>',
    "save": '<path d="M4 3h14l3 3v15H3V3z M7 3v6h10V3 M7 21v-8h10v8"/>',
    "undo": '<path d="M9 5L3 11l6 6 M3 11h11a6 6 0 0 1 6 6"/>',
    "redo": '<path d="M15 5l6 6-6 6 M21 11H10a6 6 0 0 0-6 6"/>',
    "box": '<path d="M12 2l9 5v10l-9 5-9-5V7z M3 7l9 5 9-5 M12 12v10"/>',
    "cylinder": '<ellipse cx="12" cy="5" rx="8" ry="3"/><path d="M4 5v14c0 4 16 4 16 0V5"/>',
    "sphere": '<circle cx="12" cy="12" r="9"/><ellipse cx="12" cy="12" rx="4" ry="9"/><path d="M3 12h18"/>',
    "joint": '<circle cx="12" cy="12" r="5"/><path d="M12 2v5 M12 17v5 M2 12h5 M17 12h5"/>',
    "load": '<path d="M3 20L20 3 M10 3h10v10 M3 15v6h6"/>',
    "select": '<path d="M5 2l14 12-7 1-3 7z"/>',
    "move": '<path d="M12 2v20 M2 12h20 M8 6l4-4 4 4 M8 18l4 4 4-4 M6 8l-4 4 4 4 M18 8l4 4-4 4"/>',
    "rotate": '<path d="M20 8a9 9 0 1 0 0 8 M20 2v6h-6"/>',
    "fit": '<path d="M8 3H3v5 M16 3h5v5 M3 16v5h5 M21 16v5h-5"/><rect x="8" y="8" width="8" height="8" rx="1"/>',
    "eye": '<path d="M2 12s4-7 10-7 10 7 10 7-4 7-10 7-10-7-10-7z"/><circle cx="12" cy="12" r="3"/>',
    "search": '<circle cx="10" cy="10" r="6"/><path d="M15 15l6 6"/>',
    "play": '<path d="M7 3l14 9-14 9z"/>',
    "grid": '<path d="M3 3h18v18H3z M3 9h18 M3 15h18 M9 3v18 M15 3v18"/>',
    "settings": '<path d="M3 6h18 M3 12h18 M3 18h18"/><circle cx="8" cy="6" r="2"/><circle cx="16" cy="12" r="2"/><circle cx="10" cy="18" r="2"/>',
}


@lru_cache(maxsize=128)
def line_icon(name, color="#bbc1ce"):
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" '
        f'fill="none" stroke="{color}" stroke-width="1.55" stroke-linecap="round" stroke-linejoin="round">'
        + PATHS[name]
        + "</svg>"
    )
    renderer = QSvgRenderer(QByteArray(svg.encode()))
    icon = QIcon()
    for size in (16, 20, 24, 32, 48, 64, 96):
        pixmap = QPixmap(size, size)
        pixmap.fill(Qt.GlobalColor.transparent)
        with QPainter(pixmap) as painter:
            renderer.render(painter)
        icon.addPixmap(pixmap)
    return icon
