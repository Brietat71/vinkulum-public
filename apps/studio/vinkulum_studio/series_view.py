"""Interactive native-sample plot and lazy table; neither alters scientific data."""

from bisect import bisect_left, bisect_right

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QWidget

from .theme import THEMES


class SampleTableModel(QAbstractTableModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.times, self.values = (), ()
        self.label, self.unit = "Value", ""

    def set_series(self, times, values, label, unit):
        self.beginResetModel()
        self.times, self.values, self.label, self.unit = times, values, label, unit
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.times)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else 3

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or not 0 <= index.row() < len(self.times):
            return None
        if role == Qt.ItemDataRole.DisplayRole:
            if index.column() == 0:
                return str(index.row())
            return format(
                float(
                    (self.times if index.column() == 1 else self.values)[index.row()]
                ),
                ".12g",
            )
        if role == Qt.ItemDataRole.TextAlignmentRole:
            return Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        return None

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if (
            role == Qt.ItemDataRole.DisplayRole
            and orientation == Qt.Orientation.Horizontal
        ):
            return ("Sample", "Time [s]", f"{self.label} [{self.unit}]")[section]
        return super().headerData(section, orientation, role)


class SeriesView(QWidget):
    sample_selected = Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(145)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setAccessibleName("Native-sample curve")
        self.setAccessibleDescription(
            "Click to select a time. Scroll to zoom, Shift-drag to pan, double-click to fit. Arrow keys move between samples."
        )
        self.setToolTip(self.accessibleDescription())
        self.times, self.values = (), ()
        self.label, self.unit = "No result", ""
        self.index = 0
        self.domain = (0.0, 1.0)
        self._cache = None
        self.reference = None
        self._pan = None
        self.colors = THEMES["dark"]

    def set_theme(self, name):
        self.colors = THEMES[name]
        self.update()

    def set_series(self, times, values, label, unit):
        self.times, self.values, self.label, self.unit = times, values, label, unit
        self.index = min(self.index, max(0, len(times) - 1))
        self.reset_view()

    def set_reference(self, reference):
        self.reference = reference
        self._cache = None
        self.update()

    def reset_view(self):
        self.domain = (
            (float(self.times[0]), float(self.times[-1]))
            if len(self.times) > 1
            else (0.0, 1.0)
        )
        self._cache = None
        self.update()

    def plot_rect(self):
        return QRectF(78, 40, max(1, self.width() - 100), max(1, self.height() - 80))

    def _geometry(self):
        key = (self.width(), self.height(), self.domain)
        if self._cache is not None and self._cache[0] == key:
            return self._cache[1:]
        rect = self.plot_rect()
        start = max(0, bisect_left(self.times, self.domain[0]) - 1)
        stop = min(len(self.times), bisect_right(self.times, self.domain[1]) + 1)
        values = self.values[start:stop]
        low, high = float(min(values)), float(max(values))
        reference_indices = (0, 0)
        if self.reference is not None:
            rt, rv, _ = self.reference
            ra = max(0, bisect_left(rt, self.domain[0]) - 1)
            rb = min(len(rt), bisect_right(rt, self.domain[1]) + 1)
            reference_indices = (ra, rb)
            if rb > ra:
                low, high = (
                    min(low, float(min(rv[ra:rb]))),
                    max(high, float(max(rv[ra:rb]))),
                )
        margin = max((high - low) * 0.1, abs(high) * 1e-9, 1e-9)
        low, high = low - margin, high + margin
        path = QPainterPath()
        for i in range(start, stop):
            point = self._point(self.times[i], self.values[i], low, high)
            if i == start:
                path.moveTo(point)
            else:
                path.lineTo(point)
        reference_path = QPainterPath()
        if self.reference is not None:
            rt, rv, _ = self.reference
            ra, rb = reference_indices
            for i in range(ra, rb):
                point = self._point(rt[i], rv[i], low, high)
                if i == ra:
                    reference_path.moveTo(point)
                else:
                    reference_path.lineTo(point)
        self._cache = (key, rect, low, high, path, reference_path)
        return self._cache[1:]

    def _point(self, t, y, low, high):
        rect = self.plot_rect()
        return QPointF(
            rect.left()
            + (float(t) - self.domain[0])
            / max(self.domain[1] - self.domain[0], 1e-15)
            * rect.width(),
            rect.bottom() - (float(y) - low) / (high - low) * rect.height(),
        )

    def _time_at(self, x):
        rect = self.plot_rect()
        return self.domain[0] + (x - rect.left()) / rect.width() * (
            self.domain[1] - self.domain[0]
        )

    def _select_at(self, x):
        if not len(self.times):
            return
        t = self._time_at(x)
        right = min(bisect_left(self.times, t), len(self.times) - 1)
        left = max(0, right - 1)
        self.sample_selected.emit(
            left if abs(self.times[left] - t) <= abs(self.times[right] - t) else right
        )

    def mousePressEvent(self, event):
        self.setFocus()
        if event.button() == Qt.MouseButton.LeftButton and self.plot_rect().contains(
            event.position()
        ):
            if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
                self._pan = (event.position().x(), self.domain)
            else:
                self._select_at(event.position().x())
            event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._pan is not None:
            x, domain = self._pan
            delta = (
                (x - event.position().x())
                / self.plot_rect().width()
                * (domain[1] - domain[0])
            )
            self._set_domain(domain[0] + delta, domain[1] + delta)
        elif event.buttons() & Qt.MouseButton.LeftButton:
            self._select_at(event.position().x())
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        self._pan = None
        super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event):
        self.reset_view()
        event.accept()

    def _set_domain(self, low, high):
        if len(self.times) < 2:
            return
        full_low, full_high = float(self.times[0]), float(self.times[-1])
        width = min(high - low, full_high - full_low)
        low = max(full_low, min(low, full_high - width))
        self.domain = (low, low + width)
        self._cache = None
        self.update()

    def wheelEvent(self, event):
        if len(self.times) < 2 or not self.plot_rect().contains(event.position()):
            event.ignore()
            return
        factor = 0.8 if event.angleDelta().y() > 0 else 1.25
        anchor = self._time_at(event.position().x())
        low, high = self.domain
        if (high - low) * factor >= (self.times[-1] - self.times[0]) * 1e-6:
            self._set_domain(
                anchor + (low - anchor) * factor, anchor + (high - anchor) * factor
            )
        event.accept()

    def keyPressEvent(self, event):
        keys = {
            Qt.Key.Key_Left: self.index - 1,
            Qt.Key.Key_Right: self.index + 1,
            Qt.Key.Key_Home: 0,
            Qt.Key.Key_End: len(self.times) - 1,
        }
        if event.key() in keys and len(self.times):
            self.sample_selected.emit(
                max(0, min(keys[event.key()], len(self.times) - 1))
            )
            event.accept()
        else:
            super().keyPressEvent(event)

    def paintEvent(self, event):
        c = self.colors
        with QPainter(self) as p:
            p.fillRect(self.rect(), QColor(c["field"]))
            p.setRenderHint(QPainter.RenderHint.Antialiasing)
            p.setPen(QColor(c["text"]))
            if len(self.times) < 2:
                p.drawText(
                    self.rect(),
                    Qt.AlignmentFlag.AlignCenter,
                    "Curves will appear after the first calculation.",
                )
                return
            rect, low, high, path, reference_path = self._geometry()
            value = float(self.values[self.index])
            title = f"{self.label} [{self.unit}] · {value:.6g}"
            p.drawText(
                12,
                23,
                p.fontMetrics().elidedText(
                    title, Qt.TextElideMode.ElideMiddle, self.width() - 24
                ),
            )
            for value in (low, (low + high) / 2, high):
                y = self._point(self.domain[0], value, low, high).y()
                p.setPen(QColor(c["grid"]))
                p.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
                p.setPen(QColor(c["muted"]))
                p.drawText(
                    QRectF(0, y - 9, 70, 22),
                    Qt.AlignmentFlag.AlignRight,
                    f"{value:.4g}",
                )
            p.save()
            p.setClipRect(rect)
            p.setPen(QPen(QColor(c["curve"]), 2))
            p.drawPath(path)
            p.setPen(QPen(QColor(c["cursor"]), 1.5, Qt.PenStyle.DashLine))
            p.drawPath(reference_path)
            cursor = self._point(
                self.times[self.index], self.values[self.index], low, high
            )
            p.setPen(QPen(QColor(c["cursor"]), 1, Qt.PenStyle.DashLine))
            p.drawLine(
                QPointF(cursor.x(), rect.top()), QPointF(cursor.x(), rect.bottom())
            )
            p.setBrush(QColor(c["cursor"]))
            p.drawEllipse(cursor, 3, 3)
            p.restore()
            p.setPen(QColor(c["muted"]))
            for fraction in (0.0, 0.5, 1.0):
                t = self.domain[0] + fraction * (self.domain[1] - self.domain[0])
                p.drawText(
                    QRectF(
                        rect.left() + fraction * rect.width() - 35,
                        rect.bottom() + 7,
                        70,
                        20,
                    ),
                    Qt.AlignmentFlag.AlignCenter,
                    f"{t:.4g} s",
                )
            if self.hasFocus():
                p.setPen(QPen(QColor(c["accent"]), 1))
                p.drawRect(self.rect().adjusted(0, 0, -1, -1))
