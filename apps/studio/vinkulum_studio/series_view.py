"""Small scientific curve drawn by Qt; no plot engine owns the samples."""

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QWidget


class SeriesView(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(170)
        self.times = ()
        self.values = ()
        self.label = "Aucun résultat"
        self.unit = ""
        self.index = 0

    def set_series(self, times, values, label, unit):
        self.times, self.values, self.label, self.unit = times, values, label, unit
        self.index = min(self.index, max(0, len(times) - 1))
        self.update()

    def paintEvent(self, event):
        with QPainter(self) as p:
            p.fillRect(self.rect(), QColor("#f8fafc"))
            p.setRenderHint(QPainter.RenderHint.Antialiasing)
            p.setPen(QColor("#24374a"))
            p.drawText(
                12, 22, f"{self.label} [{self.unit}]" if self.unit else self.label
            )
            if len(self.times) < 2:
                return
            rect = QRectF(66, 36, self.width() - 90, self.height() - 72)
            low, high = float(min(self.values)), float(max(self.values))
            margin = max((high - low) * 0.1, 1e-6)
            low -= margin
            high += margin
            duration = max(float(self.times[-1] - self.times[0]), 1e-12)

            def point(t, y):
                return QPointF(
                    rect.left() + (t - self.times[0]) / duration * rect.width(),
                    rect.bottom() - (y - low) / (high - low) * rect.height(),
                )

            for value in (low, (low + high) / 2, high):
                y = point(self.times[0], value).y()
                p.setPen(QColor("#d4dde5"))
                p.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
                p.setPen(QColor("#42556a"))
                p.drawText(
                    QRectF(0, y - 8, 60, 20),
                    Qt.AlignmentFlag.AlignRight,
                    f"{value:.3g}",
                )
            path = QPainterPath(point(self.times[0], self.values[0]))
            for t, y in zip(self.times[1:], self.values[1:]):
                path.lineTo(point(t, y))
            p.setPen(QPen(QColor("#187b8d"), 2))
            p.drawPath(path)
            cursor = point(self.times[self.index], self.values[self.index])
            p.setPen(QPen(QColor("#bd6425"), 1, Qt.PenStyle.DashLine))
            p.drawLine(
                QPointF(cursor.x(), rect.top()), QPointF(cursor.x(), rect.bottom())
            )
            p.setPen(QColor("#24374a"))
            p.drawText(QPointF(rect.left(), rect.bottom() + 22), f"{self.times[0]:g}")
            p.drawText(
                QPointF(rect.right() - 100, rect.bottom() + 22), f"{self.times[-1]:g} s"
            )
