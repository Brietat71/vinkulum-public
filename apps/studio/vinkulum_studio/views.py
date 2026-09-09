"""Native Qt drawings; frames are real samples, never synthetic solver data."""

import math

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import QWidget


class PendulumView(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(320, 230)
        self.length = 1.0
        self.sample = (0.0, 0.5, -math.sqrt(0.75), math.pi / 6)
        self.setAccessibleName("Pendulum view in the x-z plane")

    def paintEvent(self, event):
        with QPainter(self) as painter:
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.fillRect(self.rect(), QColor("#f4f7fb"))
            origin = QPointF(self.width() / 2, self.height() / 2)
            scale = min(self.width() * 0.39, (self.height() - 80) / 2) / self.length
            t, x, z, angle = self.sample
            tip = origin + QPointF(x * scale, -z * scale)
            painter.setPen(QPen(QColor("#cad4e1"), 1, Qt.PenStyle.DashLine))
            painter.drawLine(origin, origin + QPointF(0, self.length * scale))
            painter.setPen(
                QPen(
                    QColor("#465569"), 6, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap
                )
            )
            painter.drawLine(origin + QPointF(-25, 0), origin + QPointF(25, 0))
            painter.setPen(QPen(QColor("#187b8d"), 3))
            painter.drawLine(origin, tip)
            painter.setBrush(QColor("#187b8d"))
            painter.drawEllipse(tip, 10, 10)
            painter.setPen(QColor("#25354b"))
            painter.drawText(18, 25, "PENDULUM  /  X–Z PLANE")
            painter.drawText(
                18,
                self.height() - 16,
                f"L = {self.length:g} m     θ = {math.degrees(angle):.2f}°     t = {t:.4f} s",
            )


class AngleView(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(320, 200)
        self.samples = ()
        self.index = 0
        self.setAccessibleName("Angle versus time curve")

    def paintEvent(self, event):
        with QPainter(self) as painter:
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.fillRect(self.rect(), QColor("#ffffff"))
            painter.setPen(QColor("#25354b"))
            painter.drawText(18, 24, "ANGLE θ [°]")
            rect = QRectF(55, 42, self.width() - 80, self.height() - 83)
            if not self.samples:
                painter.drawText(
                    rect,
                    Qt.AlignmentFlag.AlignCenter,
                    "Run a calculation to display the trajectory.",
                )
                return
            duration = self.samples[-1][0]
            limit = max(
                1.0, max(abs(math.degrees(row[3])) for row in self.samples) * 1.12
            )

            def point(row):
                return QPointF(
                    rect.left() + row[0] / duration * rect.width(),
                    rect.center().y()
                    - math.degrees(row[3]) / limit * rect.height() / 2,
                )

            painter.setPen(QPen(QColor("#dce3ed"), 1))
            for value in (-limit, 0.0, limit):
                y = rect.center().y() - value / limit * rect.height() / 2
                painter.setPen(QPen(QColor("#dce3ed"), 1))
                painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
                painter.setPen(QColor("#526174"))
                painter.drawText(
                    QRectF(0, y - 9, 47, 20),
                    Qt.AlignmentFlag.AlignRight,
                    f"{value:.1f}",
                )
            path = QPainterPath(point(self.samples[0]))
            for row in self.samples[1:]:
                path.lineTo(point(row))
            painter.setPen(QPen(QColor("#187b8d"), 2))
            painter.drawPath(path)
            selected = point(self.samples[self.index])
            painter.setPen(QPen(QColor("#bd6425"), 1, Qt.PenStyle.DashLine))
            painter.drawLine(
                QPointF(selected.x(), rect.top()), QPointF(selected.x(), rect.bottom())
            )
            painter.setBrush(QColor("#bd6425"))
            painter.drawEllipse(selected, 4, 4)
            painter.setPen(QColor("#25354b"))
            painter.drawText(QPointF(rect.left(), rect.bottom() + 23), "0")
            painter.drawText(
                QPointF(rect.right() - 105, rect.bottom() + 23),
                f"{duration:g} s  /  TIME",
            )
