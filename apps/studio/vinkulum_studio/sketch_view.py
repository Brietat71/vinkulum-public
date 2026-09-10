"""Precise planar sketch interaction in logical screen pixels; geometry stays in mm."""

import math

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import QWidget

from .sketch import settled, solve


class SketchCanvas(QWidget):
    picked = Signal(object, bool)
    add_point = Signal(object, object)
    moved = Signal(object)
    diagnostic = Signal(str)
    tool_changed = Signal(str)
    dimension_activated = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAccessibleName("Planar sketch canvas, millimetres")
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMouseTracking(True)
        self.setMinimumSize(440, 350)
        self.sketch = self.preview = None
        self.selection = ()
        self.center, self.zoom = (40.0, 30.0), 6.0
        self.tool = "select"
        self.snap_mm = 5.0
        self.last_point = None
        self.cursor_world = (0, 0)
        self._pan = self._drag = None
        self._fitted = False
        self._dimension_hits = []
        self.conflicts = ()

    def set_sketch(self, sketch, selection=()):
        if sketch != self.sketch:
            self._dimension_hits = []
        self.sketch, self.selection = sketch, tuple(selection)
        if self.last_point not in {p.id for p in sketch.points}:
            self.last_point = None
        self.update()

    def set_tool(self, tool):
        self.tool = tool
        self.last_point = None
        self.preview = None
        self.setCursor(
            Qt.CursorShape.CrossCursor
            if tool == "polyline"
            else Qt.CursorShape.ArrowCursor
        )
        self.tool_changed.emit(tool)
        self.setFocus()
        self.update()

    def screen(self, point):
        return QPointF(
            self.width() / 2 + (point[0] - self.center[0]) * self.zoom,
            self.height() / 2 - (point[1] - self.center[1]) * self.zoom,
        )

    def world(self, point):
        return (
            self.center[0] + (point.x() - self.width() / 2) / self.zoom,
            self.center[1] - (point.y() - self.height() / 2) / self.zoom,
        )

    def fit(self):
        if self.sketch and self.sketch.points:
            lo = [min(p.xy_mm[i] for p in self.sketch.points) for i in (0, 1)]
            hi = [max(p.xy_mm[i] for p in self.sketch.points) for i in (0, 1)]
        else:
            lo, hi = (0, 0), (100, 70)
        self.center = tuple((a + b) / 2 for a, b in zip(lo, hi))
        self.zoom = min(
            (self.width() - 140) / max(20, hi[0] - lo[0]),
            (self.height() - 140) / max(20, hi[1] - lo[1]),
        )
        self.zoom = min(300, max(0.001, self.zoom))
        self._fitted = True
        self.update()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if not self._fitted:
            self.fit()

    def hit(self, position):
        if self.sketch is None:
            return None
        for rectangle, identifier in reversed(self._dimension_hits):
            if rectangle.contains(position):
                return "constraint", identifier
        points = {p.id: self.screen(p.xy_mm) for p in self.sketch.points}
        nearest = sorted(
            (
                (math.hypot(q.x() - position.x(), q.y() - position.y()), identifier)
                for identifier, q in points.items()
            )
        )
        if nearest and nearest[0][0] <= 9:
            return "point", nearest[0][1]
        candidates = []
        for edge in self.sketch.edges:
            a, b = points[edge.start], points[edge.end]
            dx, dy = b.x() - a.x(), b.y() - a.y()
            length = dx * dx + dy * dy
            t = (
                max(
                    0,
                    min(
                        1,
                        ((position.x() - a.x()) * dx + (position.y() - a.y()) * dy)
                        / length,
                    ),
                )
                if length
                else 0
            )
            distance = math.hypot(
                position.x() - a.x() - t * dx, position.y() - a.y() - t * dy
            )
            candidates.append((distance, edge.id))
        if candidates and min(candidates)[0] <= 6:
            return "edge", min(candidates)[1]
        return None

    def target(self, position):
        hit = self.hit(position)
        if hit and hit[0] == "point":
            return next(p.xy_mm for p in self.sketch.points if p.id == hit[1]), hit[1]
        point = self.world(position)
        if self.snap_mm:
            point = tuple(round(v / self.snap_mm) * self.snap_mm for v in point)
        return point, None

    def mousePressEvent(self, event):
        self.setFocus()
        if event.button() == Qt.MouseButton.MiddleButton:
            self._pan = (event.position(), self.center)
        elif event.button() == Qt.MouseButton.LeftButton:
            if self.tool == "polyline":
                point, existing = self.target(event.position())
                self.add_point.emit(point, existing)
                return
            key = self.hit(event.position())
            additive = bool(
                event.modifiers()
                & (
                    Qt.KeyboardModifier.ControlModifier
                    | Qt.KeyboardModifier.ShiftModifier
                )
            )
            self.picked.emit(key, additive)
            if key and key[0] == "point" and not additive:
                try:
                    result = solve(self.sketch)
                except ValueError as error:
                    self.diagnostic.emit(str(error))
                    return
                if result.conflicts:
                    self.diagnostic.emit(
                        "Resolve the highlighted constraints before dragging a point."
                    )
                elif any(result.free_axes[key[1]]):
                    self._drag = (key[1], event.position())
                else:
                    self.diagnostic.emit(
                        "This point is driven by dimensions. Select a dimension to edit its value."
                    )

    def mouseMoveEvent(self, event):
        self.cursor_world = self.world(event.position())
        if self._pan:
            start, center = self._pan
            delta = event.position() - start
            self.center = (
                center[0] - delta.x() / self.zoom,
                center[1] + delta.y() / self.zoom,
            )
        elif self._drag and (event.position() - self._drag[1]).manhattanLength() > 3:
            try:
                point, _ = self.target(event.position())
                self.preview = settled(self.sketch, drag=(self._drag[0], point))
            except ValueError as error:
                self.diagnostic.emit(str(error))
        self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.MiddleButton:
            self._pan = None
        elif event.button() == Qt.MouseButton.LeftButton:
            if self.preview is not None and self.preview != self.sketch:
                self.moved.emit(self.preview)
            self.preview = self._drag = None
        self.update()

    def wheelEvent(self, event):
        point = self.world(event.position())
        self.zoom = min(
            300, max(0.001, self.zoom * 1.2 ** (event.angleDelta().y() / 120))
        )
        offset = event.position() - QPointF(self.width() / 2, self.height() / 2)
        self.center = (
            point[0] - offset.x() / self.zoom,
            point[1] + offset.y() / self.zoom,
        )
        event.accept()
        self.update()

    def mouseDoubleClickEvent(self, event):
        key = self.hit(event.position())
        if key and key[0] == "constraint":
            self.dimension_activated.emit(key[1])

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_F:
            self.fit()
        elif event.key() == Qt.Key.Key_L:
            self.set_tool("polyline")
        elif event.key() in (Qt.Key.Key_S, Qt.Key.Key_Escape) and self.tool != "select":
            self.set_tool("select")
        else:
            super().keyPressEvent(event)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#181d24"))
        span = 60 / self.zoom
        power = 10 ** math.floor(math.log10(span))
        step = next(v * power for v in (1, 2, 5, 10) if v * power >= span)
        lo, hi = (
            self.world(QPointF(0, self.height())),
            self.world(QPointF(self.width(), 0)),
        )
        painter.setPen(QPen(QColor("#28303b"), 1))
        for axis in (0, 1):
            for index in range(
                math.ceil(lo[axis] / step), math.floor(hi[axis] / step) + 1
            ):
                value = index * step
                a, b = (
                    ((value, lo[1]), (value, hi[1]))
                    if axis == 0
                    else ((lo[0], value), (hi[0], value))
                )
                painter.drawLine(self.screen(a), self.screen(b))
        painter.setPen(QPen(QColor("#643a43"), 1))
        painter.drawLine(self.screen((lo[0], 0)), self.screen((hi[0], 0)))
        painter.setPen(QPen(QColor("#335f51"), 1))
        painter.drawLine(self.screen((0, lo[1])), self.screen((0, hi[1])))
        sketch = self.preview or self.sketch
        if sketch is None:
            return
        points = {p.id: self.screen(p.xy_mm) for p in sketch.points}
        edges = {e.id: e for e in sketch.edges}
        selected = set(self.selection)
        for edge in sketch.edges:
            color = "#bac6ff" if ("edge", edge.id) in selected else "#91b9e8"
            painter.setPen(
                QPen(QColor(color), 2.8 if ("edge", edge.id) in selected else 1.6)
            )
            painter.drawLine(points[edge.start], points[edge.end])
        dim_indices = [0, 0]
        self._dimension_hits = []
        for c in sketch.constraints:
            painter.setPen(
                QPen(
                    QColor(
                        "#f39898"
                        if c.id in self.conflicts
                        else "#c4b9ef"
                        if ("constraint", c.id) in selected
                        else "#8293ac"
                    ),
                    1,
                )
            )
            if c.kind in ("horizontal", "vertical"):
                e = edges[c.entities[0]]
                middle = (points[e.start] + points[e.end]) / 2
                painter.drawText(
                    middle + QPointF(7, -8), "H" if c.kind == "horizontal" else "V"
                )
            elif c.kind == "fixed":
                q = points[c.entities[0]]
                painter.drawRect(QRectF(q.x() - 7, q.y() - 7, 14, 14))
            else:
                a, b = (points[key] for key in c.entities)
                axis = int(c.kind == "distance_y")
                offset = 28 + dim_indices[axis] * 22
                dim_indices[axis] += 1
                if axis == 0:
                    level = min(a.y(), b.y()) - offset
                    ends = QPointF(a.x(), level), QPointF(b.x(), level)
                    tick = QPointF(0, 4)
                else:
                    level = min(a.x(), b.x()) - offset
                    ends = QPointF(level, a.y()), QPointF(level, b.y())
                    tick = QPointF(4, 0)
                for point, end in zip((a, b), ends):
                    painter.drawLine(point, end)
                    painter.drawLine(end - tick, end + tick)
                painter.drawLine(*ends)
                middle = (ends[0] + ends[1]) / 2
                text = f"Δ{'XY'[axis]} {c.values_mm[0]} mm"
                width = painter.fontMetrics().horizontalAdvance(text)
                label = middle + QPointF(
                    -width / 2 if axis == 0 else 6, -5 if axis == 0 else 0
                )
                rectangle = QRectF(
                    label.x() - 3,
                    label.y() - painter.fontMetrics().ascent(),
                    width + 6,
                    painter.fontMetrics().height(),
                )
                painter.fillRect(rectangle, QColor("#181d24"))
                painter.drawText(label, text)
                self._dimension_hits.append((rectangle, c.id))
        for point in sketch.points:
            q = points[point.id]
            active = ("point", point.id) in selected
            painter.setPen(QPen(QColor("#d6deff" if active else "#90abca"), 1.5))
            painter.setBrush(QColor("#bac6ff" if active else "#181d24"))
            painter.drawEllipse(q, 5 if active else 3.5, 5 if active else 3.5)
        if self.tool == "polyline" and self.last_point in points:
            painter.setPen(QPen(QColor("#d3c593"), 1, Qt.PenStyle.DashLine))
            painter.drawLine(points[self.last_point], self.screen(self.cursor_world))
        painter.setPen(QColor("#a1adbf"))
        painter.drawText(
            12,
            self.height() - 12,
            f"XY · mm    X {self.cursor_world[0]:.6g}    Y {self.cursor_world[1]:.6g}    Grid {step:g} mm",
        )
