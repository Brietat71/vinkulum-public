"""Bound Cocoa paint/render feedback without replacing the native VTK window.

Qt 6.10+ / VTK Cocoa can request a paint as a consequence of rendering itself:
https://discourse.vtk.org/t/16225
Only an explicit update, resize or show should request another paint render.
Camera/manipulator and simulation updates still render directly through VTK.
"""

import math
import sys

from PySide6.QtCore import QEvent, QPoint, Qt, QTimer, Signal
from PySide6.QtWidgets import QApplication, QMenu
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class RenderInteractor(QVTKRenderWindowInteractor):
    clicked = Signal(QPoint)
    pressed = Signal()

    def __init__(self, *args, **kwargs):
        self._coalesce_paints = sys.platform == "darwin"
        self._paint_requested = True
        self._painting = False
        super().__init__(*args, **kwargs)
        self._press_point = None
        self._pointer_dragged = False
        self._context_press = False
        self.context_handler = self._navigation_menu
        self._navigation_timer = QTimer(self)
        self._navigation_timer.setSingleShot(True)
        self._navigation_timer.setInterval(16)
        self._navigation_timer.timeout.connect(self._render_navigation)

    def mousePressEvent(self, event):
        self._press_point = event.position().toPoint()
        self._pointer_dragged = False
        self._context_press = event.button() == Qt.MouseButton.RightButton or (
            sys.platform == "darwin"
            and event.button() == Qt.MouseButton.LeftButton
            and bool(event.modifiers() & Qt.KeyboardModifier.MetaModifier)
        )  # Qt maps the physical macOS Control key to MetaModifier.
        self.setFocus(Qt.FocusReason.MouseFocusReason)
        if self._context_press:
            event.accept()
            return
        self.pressed.emit()
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if (
            self._press_point is not None
            and (event.position().toPoint() - self._press_point).manhattanLength()
            >= QApplication.startDragDistance()
        ):
            self._pointer_dragged = True
        if self._context_press:
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        point = event.position().toPoint()
        click = self._press_point is not None and not self._pointer_dragged
        self._press_point = None
        if self._context_press:
            self._context_press = False
            if click:
                self.context_handler(self.mapToGlobal(point))
            event.accept()
            return
        super().mouseReleaseEvent(event)
        if (
            click
            and event.button() == Qt.MouseButton.LeftButton
            and not event.modifiers()
        ):
            self.clicked.emit(point)

    def contextMenuEvent(self, event):
        # Pointer menus are dispatched after release on all platforms, including
        # QTest and macOS Control-click. Do not also open Qt's synthesized menu.
        if event.reason() == event.Reason.Keyboard:
            self.context_handler(self.mapToGlobal(self.rect().center()))
        event.accept()

    def _renderer(self):
        return self._RenderWindow.GetRenderers().GetFirstRenderer()

    def _navigation_menu(self, global_point):
        renderer = self._renderer()
        if renderer is None:
            return
        menu = QMenu(self)
        menu.addAction("Fit visible scene", lambda: self._fit_visible(renderer))
        menu.addAction("Zoom in", lambda: self._zoom(1.2))
        menu.addAction("Zoom out", lambda: self._zoom(1 / 1.2))
        menu.aboutToHide.connect(menu.deleteLater)
        menu.popup(global_point)

    def _fit_visible(self, renderer):
        renderer.ResetCamera()
        self._queue_navigation()

    def _queue_navigation(self):
        if not self._navigation_timer.isActive():
            self._navigation_timer.start()

    def _render_navigation(self):
        if self.__dict__.get("_finalized", False) or not self.isVisible():
            return
        renderer = self._renderer()
        if renderer is not None:
            renderer.ResetCameraClippingRange()
            self._RenderWindow.Render()

    def _zoom(self, factor):
        renderer = self._renderer()
        if renderer is None or not math.isfinite(factor) or factor <= 0:
            return
        camera = renderer.GetActiveCamera()
        if camera.GetParallelProjection():
            camera.SetParallelScale(camera.GetParallelScale() / factor)
        else:
            camera.Dolly(factor)
        self._queue_navigation()

    def _pan_pixels(self, delta):
        renderer = self._renderer()
        if renderer is None:
            return
        camera = renderer.GetActiveCamera()
        renderer.SetWorldPoint(*camera.GetFocalPoint(), 1.0)
        renderer.WorldToDisplay()
        x, y, z = renderer.GetDisplayPoint()
        ratio = self._getPixelRatio()
        renderer.SetDisplayPoint(x - delta.x() * ratio, y + delta.y() * ratio, z)
        renderer.DisplayToWorld()
        world = renderer.GetWorldPoint()
        if not world[3]:
            return
        shift = tuple(world[i] / world[3] - camera.GetFocalPoint()[i] for i in range(3))
        camera.SetPosition(*(a + b for a, b in zip(camera.GetPosition(), shift)))
        camera.SetFocalPoint(*(a + b for a, b in zip(camera.GetFocalPoint(), shift)))
        self._queue_navigation()

    def wheelEvent(self, event):
        # Cocoa trackpads provide continuous pixel motion. The upstream VTK
        # bridge discards it and waits for whole 120-unit wheel notches.
        if sys.platform == "darwin" and not event.pixelDelta().isNull():
            self._pan_pixels(event.pixelDelta())
        else:
            self._zoom(
                math.exp(max(-5.0, min(5.0, event.angleDelta().y() / 120 * 0.18)))
            )
        event.accept()

    def event(self, event):
        if event.type() == QEvent.Type.NativeGesture:
            kind = event.gestureType()
            if kind == Qt.NativeGestureType.ZoomNativeGesture:
                self._zoom(1.0 + event.value())
            elif kind not in (
                Qt.NativeGestureType.BeginNativeGesture,
                Qt.NativeGestureType.EndNativeGesture,
            ):
                return super().event(event)
            event.accept()
            return True
        return super().event(event)

    def __getattr__(self, name):
        # Qt can emit parent.destroyed after cyclic GC has cleared this Python
        # wrapper. The upstream forwarding implementation recursively accesses
        # _Iren in that state. Missing state must be an ordinary AttributeError.
        interactor = self.__dict__.get("_Iren")
        if interactor is None:
            raise AttributeError(name)
        if name == "__vtk__":
            return lambda: interactor
        return getattr(interactor, name)

    def Finalize(self):
        # Explicit shutdown and the parent's destroyed -> close connection can
        # both reach here. Release the GL context once, before wrapper teardown.
        window = self.__dict__.get("_RenderWindow")
        if window is not None and not self.__dict__.get("_finalized", False):
            self._finalized = True
            timer = self.__dict__.get("_Timer")
            if timer is not None:
                timer.stop()
            navigation = self.__dict__.get("_navigation_timer")
            if navigation is not None:
                navigation.stop()
            window.Finalize()

    def update(self, *args):
        self._paint_requested = True
        return super().update(*args)

    def showEvent(self, event):
        self._paint_requested = True
        super().showEvent(event)

    def paintEvent(self, event):
        if not self._coalesce_paints:
            return super().paintEvent(event)
        if self._painting or not self._paint_requested:
            return
        self._paint_requested = False
        self._painting = True
        try:
            super().paintEvent(event)
        finally:
            self._painting = False
