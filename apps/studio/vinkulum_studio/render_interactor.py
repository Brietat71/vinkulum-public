"""Bound Cocoa paint/render feedback without replacing the native VTK window.

Qt 6.10+ / VTK Cocoa can request a paint as a consequence of rendering itself:
https://discourse.vtk.org/t/16225
Only an explicit update, resize or show should request another paint render.
Camera/manipulator and simulation updates still render directly through VTK.
"""

import sys
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class RenderInteractor(QVTKRenderWindowInteractor):
    def __init__(self, *args, **kwargs):
        self._coalesce_paints = sys.platform == "darwin"
        self._paint_requested = True
        self._painting = False
        super().__init__(*args, **kwargs)

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
