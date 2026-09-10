"""Vinkulum FreeCAD extension, loaded without changing the active workbench."""

from PySide6.QtCore import QTimer
from vinkulum_freecad.host import install

QTimer.singleShot(0, install)
