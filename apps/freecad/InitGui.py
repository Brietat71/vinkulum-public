"""FreeCAD workbench entry point. SPDX-License-Identifier: Apache-2.0"""

import FreeCADGui as Gui


class VinkulumWorkbench(Workbench):
    MenuText = "Vinkulum"
    ToolTip = "Mécanique et résultats Vinkulum dans FreeCAD"

    def __init__(self):
        import vinkulum_freecad

        self.Icon = vinkulum_freecad.ICON

    def Initialize(self):
        import vinkulum_freecad

        commands = vinkulum_freecad.register_commands()
        self.appendToolbar("Vinkulum", commands)
        self.appendMenu("Vinkulum", commands)

    def GetClassName(self):
        return "Gui::PythonWorkbench"


Gui.addWorkbench(VinkulumWorkbench())
