#!/usr/bin/env python3
"""Time real workspace commands with explicit boundaries and installed-source hashes.

Run in the installed Studio environment, under Xvfb on Linux or Cocoa on macOS:
    python ci/studio_interaction_recipe.py /tmp/interaction.json
No timing threshold is used as a CI gate.
"""

import cProfile
import hashlib
import json
import platform
import pstats
import sys
import time
from importlib.metadata import version
from pathlib import Path

from PySide6.QtCore import QLocale, QSettings, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio.editor import EditorWindow

output = Path(sys.argv[1]).resolve()
output.parent.mkdir(parents=True, exist_ok=True)
app = QApplication.instance() or QApplication([])
QLocale.setDefault(QLocale("en_GB"))
app.setApplicationName("Vinkulum Studio")
app.setStyle("Fusion")
app.setQuitOnLastWindowClosed(False)
root = Path(__file__).resolve().parents[1]
window = EditorWindow(
    QSettings(str(output.with_suffix(".ini")), QSettings.Format.IniFormat)
)
window._discard_allowed = lambda: True
window.resize(1280, 800)
window.move(0, 0)
window.show()
QTest.qWait(100)
rows = []
profiler = cProfile.Profile()


def measure(name, action):
    entered = time.perf_counter()
    responsive = []
    QTimer.singleShot(0, lambda: responsive.append(time.perf_counter()))
    profiler.enable()
    action()
    returned = time.perf_counter()
    QTest.qWait(100)
    profiler.disable()
    rows.append(
        {
            "action": name,
            "synchronous_ms": 1000 * (returned - entered),
            "next_event_ms": 1000 * (responsive[0] - entered),
        }
    )


try:
    measure("open-articulated-first", window.commands["articulated"].trigger)
    for i in range(4):
        measure("model-to-simulate", lambda: window.workspace_tabs.setCurrentIndex(1))
        measure("simulate-to-model", lambda: window.workspace_tabs.setCurrentIndex(0))
    measure("open-articulated-again", window.commands["articulated"].trigger)
    measure("open-statics-first", window.commands["static_study"].trigger)
    measure("open-statics-again", window.commands["static_study"].trigger)
    measure(
        "load-parametric-plate",
        lambda: window.load(
            root / "examples/studio/platine-parametrique.vinkulum.json"
        ),
    )
    window.select_object(window.project.bodies[0].id)
    measure("open-cad-study-first", window.commands["cad_study"].trigger)
    measure("open-cad-study-again", window.commands["cad_study"].trigger)
    import vinkulum_studio

    module_root = Path(vinkulum_studio.__file__).parent
    report = {
        "boundary": "QAction/tab action to synchronous return and next event-loop opportunity; excludes time until the final frame is ready",
        "platform": platform.platform(),
        "machine": platform.machine(),
        "window_size": [window.width(), window.height()],
        "viewport_size": [window.viewport.width(), window.viewport.height()],
        "frozen": bool(getattr(sys, "frozen", False)),
        "packages": {
            name: version(name)
            for name in ("vinkulum-studio", "vinkulum", "PySide6", "vtk")
        },
        "application_modules_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(module_root.glob("*.py"))
        },
        "opengl": window.viewport.view.GetRenderWindow().ReportCapabilities(),
        "rows": rows,
    }
    output.write_text(json.dumps(report, indent=2))
    for study in (
        window._articulated_window,
        window._static_window,
        window._mesh_window,
    ):
        if study is not None:
            study.hide()
    window.show()
    window.raise_()
    QTest.qWait(100)
    # QWidget.grab() cannot capture the native OpenGL child. Capture the actual
    # screen window, as the frozen-bundle workspace recipe does.
    if (
        not window.screen()
        .grabWindow(window.winId())
        .save(str(output.with_suffix(".png")))
    ):
        raise RuntimeError("Could not capture the visible workspace.")
    with output.with_suffix(".profile.txt").open("w") as stream:
        pstats.Stats(profiler, stream=stream).sort_stats("cumulative").print_stats(45)
    print(json.dumps(rows), flush=True)
finally:
    for widget in (
        window._articulated_window,
        window._static_window,
        window._mesh_window,
    ):
        if widget is not None:
            widget._discard_allowed = lambda: True
    window.close()
    QTest.qWait(50)
