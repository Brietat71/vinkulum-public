"""Exercise the installed application and its own worker, including native GL."""

import json
from dataclasses import replace
from pathlib import Path
import platform
import sys

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QApplication

from .document import History
from .editor import EditorWindow
from .examples3d import double_pendulum


def main(directory):
    output = Path(directory)
    output.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])
    window = EditorWindow()
    window._discard_allowed = lambda: True
    project = replace(double_pendulum(), duration=0.1)
    window.history = History(project)
    window._refresh(True)
    report = {
        "status": "failed",
        "platform": platform.platform(),
        "machine": platform.machine(),
        "frozen": bool(getattr(sys, "frozen", False)),
    }

    def finish(error=None):
        if error:
            report["error"] = str(error)
        (output / "bundle-check.json").write_text(json.dumps(report, indent=2))
        window.close()
        app.exit(0 if report["status"] == "passed" else 1)

    def completed(result):
        try:
            assert len(result.project.bodies) == 2
            assert result.time[-1] >= 0.1 - 1e-10
            window.slider.setValue(window.slider.maximum())
            window.viewport.screenshot(output / "scene.png")
            report.update(
                status="passed",
                samples=len(result.time),
                manifest=json.loads(result.manifest_json),
                opengl=window.viewport.view.GetRenderWindow().ReportCapabilities(),
            )
            QTimer.singleShot(0, finish)
        except Exception as exc:
            finish(exc)

    sys.excepthook = lambda kind, value, tb: finish(value)
    window.controller.completed.connect(completed)
    window.controller.problem.connect(finish)
    window.show()
    QTimer.singleShot(0, window.run)
    QTimer.singleShot(60000, lambda: finish("Worker timeout"))
    return app.exec()
