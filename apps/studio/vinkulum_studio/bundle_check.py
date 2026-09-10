"""Exercise the installed application and its own worker, including native GL."""

import json
import platform
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication

from .document import History
from .editor import EditorWindow
from .examples3d import double_pendulum


def check_scene_image(path):
    """Detect an empty main renderer even when the orientation overlay still draws.

    This fixed recipe has a non-black graphite background. Probe its interior,
    excluding the corner overlay; this is a render health check, not image parity.
    """
    image = QImage(str(path))
    if image.isNull() or min(image.width(), image.height()) < 100:
        raise RuntimeError("Missing or undersized scene capture")
    visible = sum(
        max(
            image.pixelColor(
                x * image.width() // 12, y * image.height() // 12
            ).getRgb()[:3]
        )
        > 8
        for x in range(1, 11)
        for y in range(1, 11)
    )
    if visible < 80:
        raise RuntimeError("Main 3D renderer is black despite a live Qt window")
    return {"nonblack_probes": visible, "total_probes": 100}


def main(directory, *, app=None, window=None):
    import faulthandler

    # Covers Cocoa initialization too, before the Qt timer can run.
    faulthandler.dump_traceback_later(90, exit=True)
    output = Path(directory)
    output.mkdir(parents=True, exist_ok=True)
    # Exercise the same separate CAD worker the installed GUI launches.
    with tempfile.TemporaryDirectory(prefix="vinkulum-bundle-cad-") as cad_directory:
        request, response = (
            Path(cad_directory) / "request.json",
            Path(cad_directory) / "response.json",
        )
        request.write_text(json.dumps({"operation": "self_check"}))
        arguments = (
            ["--cad-worker"]
            if getattr(sys, "frozen", False)
            else ["-m", "vinkulum_studio.cad_worker"]
        )
        completed = subprocess.run(
            [sys.executable, *arguments, str(request), str(response)],
            capture_output=True,
            timeout=45,
            check=False,
        )
        if completed.returncode != 0 or not response.exists():
            raise RuntimeError(
                "The bundled CAD worker failed: "
                + completed.stderr.decode(errors="replace")[-2048:]
            )
        cad_report = json.loads(response.read_text())
        if cad_report.get("status") != "completed":
            raise RuntimeError(str(cad_report))
    app = app or QApplication.instance() or QApplication([])
    window = window or EditorWindow()
    window._discard_allowed = lambda: True
    project = replace(double_pendulum(), duration=0.1)
    window.history = History(project)
    window._refresh(True)
    report = {
        "status": "failed",
        "platform": platform.platform(),
        "machine": platform.machine(),
        "frozen": bool(getattr(sys, "frozen", False)),
        "cad": cad_report["cad_check"],
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
            image_check = check_scene_image(output / "scene.png")
            report.update(
                status="passed",
                samples=len(result.time),
                image_check=image_check,
                manifest=json.loads(result.manifest_json),
                opengl=window.viewport.view.GetRenderWindow().ReportCapabilities(),
            )
            QTimer.singleShot(0, finish)
        except Exception as exc:  # noqa: BLE001 -- report failures raised inside the Qt callback.
            finish(exc)

    sys.excepthook = lambda kind, value, tb: finish(value)
    window.controller.completed.connect(completed)
    window.controller.problem.connect(finish)
    window.show()
    QTimer.singleShot(0, window.run)
    QTimer.singleShot(60000, lambda: finish("Worker timeout"))
    return app.exec()
