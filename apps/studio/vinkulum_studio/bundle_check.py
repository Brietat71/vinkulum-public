"""Exercise the installed application and its own worker, including native GL."""

import json
import os
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
    faulthandler.dump_traceback_later(300, exit=True)
    output = Path(directory)
    output.mkdir(parents=True, exist_ok=True)
    # Exercise the same separate CAD worker the installed GUI launches.
    with tempfile.TemporaryDirectory(prefix="vinkulum-bundle-cad-") as cad_directory:
        request, response = (
            Path(cad_directory) / "request.json",
            Path(cad_directory) / "response.json",
        )
        from .engine_threads import process_environment

        request.write_text(
            json.dumps(
                {"operation": "self_check", "execution": {"threads": 1, "budget": 1}}
            )
        )
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
            env=process_environment(1, engine="occt"),
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
        "cad_execution": cad_report["execution"],
        "calculix": {"status": "not_requested"},
    }
    finished = False

    def finish(error=None):
        nonlocal finished
        if finished:
            return
        finished = True
        if error:
            report["status"] = "failed"
            report["error"] = str(error)
        (output / "bundle-check.json").write_text(json.dumps(report, indent=2))
        window.close()
        app.exit(0 if report["status"] == "passed" else 1)

    def static_completed(result):
        import numpy as np

        from .calculix import load_static_result

        try:
            study, values, directory = result
            np.testing.assert_allclose(
                values["strain_energy_J"], 1000**2 / (2 * 210e9 * 0.01), rtol=5e-7
            )
            np.testing.assert_allclose(
                np.array(values["stress"])[:, 0], 100000, rtol=5e-7
            )
            reopened = load_static_result(directory)
            assert reopened[0] == study and reopened[1] == values
            static = window._static_window
            static._auto_scale()
            static.viewport.screenshot(output / "static-scene.png")
            image_check = check_scene_image(output / "static-scene.png")
            report["calculix"] = {
                "status": "passed",
                "engine_version": values["engine_version"],
                "engine_sha256": values["engine_sha256"],
                "input_sha256": values["input_sha256"],
                "strain_energy_J": values["strain_energy_J"],
                "archive_reopened": True,
                "image_check": image_check,
                "bundled_engine": False,
            }
            report["status"] = "passed"
            QTimer.singleShot(0, finish)
        except Exception as exc:  # noqa: BLE001 -- report failures raised inside the Qt callback.
            finish(exc)

    def completed(result):
        try:
            assert len(result.project.bodies) == 2
            assert result.time[-1] >= 0.1 - 1e-10
            window.slider.setValue(window.slider.maximum())
            window.viewport.screenshot(output / "scene.png")
            image_check = check_scene_image(output / "scene.png")
            report.update(
                samples=len(result.time),
                image_check=image_check,
                manifest=json.loads(result.manifest_json),
                opengl=window.viewport.view.GetRenderWindow().ReportCapabilities(),
            )
            examples = (
                Path(sys.executable).parent / "Examples"
                if getattr(sys, "frozen", False) and sys.platform == "linux"
                else Path(os.environ["VINKULUM_BUNDLE_EXAMPLES"])
                if os.environ.get("VINKULUM_BUNDLE_EXAMPLES")
                else None
            )
            if examples is not None:
                if not examples.is_dir():
                    raise RuntimeError("The delivered Examples directory is missing.")
                from .bundle_previews import check_previews

                report["source_previews"] = check_previews(window, output, examples)
            if os.environ.get("VINKULUM_BUNDLE_CALCULIX") == "1":
                window.open_static_study()
                static = window._static_window
                static._discard_allowed = lambda: True
                static.output.setText(str(output / "static"))
                static.controller.completed.connect(static_completed)
                static.controller.problem.connect(finish)
                static.start()
                if static.controller.process is None:
                    raise RuntimeError(static.status.text())
            else:
                report["status"] = "passed"
                QTimer.singleShot(0, finish)
        except Exception as exc:  # noqa: BLE001 -- report failures raised inside the Qt callback.
            finish(exc)

    sys.excepthook = lambda kind, value, tb: finish(value)
    window.controller.completed.connect(completed)
    window.controller.problem.connect(finish)
    window.show()
    QTimer.singleShot(0, window.run)
    QTimer.singleShot(270000, lambda: finish("Bundle qualification timeout"))
    return app.exec()
