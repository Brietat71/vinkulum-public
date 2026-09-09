"""Capture and measure the real Studio workspace under a desktop Qt session.

Example: xvfb-run -a -s '-screen 0 3200x2200x24' env QT_SCALE_FACTOR=2 \
    python ci/studio_gui_recipe.py /tmp/studio-hidpi
"""

import json
import platform
import statistics
import sys
import time
from pathlib import Path

from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio.editor import EditorWindow
from vinkulum_studio.examples3d import slider_crank


def main(directory):
    output = Path(directory)
    output.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])
    app.setStyle("Fusion")
    window = EditorWindow()
    window._discard_allowed = lambda: True
    unhandled = []
    old_hook = sys.excepthook
    sys.excepthook = lambda kind, value, traceback: unhandled.append(str(value))
    report = {
        "platform": platform.platform(),
        "scale": app.primaryScreen().devicePixelRatio(),
    }

    def capture(name):
        QTest.qWait(80)
        window.viewport.render()
        pixmap = app.primaryScreen().grabWindow(window.winId())
        assert not pixmap.isNull()
        assert pixmap.save(str(output / (name + ".png")))
        report[name] = {
            "window_logical": list(window.size().toTuple()),
            "image_pixels": [pixmap.width(), pixmap.height()],
            "scene_logical": list(window.viewport.size().toTuple()),
        }

    def calculate():
        window.run()
        deadline = time.monotonic() + 20
        while window.controller.process is not None and time.monotonic() < deadline:
            QTest.qWait(10)
        assert window.controller.process is None, "Worker timed out"
        assert window.result is not None, window.status.text()
        return window.result

    try:
        window.show()
        window._commit(slider_crank(), fit=True)
        window.select_object(window.project.bodies[0].id)
        capture("modeling-dark")
        window.workspace_tabs.setCurrentIndex(1)
        capture("simulation-dark")
        first = calculate()
        window.mode.setCurrentIndex(0)
        window.step.setText("0.0025")
        second = calculate()
        assert second.run_id != first.run_id
        window.series_combo.setCurrentIndex(window.series_combo.count() - 1)
        window.compare_combo.setCurrentIndex(
            window.compare_combo.findData(first.run_id)
        )
        window.slider.setValue(len(second.time) // 2)
        capture("comparison-dark")
        # A laptop window must keep every playback control inside its panel.
        window.resize(1280, 800)
        window.set_theme("light")
        capture("comparison-light-1280")
        for control in (window.play_button, window.slider, window.time_input):
            point = control.mapTo(window, control.rect().bottomRight())
            assert window.rect().contains(point), (
                f"Clipped control: {control.objectName()}"
            )
        report["run_ids"] = [first.run_id, second.run_id]
        report["native_samples"] = [len(first.time), len(second.time)]
        window.resize(1440, 950)
        QTest.qWait(100)
        durations = []
        for index in range(60):
            start = time.perf_counter()
            window.slider.setValue(index)
            app.processEvents()
            durations.append(1000 * (time.perf_counter() - start))
        report["frame_ms_3_bodies_2_curves"] = {
            "median": statistics.median(durations),
            "p95": sorted(durations)[56],
        }
        report["renderer"] = window.viewport.view.GetRenderWindow().ReportCapabilities()
        window.close()
        QTest.qWait(20)
        assert not unhandled, unhandled
        report["status"] = "passed"
        (output / "recipe.json").write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n"
        )
        print(
            json.dumps(
                {k: v for k, v in report.items() if k != "renderer"}, ensure_ascii=False
            )
        )
    finally:
        window.close()
        sys.excepthook = old_hook


if __name__ == "__main__":
    main(sys.argv[1])
