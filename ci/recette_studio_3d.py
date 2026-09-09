"""Capture an installed Studio and measure rendering; run on an actual Qt desktop.

python ci/recette_studio_3d.py --output /tmp/studio-3d
Linux automation: QT_QPA_PLATFORM=xcb xvfb-run -s '-screen 0 1920x1080x24' ...
No pass/fail FPS threshold: record the renderer and the measured latency.
"""

import argparse
from dataclasses import replace
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import tempfile
import time

import numpy as np
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

import vinkulum
import vinkulum._vinkulum as native
import vinkulum_studio
from vinkulum_studio.document import Body, History, Project, new_id
from vinkulum_studio.editor import EditorWindow
from vinkulum_studio.examples3d import double_pendulum, slider_crank
from vinkulum_studio.mechanism import MechanicalResult, simulate_project
from vinkulum_studio.viewport import Viewport


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])
    window = EditorWindow()
    window._discard_allowed = lambda: True
    report = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "python": platform.python_version(),
        "versions": {
            name: version(name)
            for name in ("vinkulum", "vinkulum-studio", "vtk", "PySide6", "numpy")
        },
        "modules": {"kernel": vinkulum.__file__, "studio": vinkulum_studio.__file__},
        "native_sha256": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
        "qt_platform": app.platformName(),
        "render_measurements": [],
        "human_desktop_recipe": "pending",
    }
    try:
        window.resize(1600, 1000)
        window.show()
        for name, factory in (
            ("double-pendulum", double_pendulum),
            ("slider-crank", slider_crank),
        ):
            project = replace(factory(), duration=0.5)
            with tempfile.TemporaryDirectory() as directory:
                data = simulate_project(project, name, directory)
                result = MechanicalResult.read(data, name, project, directory)
            window.history = History(project)
            window.selection = project.bodies[0].id
            window._refresh(True)
            window._completed(result)
            window.slider.setValue(len(result.time) // 2)
            window.viewport.camera("iso")
            QTest.qWait(100)
            app.primaryScreen().grabWindow(window.winId()).save(
                str(args.output / f"{name}.png")
            )
            window.viewport.screenshot(str(args.output / f"{name}-viewport.png"))
        report["opengl"] = window.viewport.view.GetRenderWindow().ReportCapabilities()
    finally:
        window.close()
    # Independent 1920x1080 viewport: includes actor updates and synchronous Render.
    view = Viewport()
    try:
        view.resize(1920, 1080)
        view.show()
        app.processEvents()
        for count in (3, 32):
            bodies = tuple(
                Body(new_id(), f"Corps {i}", position=(i % 8 * 0.3, i // 8 * 0.3, 0.0))
                for i in range(count)
            )
            project = Project(new_id(), bodies=bodies)
            view.set_project(project, fit=True)
            view.set_editable(False)
            elapsed = []
            for frame in range(140):
                start = time.perf_counter()
                poses = {
                    b.id: (
                        (b.position[0], b.position[1], 0.1 * np.sin(frame * 0.05 + i)),
                        b.orientation,
                    )
                    for i, b in enumerate(bodies)
                }
                view.set_poses(poses, frame / 60)
                app.processEvents()
                if frame >= 20:
                    elapsed.append((time.perf_counter() - start) * 1000)
            report["render_measurements"].append(
                {
                    "bodies": count,
                    "render_size": list(view.view.GetRenderWindow().GetSize()),
                    "samples": len(elapsed),
                    "median_ms": float(np.median(elapsed)),
                    "p95_ms": float(np.percentile(elapsed, 95)),
                    "mean_fps": float(1000 / np.mean(elapsed)),
                    "includes": "moving actor updates, synchronous render and Qt events; no joints or loads",
                }
            )
    finally:
        view.shutdown()
        view.close()
    (args.output / "desktop.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "versions": report["versions"],
                "render_measurements": report["render_measurements"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
