#!/usr/bin/env python3
"""Capture an actual articulated workspace and its separately computed operators."""

import argparse
import json
import time
from pathlib import Path

from PySide6.QtCore import QLocale, QSettings
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio import __version__
from vinkulum_studio.articulated_window import ArticulatedWindow
from vinkulum_studio.examples3d import double_pendulum


def main(output, interpreter):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    app = QApplication.instance() or QApplication([])
    QLocale.setDefault(QLocale("en_GB"))
    app.setStyle("Fusion")
    window = ArticulatedWindow(
        double_pendulum(),
        settings=QSettings(str(output / "settings.ini"), QSettings.Format.IniFormat),
    )
    window.interpreter.setText(str(interpreter.absolute()))
    window.output.setText(str(output / "runs"))
    state = {
        "q": [0.2, -0.3],
        "velocity": [0.1, 0.4],
        "acceleration": [0.7, 0.2],
        "effort": [0.9, 0.1],
        "time_s": 0.0,
    }
    # This is the recipe's starting input, so a failed run cannot cause a modal
    # unsaved-input question in its cleanup. Interactive editing is tested apart.
    window.set_project(window.project, state)
    window.show()
    try:
        QTest.qWait(100)
        window.start()
        deadline = time.monotonic() + 40
        while window.controller.process is not None and time.monotonic() < deadline:
            QTest.qWait(20)
        if window.result is None:
            raise RuntimeError(window.status.text())
        window.channel.setCurrentIndex(5)
        window.body_choice.setCurrentIndex(1)
        QTest.qWait(100)
        window.viewport.render()
        if (
            not app.primaryScreen()
            .grabWindow(window.winId())
            .save(str(output / "studio-pinocchio.png"))
        ):
            raise RuntimeError("Could not capture the articulated workspace.")
        window.viewport.screenshot(output / "mechanism.png")
        result = window.result
        report = {
            "status": "passed",
            "studio_version": __version__,
            "engine": result.report["engine"],
            "engine_version": result.report["engine_version"],
            "state": result.state,
            "nq": result.report["nq"],
            "nv": result.report["nv"],
            "scientific_status": result.report["scientific_status"],
            "project_sha256": result.report["project_sha256"],
            "run_directory": str(result.directory.relative_to(output)),
        }
        (output / "recipe.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    finally:
        window.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--python",
        type=Path,
        required=True,
        help="Separate Pinocchio environment Python",
    )
    options = parser.parse_args()
    main(options.output, options.python)
