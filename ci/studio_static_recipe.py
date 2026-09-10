#!/usr/bin/env python3
"""Capture the real CalculiX workspace, including its preserved run directory."""

import argparse
import json
import time
from pathlib import Path

from PySide6.QtCore import QLocale
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio.calculix import load_study
from vinkulum_studio.static_window import StaticWindow


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    app = QApplication.instance() or QApplication([])
    QLocale.setDefault(QLocale("en_GB"))
    app.setStyle("Fusion")
    window = StaticWindow()
    window.show()
    try:
        example = (
            Path(__file__).resolve().parent.parent
            / "examples/studio/fem/tension.ccx.json"
        )
        window.set_study(load_study(example), "Tension specimen · 8-element patch")
        window.output.setText(str(output))
        QTest.qWait(100)
        window.start()
        deadline = time.monotonic() + 75
        while window.controller.process is not None and time.monotonic() < deadline:
            QTest.qWait(20)
        if window.result is None:
            raise RuntimeError(window.status.text())
        window._auto_scale()
        window.channels.setCurrentIndex(1)
        QTest.qWait(150)
        window.viewport.render()
        app.primaryScreen().grabWindow(window.winId()).save(
            str(output / "studio-static.png")
        )
        window.viewport.screenshot(output / "displacement-field.png")
        study, result, run_directory = window.result
        report = {
            "status": "passed",
            "engine": result["engine"],
            "engine_version": result["engine_version"],
            "nodes": len(study.nodes),
            "elements": len(study.elements),
            "strain_energy_J": result["strain_energy_J"],
            "scientific_status": result["scientific_status"],
            "display_deformation_multiplier": window.scale.value(),
            "input_sha256": result["input_sha256"],
            "run_directory": run_directory.name,
        }
        (output / "recipe.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    finally:
        window.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
