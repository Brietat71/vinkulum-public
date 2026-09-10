#!/usr/bin/env python3
"""Capture the real quadratic-tetrahedron workspace and independent bending study."""

import argparse
import hashlib
import json
import time
from pathlib import Path

from PySide6.QtCore import QLocale
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio import __version__
from vinkulum_studio.calculix import load_static_result, load_study
from vinkulum_studio.static_window import StaticWindow


def main(output):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    example = (
        Path(__file__).resolve().parents[1]
        / "examples/studio/fem/pure-bending-c3d10.ccx.json"
    )
    app = QApplication.instance() or QApplication([])
    QLocale.setDefault(QLocale("en_GB"))
    app.setStyle("Fusion")
    window = StaticWindow()
    window.show()
    try:
        window.set_study(load_study(example), "Pure bending · quadratic tetrahedra")
        window.output.setText(str(output))
        window.start()
        deadline = time.monotonic() + 65
        while window.controller.process is not None and time.monotonic() < deadline:
            QTest.qWait(20)
        if window.controller.process is not None or window.result is None:
            raise RuntimeError(window.status.text())
        study, result, root = window.result
        window._auto_scale()
        window.channels.setCurrentIndex(1)
        QTest.qWait(100)
        window.viewport.render()
        if (
            not app.primaryScreen()
            .grabWindow(window.winId())
            .save(str(output / "tetra-bending.png"))
        ):
            raise RuntimeError("Could not capture the tetrahedral workspace.")
        reopened, checked, _ = load_static_result(root)
        if reopened != study or checked != result:
            raise RuntimeError("The captured calculation could not be rechecked.")
        reference = 1e12 * (0.3 * 0.2**3 / 12) / (2 * 210e9)
        relative = abs(result["strain_energy_J"] / reference - 1)
        if relative > 6e-7:
            raise RuntimeError("The quadratic bending reference failed.")
        report = {
            "status": "passed",
            "studio_version": __version__,
            "engine_version": result["engine_version"],
            "element_type": study.element_type,
            "elements": len(study.elements),
            "nodes": len(study.nodes),
            "integration_points_per_element": 4,
            "calculation": root.name,
            "study_sha256": hashlib.sha256(example.read_bytes()).hexdigest(),
            "strain_energy_J": result["strain_energy_J"],
            "reference_energy_J": reference,
            "relative_energy_error": relative,
            "display_deformation_multiplier": window.scale.value(),
            "scientific_status": result["scientific_status"],
        }
        (output / "recipe.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    finally:
        window.controller.shutdown()
        window.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
