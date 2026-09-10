"""Capture signed pressure, total force and support symbols from a saved CAD study.

Run with an installed Studio wheel under a real Qt/OpenGL desktop or Xvfb.
The input must contain one pressure condition and at least one support. No
external engine is executed. Output images are unmodified window captures.
"""

import argparse
import hashlib
import json
import time
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog
from vinkulum_studio import __version__
from vinkulum_studio.boundary_display import effective_values
from vinkulum_studio.calculix import load_study
from vinkulum_studio.mesh_window import MeshWindow


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(study_path, output):
    study_path, output = study_path.resolve(), output.resolve()
    study = load_study(study_path)
    binding = study.mesh_binding
    if binding is None:
        raise ValueError("Open a captured CAD-face study.")
    pressure = [c for c in binding.conditions if c.kind == "pressure"]
    supports = tuple(c for c in binding.conditions if c.kind == "support")
    if len(pressure) != 1 or not supports:
        raise ValueError("This recipe needs one pressure condition and a support.")
    before = digest(study_path)
    output.mkdir(parents=True, exist_ok=False)
    app = QApplication.instance() or QApplication([])
    window = MeshWindow(binding.solid.request.body)
    window.show()
    with patch.object(
        QFileDialog, "getOpenFileName", return_value=(str(study_path), "JSON")
    ):
        window.open_study()
    deadline = time.monotonic() + 180
    while window.busy and time.monotonic() < deadline:
        QTest.qWait(20)
    if window.busy or window.conditions != binding.conditions:
        raise RuntimeError(window.status.text())
    window.edges.setChecked(False)
    camera = window.viewport.renderer.GetActiveCamera()
    camera.SetPosition(-1, -1, 0.8)
    camera.SetFocalPoint(0, 0, 0)
    window.viewport.fit_scene()
    force = replace(
        pressure[0],
        kind="total_force",
        name="Upward resultant",
        values=(0.0, 0.0, 120.0),
    )
    cases = []
    for name, factor, conditions in (
        ("inward-pressure", 1.0, supports + tuple(pressure)),
        ("outward-pressure", -1.0, supports + tuple(pressure)),
        ("world-force", 1.0, supports + (force,)),
    ):
        window.load_factor.setText(repr(factor))
        window._set_conditions(conditions)
        window.set_selected_faces(())
        window.viewport.render()
        QTest.qWait(150)
        path = output / f"{name}.png"
        if not window.screen().grabWindow(window.winId()).save(str(path)):
            raise RuntimeError("Could not capture the rendered workspace.")
        cases.append(
            {
                "name": name,
                "load_factor": factor,
                "conditions": [
                    {**asdict(c), "applied_values": effective_values(c, factor)}
                    for c in conditions
                ],
                "glyphs": [asdict(g) for g in window.viewport.glyph_layer.glyphs],
                "image": path.name,
                "image_sha256": digest(path),
            }
        )
    window.load_factor.setText(repr(binding.load_factor))
    window._set_conditions(binding.conditions)
    if digest(study_path) != before:
        raise RuntimeError("The saved study changed during display qualification.")
    report = {
        "studio_version": __version__,
        "study_sha256": before,
        "input_unchanged": True,
        "nodes": len(binding.solid.mesh.nodes),
        "elements": len(binding.solid.mesh.elements),
        "device_pixel_ratio": window.devicePixelRatioF(),
        "cases": cases,
        "scope": "Sampled directions on the initial captured mesh; not nodal-force or stress validation.",
    }
    (output / "recipe.json").write_text(json.dumps(report, indent=2) + "\n")
    if not window.close():
        raise RuntimeError("The restored study did not close cleanly.")
    app.processEvents()
    print(json.dumps({k: v for k, v in report.items() if k != "cases"}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--study", required=True, type=Path)
    args = parser.parse_args()
    main(args.study, args.output)
