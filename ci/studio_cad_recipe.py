"""Create a machined part through real CAD workers, save and capture Studio."""

import argparse
import json
import time
from pathlib import Path

from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio import __version__
from vinkulum_studio.cad_dialog import CadDialog
from vinkulum_studio.document import Body, replace_cad_body, save_project
from vinkulum_studio.editor import EditorWindow


def main(directory, width=1440, height=950):
    output = Path(directory)
    output.mkdir(parents=True, exist_ok=True)
    app = QApplication.instance() or QApplication([])
    app.setStyle("Fusion")
    window = EditorWindow()
    window.resize(width, height)
    window._discard_allowed = lambda: True
    window.show()
    window.new_project()
    timings = []

    def operation(kind, dimensions=(), position=(0, 0, 0), target=None, tool=None):
        dialog = CadDialog(window.project, target, window)
        dialog.operation.setCurrentIndex(dialog.operation.findData(kind))
        for field, value in zip(dialog.dimensions, dimensions):
            field.setValue(value)
        for field, value in zip(dialog.position, position):
            field.setValue(value)
        if tool:
            dialog.b.setCurrentIndex(dialog.b.findData(tool))
        dialog.radius.setValue(2)
        dialog.name.setText("Machined plate" if kind == "box" else "Hole tool")
        dialog.show()
        start = time.perf_counter()
        dialog.start()
        deadline = time.monotonic() + 60
        while dialog.process is not None and time.monotonic() < deadline:
            QTest.qWait(10)
        assert dialog.result_data is not None, dialog.status.text()
        timings.append(
            {
                "operation": kind,
                "wall_ms": 1000 * (time.perf_counter() - start),
                **dialog.result_data.get("timings_ms", {}),
            }
        )
        body = Body.from_dict(dialog.result_data["body"])
        window._commit(replace_cad_body(window.project, body), fit=True)
        window.select_object(body.id)
        dialog.deleteLater()
        return body.id

    plate = operation("box", (120, 75, 20))
    cutter = operation("cylinder", (12, 40), position=(25, 0, 0))
    operation("cut", target=plate, tool=cutter)
    window._commit(window.project.remove(cutter), fit=True)
    operation("fillet", target=plate)
    window.viewport.set_parallel_projection(True)
    window.select_object(plate)
    QTest.qWait(150)
    window.viewport.render()
    app.primaryScreen().grabWindow(window.winId()).save(str(output / "studio-cad.png"))
    window.viewport.screenshot(output / "cad-scene.png")
    save_project(output / "platine-percee.vinkulum.json", window.project)
    body = window.project.bodies[0]
    report = {
        "status": "passed",
        "app": __version__,
        "occt": body.cad.occt_version,
        "build123d": body.cad.build123d_version,
        "mass_kg": body.mass,
        "volume_m3": body.cad.volume_m3,
        "triangles": len(body.cad.triangles),
        "operations": timings,
    }
    (output / "recipe.json").write_text(json.dumps(report, indent=2) + "\n")
    window.close()
    print(json.dumps(report))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory")
    parser.add_argument("--width", type=int, default=1440)
    parser.add_argument("--height", type=int, default=950)
    args = parser.parse_args()
    main(args.directory, args.width, args.height)
