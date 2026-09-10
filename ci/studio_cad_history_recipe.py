#!/usr/bin/env python3
"""Create a real plate through the CAD worker, edit an upstream feature and capture."""

import argparse
import json
import time
from dataclasses import asdict
from pathlib import Path

from PySide6.QtCore import QLocale, QSettings
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio import __version__
from vinkulum_studio.cad_controller import CadController
from vinkulum_studio.cad_history_dialog import CadHistoryDialog
from vinkulum_studio.document import Body, Project, new_id, save_project
from vinkulum_studio.editor import EditorWindow


def wait(controller):
    deadline = time.monotonic() + 65
    while controller.process is not None and time.monotonic() < deadline:
        QTest.qWait(20)
    if controller.process is not None or controller.failure_kind:
        raise RuntimeError(str(controller.failure_kind) + "\n" + controller.last_log)


def main(output):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    app = QApplication.instance() or QApplication([])
    QLocale.setDefault(QLocale("en_GB"))
    app.setStyle("Fusion")
    controller = CadController()

    def run(**request):
        controller.start(request)
        wait(controller)
        return Body.from_dict(controller.last_result["body"])

    plate = run(
        operation="box",
        name="Machined plate",
        dimensions_mm=(120, 75, 20),
        density=7800,
    )
    tool = run(
        operation="cylinder",
        name="Hole cutter",
        dimensions_mm=(12, 40),
        position_mm=(30, 5, 0),
        density=7800,
    )
    plate = run(operation="cut", a=asdict(plate), b=asdict(tool), density=7800)
    plate = run(operation="fillet", a=asdict(plate), radius_mm=2, density=7800)
    project = Project(new_id(), "Parametric machined plate", bodies=(plate,))
    editor = EditorWindow(
        QSettings(str(output / "settings.ini"), QSettings.Format.IniFormat)
    )
    editor.show()
    dialog = CadHistoryDialog(project, plate.id, editor)
    dialog.show()
    try:
        QTest.qWait(100)
        field = dialog.dimensions[0]
        field.setFocus()
        field.selectAll()
        QTest.keyClicks(field, "150")
        dialog.start_preview()
        wait(dialog.controller)
        if not dialog.apply_button.isEnabled():
            raise RuntimeError(dialog.status.text())
        QTest.qWait(100)
        dialog.viewport.render()
        if (
            not app.primaryScreen()
            .grabWindow(dialog.winId())
            .save(str(output / "cad-history.png"))
        ):
            raise RuntimeError("Could not capture the CAD feature editor.")
        dialog.viewport.screenshot(output / "plate.png")
        save_project(output / "before.vinkulum.json", project)
        save_project(
            output / "after.vinkulum.json", project.replace_object(dialog.preview_body)
        )
        report = {
            "status": "passed",
            "studio_version": __version__,
            "occt_version": dialog.preview_body.cad.occt_version,
            "edited_parameter": "Stock length: 120 mm → 150 mm",
            "features": len(dialog.preview_body.cad.recipe.features),
            "mass_before_kg": plate.mass,
            "mass_after_kg": dialog.preview_body.mass,
            "recipe_sha256": dialog.preview_body.cad.recipe.fingerprint(),
        }
        (output / "recipe.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    finally:
        dialog.reject()
        controller.shutdown()
        editor.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
