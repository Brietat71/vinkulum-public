"""Qualify the delivered CAD, mesh and operator workspaces in the frozen app."""

import hashlib
import os
import time
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QEventLoop, QTimer
from PySide6.QtWidgets import QFileDialog

from .cad_history_dialog import CadHistoryDialog
from .calculix import load_static_result, load_study
from .document import History, load_project
from .meshing import load_mesh


def delay(milliseconds=20):
    loop = QEventLoop()
    QTimer.singleShot(milliseconds, loop.quit)
    loop.exec()


def wait(predicate, diagnostic, seconds=180):
    deadline = time.monotonic() + seconds
    while not predicate() and time.monotonic() < deadline:
        delay()
    if not predicate():
        raise RuntimeError(diagnostic())


def check_previews(editor, output, examples):
    """Use shipped files and actual windows; never import an optional engine."""
    from .bundle_check import check_scene_image

    examples, output = Path(examples), Path(output)
    before = {
        p: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
        for p in examples.rglob("*")
        if p.is_file()
    }

    def capture(window, name):
        delay(150)
        if not window.screen().grabWindow(window.winId()).save(str(output / name)):
            raise RuntimeError(f"Could not capture {name}.")
        path = output / (Path(name).stem + "-scene.png")
        window.viewport.screenshot(path)
        return check_scene_image(path)

    # Captured operators open without a Pinocchio installation or worker.
    editor.open_articulated_study()
    articulated = editor._articulated_window
    with patch.object(
        QFileDialog,
        "getOpenFileName",
        return_value=(str(examples / "pinocchio/result.json"), "JSON"),
    ):
        articulated.open_result()
    wait(lambda: not articulated.archive.busy, articulated.status.text)
    if articulated.result is None:
        raise RuntimeError(articulated.status.text())
    report = {
        "pinocchio_capture": {
            "status": "passed",
            "coordinates": articulated.result.report["nv"],
            "engine_executed": False,
            "image_check": capture(articulated, "pinocchio-workspace.png"),
        }
    }
    articulated._discard_allowed = lambda: True
    if not articulated.close():
        raise RuntimeError("The idle articulated window did not close.")

    project = load_project(examples / "parametric-plate.vinkulum.json")
    body = project.bodies[0]
    editor.workspace_tabs.setCurrentIndex(0)
    editor.mode.setCurrentIndex(0)
    editor.history = History(project)
    editor._saved = project
    editor._refresh(fit=True)
    editor.select_object(body.id)
    dialog = CadHistoryDialog(project, body.id, editor)
    dialog.show()
    dialog.dimensions[0].setText("150")
    dialog.start_preview()
    wait(lambda: dialog.controller.process is None, dialog.status.text)
    if dialog.preview_body is None or not dialog.apply_button.isEnabled():
        raise RuntimeError(dialog.status.text())
    assert dialog.project == project and dialog.result_body is None
    assert dialog.preview_body.mass > body.mass
    assert dialog.preview_body.cad.recipe.features[0].dimensions_mm[0] == 150
    dialog.apply_button.click()
    assert dialog.result_body is not None and dialog.result_body.id == body.id
    report["cad_history"] = {
        "status": "passed",
        "stock_length_mm": 150,
        "original_mass_kg": body.mass,
        "regenerated_mass_kg": dialog.result_body.mass,
        "separate_cad_worker": True,
    }

    editor.commands["cad_study"].trigger()
    mesh = editor._mesh_window
    if mesh is None:
        raise RuntimeError(editor.status.text())
    executable = os.environ.get("VINKULUM_BUNDLE_GMSH")
    report["gmsh"] = {"status": "not_requested", "bundled_engine": False}
    if executable:
        mesh.executable.setText(executable)
        mesh.output.setText(str(output / "generated-mesh"))
        mesh.size.setText("40")
        mesh.generate.click()
        wait(lambda: not mesh.busy, mesh.status.text)
        if mesh.solid is None or mesh.mesh_root is None:
            raise RuntimeError(mesh.status.text())
        assert load_mesh(mesh.mesh_root)[0] == mesh.solid
        report["gmsh"] = {
            "status": "passed",
            "bundled_engine": False,
            "nodes": len(mesh.solid.mesh.nodes),
            "elements": len(mesh.solid.mesh.elements),
            "result_sha256": hashlib.sha256(
                (mesh.mesh_root / "result.json").read_bytes()
            ).hexdigest(),
            "gmsh_version": mesh.solid.gmsh_version,
            "occt_version": mesh.solid.occ_version,
            "executable_sha256": mesh.solid.executable_sha256,
            "archive_reopened": True,
        }
    mesh._discard_allowed = lambda: True
    path = examples / "cad-statics/calculation/study.json"
    with patch.object(QFileDialog, "getOpenFileName", return_value=(str(path), "JSON")):
        mesh.open_study()
        wait(lambda: not mesh.busy, mesh.status.text)
    saved = load_study(path)
    assert mesh.solid == saved.mesh_binding.solid
    assert mesh.conditions == saved.mesh_binding.conditions
    mesh.set_selected_faces((mesh.conditions[0].surfaces[0],))
    report["cad_mesh_capture"] = {
        "status": "passed",
        "nodes": len(mesh.solid.mesh.nodes),
        "elements": len(mesh.solid.mesh.elements),
        "face_conditions": len(mesh.conditions),
        "image_check": capture(mesh, "mesh-workspace.png"),
        "engine_executed": False,
    }
    mesh.open_calculix.click()
    wait(lambda: not mesh.busy, mesh.status.text)
    static = mesh._static_window
    if static is None or static.study != saved:
        raise RuntimeError("Captured CAD study changed during solver handoff.")
    result_path = examples / "cad-statics/calculation/result.json"
    with patch.object(
        QFileDialog, "getOpenFileName", return_value=(str(result_path), "JSON")
    ):
        static.open_result_dialog()
    wait(lambda: not static.archive.busy, static.status.text)
    if static.result is None:
        raise RuntimeError(static.status.text())
    assert static.result[1] == load_static_result(result_path)[1]
    static._auto_scale()
    report["quadratic_static_capture"] = {
        "status": "passed",
        "energy_J": static.result[1]["strain_energy_J"],
        "image_check": capture(static, "cad-static-workspace.png"),
        "engine_executed": False,
    }
    static._discard_allowed = lambda: True
    if not mesh.close():
        raise RuntimeError("The idle mesh workspace did not close.")
    assert before == {
        p: (hashlib.sha256(p.read_bytes()).hexdigest(), p.stat().st_mtime_ns)
        for p in before
    }, "The packaged examples were modified during inspection."
    report["example_files_unchanged"] = True
    return report
