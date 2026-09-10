"""Capture actual Studio imports; the retained STEP inputs are never rewritten.

SPDX-License-Identifier: Apache-2.0
Requires the installed Studio CAD environment and a real Qt display.
"""

import hashlib
import json
import platform
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
from PySide6.QtCore import QSettings, QTimer, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QFileDialog
import vinkulum_studio
from vinkulum_studio.cad_dialog import CadDialog
from vinkulum_studio.document import save_project
from vinkulum_studio.editor import EditorWindow


def main(output):
    output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).parent
    reference = json.loads((root / "reference.json").read_text())
    app = QApplication.instance() or QApplication([])
    app.setStyle("Fusion")
    reports, errors = {}, []
    for name, expected in reference.items():
        source = root / f"{name}.step"
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        window = EditorWindow(
            QSettings(str(output / f"{name}.ini"), QSettings.Format.IniFormat)
        )
        window._discard_allowed = lambda: True
        window.resize(1440, 950)
        window.show()
        window.new_project()
        watchdog = QTimer(window)
        watchdog.setSingleShot(True)

        def drive():
            dialog = app.activeModalWidget()
            try:
                assert isinstance(dialog, CadDialog)
                dialog.name.setText(name)
                dialog.density.setValue(expected["density_kg_m3"])
                QTest.mouseClick(dialog.apply_button, Qt.MouseButton.LeftButton)
            except BaseException as error:
                errors.append(repr(error))
                if dialog is not None:
                    dialog.reject()

        def timed_out():
            dialog = app.activeModalWidget()
            errors.append(
                dialog.status.text()
                if isinstance(dialog, CadDialog)
                else "Import timed out"
            )
            if dialog is not None:
                dialog.reject()

        watchdog.timeout.connect(timed_out)
        try:
            watchdog.start(30000)
            QTimer.singleShot(0, drive)
            with patch.object(
                QFileDialog, "getOpenFileName", return_value=(str(source), "STEP")
            ):
                window.open_cad("import_step")
            watchdog.stop()
            assert not errors, errors
            body = window.object()
            assert body is not None and body.cad is not None, window.status.text()
            assert len(window.history.past) == 1
            np.testing.assert_allclose(body.mass, expected["mass_kg"], rtol=1e-8)
            np.testing.assert_allclose(
                body.inertia(),
                np.array(expected["inertia_kg_m2"]).ravel(),
                rtol=1e-8,
                atol=1e-14,
            )
            window.viewport.set_parallel_projection(True)
            q = np.array(expected["rotation"])
            camera = window.viewport.renderer.GetActiveCamera()
            camera.SetFocalPoint(*body.position)
            camera.SetPosition(*(np.array(body.position) + q @ (0.10, -0.13, 0.16)))
            camera.SetViewUp(*(q @ (0, 1, 0)))
            window.viewport.renderer.ResetCamera()
            window.viewport.render()
            QTest.qWait(100)
            assert (
                app.primaryScreen()
                .grabWindow(window.winId())
                .save(str(output / f"{name}.png"))
            )
            if name == "offset-closed-void-mm":
                # Diagnostic display only: show the enclosed shell through the
                # outer one. Neither the solid nor the product's default view changes.
                actor = window.viewport._actors[body.id][0]
                actor.GetProperty().SetRepresentationToWireframe()
                window.viewport.render()
                QTest.qWait(100)
                assert (
                    app.primaryScreen()
                    .grabWindow(window.winId())
                    .save(str(output / f"{name}-wireframe.png"))
                )
                actor.GetProperty().SetRepresentationToSurface()
            save_project(output / f"{name}.vinkulum.json", window.project)
            reports[name] = {
                "source_sha256": digest,
                "volume_m3": body.cad.volume_m3,
                "mass_kg": body.mass,
                "centre_m": body.position,
                "inertia_kg_m2": np.array(body.inertia()).reshape(3, 3).tolist(),
                "occt": body.cad.occt_version,
                "build123d": body.cad.build123d_version,
                "undo_transactions": len(window.history.past),
            }
            assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
        finally:
            watchdog.stop()
            window.close()
            window.deleteLater()
            QTest.qWait(30)
    report = {
        "studio": vinkulum_studio.__version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "parts": reports,
    }
    (output / "observed.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
