"""Import the retained STEP through the real CAD worker and capture Studio.

SPDX-License-Identifier: Apache-2.0
Usage: python capture.py OUTPUT_DIRECTORY (requires the Studio CAD environment).
"""

import hashlib
import json
import platform
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from vinkulum_studio import __version__
from vinkulum_studio.cad_dialog import CadDialog
from vinkulum_studio.document import Body, replace_cad_body
from vinkulum_studio.editor import EditorWindow


def main(output):
    output.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_name("l-bracket.step")
    app = QApplication.instance() or QApplication([])
    app.setStyle("Fusion")
    with tempfile.TemporaryDirectory() as settings:
        QSettings.setDefaultFormat(QSettings.Format.IniFormat)
        QSettings.setPath(
            QSettings.Format.IniFormat, QSettings.Scope.UserScope, settings
        )
        window = EditorWindow()
        window._discard_allowed = lambda: True
        window.resize(1440, 950)
        window.show()
        window.new_project()
        dialog = CadDialog(window.project, None, window)
        try:
            dialog.operation.setCurrentIndex(dialog.operation.findData("import_step"))
            dialog.name.setText("Independent L bracket")
            dialog.density.setValue(7800)
            dialog.show()
            # Supply the known fixture path in place of the interactive file chooser.
            dialog.controller.start({**dialog.request(), "path": str(source)})
            deadline = time.monotonic() + 60
            while dialog.process is not None and time.monotonic() < deadline:
                QTest.qWait(10)
            if dialog.result_data is None:
                raise RuntimeError(dialog.status.text() or "CAD import timed out")
            body = Body.from_dict(dialog.result_data["body"])
            window._commit(replace_cad_body(window.project, body), fit=True)
            window.select_object(body.id)
            window.viewport.set_parallel_projection(True)
            q = np.array(((-3, -4, 12), (12, 3, 4), (-4, 12, 3))) / 13
            camera = window.viewport.renderer.GetActiveCamera()
            camera.SetFocalPoint(*body.position)
            camera.SetPosition(*(np.array(body.position) + q @ (0.10, -0.13, 0.16)))
            camera.SetViewUp(*(q @ (0, 1, 0)))
            window.viewport.renderer.ResetCamera()
            window.viewport.render()
            QTest.qWait(200)
            if (
                not app.primaryScreen()
                .grabWindow(window.winId())
                .save(str(output / "studio.png"))
            ):
                raise RuntimeError("Screenshot could not be saved")
            report = {
                "status": "imported",
                "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "python": platform.python_version(),
                "studio": __version__,
                "occt": body.cad.occt_version,
                "build123d": body.cad.build123d_version,
                "volume_m3": body.cad.volume_m3,
                "mass_kg": body.mass,
                "centre_m": body.position,
                "inertia_kg_m2": np.array(body.inertia()).reshape(3, 3).tolist(),
            }
            (output / "observed.json").write_text(json.dumps(report, indent=2) + "\n")
            print(json.dumps(report))
        finally:
            dialog.reject()
            window.close()


if __name__ == "__main__":
    main(Path(sys.argv[1]))
