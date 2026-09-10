"""Check the actual launcher window after normal application initialisation."""

import json
import sys
from pathlib import Path

from PySide6.QtCore import QLocale, QTimer

from . import __version__
from .bundle_check import check_scene_image


def schedule_check(app, window, directory):
    output = Path(directory)
    output.mkdir(parents=True, exist_ok=True)
    finished = False

    def check():
        nonlocal finished
        if finished:
            return
        finished = True
        report = {
            "status": "failed",
            "app_version": __version__,
            "frozen": bool(getattr(sys, "frozen", False)),
        }
        try:
            if not window.isVisible() or not window.windowHandle().isExposed():
                raise RuntimeError("Normal launcher window was not exposed.")
            if QLocale().name() != "en_GB":
                raise RuntimeError(
                    "Normal launcher did not establish its English locale."
                )
            window.viewport.screenshot(output / "startup-scene.png")
            image_check = check_scene_image(output / "startup-scene.png")
            report.update(
                status="passed",
                locale=QLocale().name(),
                title=window.windowTitle(),
                image_check=image_check,
            )
        except Exception as error:  # noqa: BLE001 -- serialize errors raised in a Qt callback.
            report["error"] = str(error)
        (output / "startup-check.json").write_text(json.dumps(report, indent=2) + "\n")
        window.close()
        app.exit(0 if report["status"] == "passed" else 1)

    QTimer.singleShot(250, check)
