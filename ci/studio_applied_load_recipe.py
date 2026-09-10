#!/usr/bin/env python3
"""Capture the actual applied-load derivative table from a retained result."""

import argparse
import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "xcb")

from PySide6.QtCore import QSettings
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from vinkulum_studio.articulated_result import load_operators
from vinkulum_studio.articulated_window import ArticulatedWindow


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    app = QApplication([])
    result = load_operators(args.result)
    window = ArticulatedWindow(
        result.project,
        settings=QSettings(
            str(args.output / "settings.ini"), QSettings.Format.IniFormat
        ),
    )
    window.resize(1440, 940)
    window.set_project(result.project, result.state)
    window.show()
    window._completed(result)
    window.channel.setCurrentIndex(6)
    QTest.qWait(350)
    capture = app.primaryScreen().grabWindow(int(window.winId()))
    if capture.isNull() or not capture.save(
        str(args.output / "external-derivative.png")
    ):
        raise RuntimeError("Could not capture the native articulated window.")
    window.close()
    QTest.qWait(20)


if __name__ == "__main__":
    main()
