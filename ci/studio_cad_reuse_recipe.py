#!/usr/bin/env python3
"""Measure the actual CAD dialog → document → render path and its process memory."""

import argparse
import functools
import hashlib
import json
import os
import platform
import runpy
import statistics
import sys
import time
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path
from unittest.mock import patch

from PySide6.QtCore import QCoreApplication, QEvent, QSettings, Qt, QTimer
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

import vinkulum_studio
from vinkulum_studio.cad_dialog import CadDialog
from vinkulum_studio.cad_process_cache import application_cad_cache
from vinkulum_studio.cpu_scheduler import application_scheduler
from vinkulum_studio.document import Body, History, Project, new_id, save_project
from vinkulum_studio.editor import EditorWindow
from vinkulum_studio.model import write_json

REFERENCE = runpy.run_path(str(Path(__file__).with_name("qualify_cad_service.py")))


def rss(pid):
    try:
        lines = Path(f"/proc/{pid}/status").read_text().splitlines()
        # A process can exit between QProcess.processId() and this read. Linux
        # retains a short-lived status entry without VmRSS until it is reaped.
        return next(
            (int(line.split()[1]) for line in lines if line.startswith("VmRSS:")), 0
        )
    except FileNotFoundError, ProcessLookupError:
        return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("fresh", "cached"), required=True)
    parser.add_argument("--repeats", type=int, default=7)
    args = parser.parse_args()
    if sys.platform != "linux" or not 1 <= args.repeats <= 20:
        parser.error("Requires Linux and 1–20 repetitions")
    os.environ["VINKULUM_STUDIO_CPUS"] = "2"
    os.environ["VINKULUM_CAD_IDLE_SECONDS"] = "15"
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    app = QApplication([])
    app.setStyle("Fusion")
    cache, scheduler = application_cad_cache(), application_scheduler()
    window = EditorWindow(
        QSettings(str(root / "settings.ini"), QSettings.Format.IniFormat)
    )
    window._discard_allowed = lambda: True
    window.resize(1440, 950)
    window.show()
    QTest.qWait(150)
    requests = REFERENCE["fixtures"]()
    a, b = (Body.from_dict(requests["cut"][key]) for key in ("a", "b"))
    base = Project(new_id(), "CAD worker reuse", bodies=(a, b))
    references, rows, errors = {}, [], []
    report = {
        "studio": vinkulum_studio.__version__,
        "platform": platform.platform(),
        "python": sys.version,
        "packages": {
            name: version(name)
            for name in (
                "vinkulum",
                "PySide6",
                "vtk",
                "build123d",
                "cadquery-ocp-novtk",
            )
        },
        "mode": args.mode,
        "dpi_scale": window.devicePixelRatioF(),
        "threads": scheduler.capacity,
        "gui_initial_rss_kib": rss(os.getpid()),
        "rows": rows,
        "source_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(Path(vinkulum_studio.__file__).parent.glob("*.py"))
        },
        "recipe_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reference_sha256": hashlib.sha256(
            Path(__file__).with_name("qualify_cad_service.py").read_bytes()
        ).hexdigest(),
    }

    def operation(name, repeat):
        window.history = History(base)
        window.selection = a.id
        window._refresh(True)
        row = {
            "case": name,
            "repeat": repeat,
            "gui_max_rss_kib": 0,
            "worker_max_rss_kib": 0,
        }
        running = {"dialog": None}
        timer, watchdog = QTimer(), QTimer()
        timer.setInterval(20)
        watchdog.setSingleShot(True)

        def sample():
            row["gui_max_rss_kib"] = max(row["gui_max_rss_kib"], rss(os.getpid()))
            dialog = running["dialog"]
            process = dialog.controller.process if dialog is not None else None
            if process is None:
                process = cache.idle_process
            if process is not None and process.processId():
                row["worker_pid"] = process.processId()
                row["worker_max_rss_kib"] = max(
                    row["worker_max_rss_kib"], rss(process.processId())
                )

        def timed_out():
            errors.append("CAD dialog qualification timed out")
            if running["dialog"] is not None:
                running["dialog"].reject()

        def drive():
            dialog = app.activeModalWidget()
            running["dialog"] = dialog
            try:
                assert isinstance(dialog, CadDialog)
                for field, value in zip(dialog.dimensions, (100, 60, 20)):
                    field.setValue(value)
                dialog.density.setValue(7800)
                dialog.radius.setValue(1)
                dialog.b.setCurrentIndex(dialog.b.findData(b.id))
                if args.mode == "fresh":
                    # Public interpreter override selects the retained one-shot
                    # path, with the same interpreter, geometry code and inputs.
                    dialog.controller.start = functools.partial(
                        dialog.controller.start, executable=sys.executable
                    )
                dialog.controller.completed.connect(
                    lambda result: row.update(
                        result=result,
                        request_data=dialog.controller.last_request,
                        worker_completed_at=time.perf_counter(),
                    )
                )
                row["clicked_at"] = time.perf_counter()
                timer.start()
                sample()
                QTest.mouseClick(dialog.apply_button, Qt.MouseButton.LeftButton)
            except BaseException as error:
                errors.append(repr(error))
                if isinstance(dialog, CadDialog):
                    dialog.reject()

        timer.timeout.connect(sample)
        watchdog.timeout.connect(timed_out)
        watchdog.start(60000)
        QTimer.singleShot(0, drive)
        window.open_cad(name)
        window.viewport.render()
        completed = time.perf_counter()
        sample()
        timer.stop()
        watchdog.stop()
        assert not errors, errors
        assert "result" in row, window.status.text()
        result = row.pop("result")
        row["request"] = f"{repeat}-{name}-request.json"
        write_json(root / row["request"], row.pop("request_data"))
        row["request_sha256"] = hashlib.sha256(
            (root / row["request"]).read_bytes()
        ).hexdigest()
        assert row["request_sha256"] == result["request_sha256"]
        row["click_to_render_ms"] = 1000 * (completed - row.pop("clicked_at"))
        row["after_worker_ms"] = 1000 * (completed - row.pop("worker_completed_at"))
        row["worker_timings_ms"] = result["timings_ms"]
        row["service_request_id"] = result.get("service_request_id")
        assert len(window.history.past) == 1
        assert scheduler.allocated == 0
        body = window.object()
        assert isinstance(body, Body) and body.cad is not None
        REFERENCE["analytic"](name, body)
        if name in references:
            REFERENCE["equivalent"](body, references[name])
        else:
            references[name] = body
        row["mass_kg"] = body.mass
        row["volume_m3"] = body.cad.volume_m3
        row["inertia_kg_m2"] = body.inertia()
        row["triangles"] = len(body.cad.triangles)
        row["body"] = f"{repeat}-{name}.json"
        (root / row["body"]).write_text(json.dumps(asdict(body), allow_nan=False))
        row["body_sha256"] = hashlib.sha256(
            (root / row["body"]).read_bytes()
        ).hexdigest()
        committed = window.project
        window.undo()
        assert window.project.bodies == base.bodies
        window.redo()
        assert window.project.bodies == committed.bodies
        QTest.qWait(20)
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        row["retained_cad_dialogs"] = len(window.findChildren(CadDialog))
        assert row["retained_cad_dialogs"] == 0
        rows.append(row)

    try:
        with patch(
            "sys.excepthook",
            side_effect=lambda kind, error, trace: errors.append(repr(error)),
        ):
            if args.mode == "cached":
                operation("box", "first-use")
            for repeat in range(args.repeats):
                for name in requests:
                    operation(name, repeat)
            assert not errors, errors
            measured = [row for row in rows if row["repeat"] != "first-use"]
            if args.mode == "cached":
                assert len({row["worker_pid"] for row in rows}) == 1
                assert len({row["service_request_id"] for row in rows}) == len(rows)
            else:
                assert all(row["service_request_id"] is None for row in rows)
            report["summary"] = {
                name: {
                    "median_ms": statistics.median(
                        row["click_to_render_ms"]
                        for row in measured
                        if row["case"] == name
                    ),
                    "min_ms": min(
                        row["click_to_render_ms"]
                        for row in measured
                        if row["case"] == name
                    ),
                    "max_ms": max(
                        row["click_to_render_ms"]
                        for row in measured
                        if row["case"] == name
                    ),
                }
                for name in requests
            }
            save_project(root / "last-model.vinkulum.json", window.project)
            capture = app.primaryScreen().grabWindow(int(window.winId()))
            assert not capture.isNull() and capture.save(str(root / "studio.png"))
            report["gui_final_rss_kib"] = rss(os.getpid())
            if args.mode == "cached":
                pid = cache.idle_process.processId()
                report["idle_worker_rss_kib"] = rss(pid)
                report["idle_remaining_ms"] = cache.timer.remainingTime()
                started = time.perf_counter()
                deadline = time.monotonic() + 20
                while (
                    cache.idle_process is not None or cache._retiring
                ) and time.monotonic() < deadline:
                    QTest.qWait(10)
                assert cache.idle_process is None and not cache._retiring
                assert rss(pid) == 0 and scheduler.allocated == 0
                report["idle_expiry_wait_ms"] = 1000 * (time.perf_counter() - started)
                report["gui_after_expiry_rss_kib"] = rss(os.getpid())
            (root / "qualification.json").write_text(
                json.dumps(report, indent=2, allow_nan=False) + "\n"
            )
            print(json.dumps(report["summary"], indent=2))
    finally:
        cache.shutdown()
        window.close()
        QTest.qWait(20)


if __name__ == "__main__":
    main()
