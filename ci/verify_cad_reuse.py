#!/usr/bin/env python3
"""Recheck captured desktop CAD measurements without importing OCCT or Qt."""

import argparse
import hashlib
import json
import runpy
import statistics
import sys
from pathlib import Path

from vinkulum_studio.document import Body

REFERENCE = runpy.run_path(str(Path(__file__).with_name("qualify_cad_service.py")))


def verify(root):
    references, source_hashes, count = {}, None, 0
    for case in ("fresh-100", "cached-100", "cached-200", "fresh-200"):
        directory = root / case
        report = json.loads((directory / "qualification.json").read_text())
        mode, scale = case.split("-")
        assert report["mode"] == mode and report["dpi_scale"] == int(scale) / 100
        assert report["threads"] == 2
        if source_hashes is None:
            source_hashes = report["source_sha256"]
        assert report["source_sha256"] == source_hashes
        rows = report["rows"]
        expected = {
            (repeat, name) for repeat in range(7) for name in ("box", "cut", "fillet")
        }
        if mode == "cached":
            expected.add(("first-use", "box"))
        assert len(rows) == len(expected)
        assert {(row["repeat"], row["case"]) for row in rows} == expected
        for row in rows:
            assert row["retained_cad_dialogs"] == 0
            assert row["click_to_render_ms"] > 0 and row["after_worker_ms"] >= 0
            assert row["gui_max_rss_kib"] > 0 and row["worker_max_rss_kib"] > 0
            for field in ("body", "request"):
                path = directory / row[field]
                assert path.parent == directory
                assert (
                    hashlib.sha256(path.read_bytes()).hexdigest()
                    == row[field + "_sha256"]
                ), path
            request = json.loads((directory / row["request"]).read_text())
            assert request["operation"] == row["case"]
            assert request["execution"] == {"threads": 2, "budget": 2}
            assert request["density"] == 7800
            body = Body.from_dict(json.loads((directory / row["body"]).read_text()))
            assert (
                body.mass == row["mass_kg"] and body.cad.volume_m3 == row["volume_m3"]
            )
            assert list(body.inertia()) == row["inertia_kg_m2"]
            assert len(body.cad.triangles) == row["triangles"]
            REFERENCE["analytic"](row["case"], body)
            if row["case"] in references:
                REFERENCE["equivalent"](body, references[row["case"]])
            else:
                references[row["case"]] = body
            count += 1
        for name, values in report["summary"].items():
            times = [
                r["click_to_render_ms"]
                for r in rows
                if r["case"] == name and r["repeat"] != "first-use"
            ]
            assert values == {
                "median_ms": statistics.median(times),
                "min_ms": min(times),
                "max_ms": max(times),
            }
        if mode == "cached":
            assert len({r["worker_pid"] for r in rows}) == 1
            assert len({r["service_request_id"] for r in rows}) == len(rows)
            assert all(r["service_request_id"] for r in rows)
            assert report["idle_worker_rss_kib"] > 0
            assert report["idle_expiry_wait_ms"] > 0
        else:
            assert all(r["service_request_id"] is None for r in rows)
    assert count == 86
    assert not any(name in sys.modules for name in ("OCP", "build123d", "PySide6"))
    print(f"Verified {count} captured GUI CAD results without importing OCCT or Qt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    verify(parser.parse_args().root.resolve())
