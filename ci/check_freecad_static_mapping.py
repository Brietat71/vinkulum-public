"""Recheck captured native runs and challenge boundary identity without a solver."""

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import build123d as bd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apps/freecad"))
from static_worker import (
    binding_for,
    import_body,
    import_faces,
    load_mesh,
    load_static_result,
    match_faces,
)


def refused(action, message):
    try:
        action()
    except ValueError as error:
        assert message in str(error), str(error)
    else:
        raise AssertionError("Invalid face mapping accepted")


def check(directory):
    request = json.loads((directory / "request.json").read_text())
    body = import_body(directory, request)
    faces = import_faces(directory, request, float(np.linalg.norm(body.dimensions)))
    solid, _ = load_mesh(directory / "mesh")
    load_static_result(directory / "static")  # Recheck raw deck and DAT evidence.
    mapping, _ = match_faces(solid, faces, request["faces"])
    renumbered = replace(
        solid,
        surfaces=tuple(replace(s, id=100000 + 37 * s.id) for s in solid.surfaces),
    )
    changed, _ = match_faces(renumbered, faces, request["faces"])
    assert changed == {i: [100000 + 37 * n for n in ids] for i, ids in mapping.items()}
    assert (
        binding_for(solid, mapping, request).arrays
        == binding_for(renumbered, changed, request).arrays
    )
    refused(
        lambda: match_faces(solid, [faces[0], faces[0]], [request["faces"][0]] * 2),
        "ambiguous",
    )
    refused(
        lambda: match_faces(
            solid,
            [faces[0].moved(bd.Location((1000, 1000, 1000)))],
            [request["faces"][0]],
        ),
        "no complete matching",
    )
    refused(
        lambda: match_faces(
            solid,
            faces,
            [{**c, "area_m2": c["area_m2"] * 2} for c in request["faces"]],
        ),
        "does not cover",
    )
    return {
        "case": directory.name,
        "raw_solver_archive_rechecked": True,
        "renumbered_surfaces_preserve_support_and_load_arrays": True,
        "duplicate_geometric_match_refused": True,
        "displaced_face_refused": True,
        "wrong_area_refused": True,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qualification", type=Path)
    args = parser.parse_args()
    print(
        json.dumps(
            [check(args.qualification / f"case-{i}") for i in range(2)], indent=2
        )
    )
