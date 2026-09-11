"""Exercise actual OCCT witness coverage, including rigid poses and false partitions.

SPDX-License-Identifier: Apache-2.0
Run with the OCCT 8 engine Python; no FreeCAD GUI or solver is required.
"""

import argparse
import hashlib
import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import build123d as bd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "apps/freecad"))
from cad_face_witness import coverage


def verify():
    checks = []
    whole = bd.Face.make_rect(20, 10)
    left = bd.Face.make_rect(10, 10).translate((-5, 0, 0))
    right = bd.Face.make_rect(10, 10).translate((5, 0, 0))
    # Independent construction of the world face: build123d axis rotations,
    # compared with the engine's row-major matrix and metre translation.
    a, b = np.deg2rad([23, 37])
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    rz = np.array([[np.cos(b), -np.sin(b), 0], [np.sin(b), np.cos(b), 0], [0, 0, 1]])
    for rotated in (False, True):
        matrix = rz @ rx if rotated else np.eye(3)
        position = (0.123, -0.041, 0.078) if rotated else (0, 0, 0)
        reference = (
            whole.rotate(bd.Axis.X, 23).rotate(bd.Axis.Z, 37).translate((123, -41, 78))
            if rotated else whole
        )
        cases = [
            ("whole", [whole], True),
            ("partition", [left, right], True),
            ("reversed partition", [right, left], True),
            ("missing half", [left], False),
            ("duplicated half", [left, left], False),
            ("overlap", [whole, left], False),
            ("wrong plane", [whole.translate((0, 0, 1))], False),
            ("outside", [whole.translate((30, 0, 0))], False),
        ]
        for name, faces, expected in cases:
            with tempfile.TemporaryDirectory(prefix="vinkulum-coverage-") as tmp:
                root = Path(tmp)
                hashes = []
                for identifier, face in enumerate(faces, 1):
                    path = root / f"cad-face-{identifier}.brep"
                    bd.export_brep(face, path)
                    hashes.append((identifier, hashlib.sha256(path.read_bytes()).hexdigest()))
                solid = SimpleNamespace(
                    request=SimpleNamespace(body=SimpleNamespace(
                        orientation=tuple(matrix.flat), position=position,
                    )), cad_faces=hashes,
                )
                try:
                    area, centre = coverage(solid, root, range(1, len(faces) + 1), reference)
                except ValueError:
                    if expected:
                        raise
                else:
                    if not expected:
                        raise AssertionError(f"Invalid coverage admitted: {name}")
                    np.testing.assert_allclose(area, 0.0002, rtol=1e-12, atol=0)
                    np.testing.assert_allclose(centre, position, rtol=0, atol=1e-14)
                checks.append(f"{'rotated/translated' if rotated else 'identity'}: {name}")
    return {"status": "passed", "checks": checks}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = verify()
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
