"""Independent signed-volume integrals for the analytic STEP corpus.

SPDX-License-Identifier: Apache-2.0
No generator, Vinkulum code or CAD library is imported. All polynomial
coefficients and transformations are rational; pi enters only the final output.
"""

import json
import math
from fractions import Fraction as F
from pathlib import Path


def box(dimensions, centre, sign):
    x, y, z = map(F, dimensions)
    volume = sign * x * y * z
    diagonal = ((y * y + z * z) / 12, (x * x + z * z) / 12, (x * x + y * y) / 12)
    return volume, tuple(map(F, centre)), diagonal


def cylinder(radius, height, centre, sign):
    r, h = F(radius), F(height)
    # A common factor pi cancels in the weighted centroid.
    volume = sign * r * r * h
    diagonal = ((3 * r * r + h * h) / 12, (3 * r * r + h * h) / 12, r * r / 2)
    return volume, tuple(map(F, centre)), diagonal


def rotation():
    # Unit quaternion (1,2,3,4)/sqrt(30), evaluated without a square root.
    w, x, y, z = map(F, (1, 2, 3, 4))
    norm = w * w + x * x + y * y + z * z
    q = (
        (w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y)),
        (2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x)),
        (2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z),
    )
    q = tuple(tuple(v / norm for v in row) for row in q)
    assert all(
        sum(q[i][k] * q[j][k] for k in range(3)) == (i == j)
        for i in range(3)
        for j in range(3)
    )
    determinant = sum(
        q[0][i]
        * (
            q[1][(i + 1) % 3] * q[2][(i + 2) % 3]
            - q[1][(i + 2) % 3] * q[2][(i + 1) % 3]
        )
        for i in range(3)
    )
    assert determinant == 1
    return q


def integrate(parts, factor, density, translation, length, unit):
    volume = sum(part[0] for part in parts)
    centre = tuple(sum(v * c[i] for v, c, d in parts) / volume for i in range(3))
    inertia = [[F(0) for _ in range(3)] for _ in range(3)]
    for v, c, diagonal in parts:
        delta = tuple(value - origin for value, origin in zip(c, centre))
        squared = sum(value * value for value in delta)
        for i in range(3):
            for j in range(3):
                inertia[i][j] += v * (
                    (diagonal[i] + squared if i == j else 0) - delta[i] * delta[j]
                )
    q = rotation()
    placed = tuple(
        sum(q[i][j] * centre[j] for j in range(3)) + translation[i] for i in range(3)
    )
    rotated = [
        [
            sum(q[i][k] * inertia[k][l] * q[j][l] for k in range(3) for l in range(3))
            for j in range(3)
        ]
        for i in range(3)
    ]
    scale = math.pi if factor == "pi" else 1
    return {
        "density_kg_m3": density,
        "step_length_unit": unit,
        "volume_mm3_coefficient": str(volume),
        "analytic_factor": factor,
        "local_centroid_mm_rational": [str(v) for v in centre],
        "local_inertia_mm5_coefficient": [[str(v) for v in row] for row in inertia],
        "rotation": [[float(v) for v in row] for row in q],
        "translation_mm": list(translation),
        "volume_m3": float(volume) * scale * 1e-9,
        "mass_kg": float(volume) * scale * 1e-9 * density,
        "centre_m": [float(v) * 1e-3 for v in placed],
        "inertia_kg_m2": [
            [float(v) * scale * density * 1e-15 for v in row] for row in rotated
        ],
        "characteristic_length_m": length,
        "relative_scale_tolerance": 1e-8,
    }


def references():
    return {
        "eccentric-bore-metres": integrate(
            [cylinder(25, 18, (0, 0, 9), 1), cylinder(6, 18, (8, -5, 9), -1)],
            "pi",
            8100,
            (300, -140, 220),
            0.08,
            "metre",
        ),
        "offset-closed-void-mm": integrate(
            [box((60, 44, 30), (30, 22, 15), 1), box((22, 18, 12), (36, 17, 16), -1)],
            "1",
            2700,
            (-80, 25, 41),
            0.09,
            "millimetre",
        ),
    }


if __name__ == "__main__":
    Path(__file__).with_name("reference.json").write_text(
        json.dumps(references(), indent=2, allow_nan=False) + "\n"
    )
