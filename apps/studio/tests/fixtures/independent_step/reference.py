"""Independent exact box integrals and parallel-axis reference; no CAD imports.

SPDX-License-Identifier: Apache-2.0
Lengths start in mm; density is kg/m^3; output quantities use SI.
"""

import json
from fractions import Fraction as F
from math import prod
from pathlib import Path


def reference():
    # Two disjoint box interiors sharing the 12 x 8 mm interface at y=12 mm.
    boxes = (((60, 12, 8), (30, 6, 4)), ((12, 36, 8), (6, 30, 4)))
    volumes = [prod(size) for size, _ in boxes]
    volume = sum(volumes)
    centre = [
        sum(F(v, volume) * c[i] for v, (_, c) in zip(volumes, boxes)) for i in range(3)
    ]
    mass = F(7800 * volume, 10**9)
    inertia = [[F(0) for _ in range(3)] for _ in range(3)]
    for v, (size, c) in zip(volumes, boxes):
        m = F(7800 * v, 10**9)
        d = [c[i] - centre[i] for i in range(3)]
        for i in range(3):
            for j in range(3):
                value = -d[i] * d[j]
                if i == j:
                    value += sum(x * x for x in d)
                    value += F(sum(size[k] ** 2 for k in range(3) if k != i), 12)
                inertia[i][j] += m * value / 10**6  # kg mm^2 -> kg m^2.

    # Independent specification of the proper rotation and placement.
    q = [[F(v, 13) for v in row] for row in ((-3, -4, 12), (12, 3, 4), (-4, 12, 3))]
    assert all(
        sum(q[i][k] * q[j][k] for k in range(3)) == (i == j)
        for i in range(3)
        for j in range(3)
    )
    placed_centre = [
        (sum(q[i][j] * centre[j] for j in range(3)) + t) / 1000
        for i, t in enumerate((17, -23, 31))
    ]
    placed_inertia = [
        [
            sum(q[i][k] * inertia[k][l] * q[j][l] for k in range(3) for l in range(3))
            for j in range(3)
        ]
        for i in range(3)
    ]
    return {
        "density_kg_m3": 7800,
        "volume_m3": float(F(volume, 10**9)),
        "mass_kg": float(mass),
        "centre_local_mm": [float(v) for v in centre],
        "centre_m": [float(v) for v in placed_centre],
        "inertia_local_kg_m2": [[float(v) for v in row] for row in inertia],
        "inertia_kg_m2": [[float(v) for v in row] for row in placed_inertia],
        "characteristic_length_m": 0.08,
        "relative_scale_tolerance": 1e-8,
    }


if __name__ == "__main__":
    Path(__file__).with_name("reference.json").write_text(
        json.dumps(reference(), indent=2) + "\n", encoding="ascii"
    )
