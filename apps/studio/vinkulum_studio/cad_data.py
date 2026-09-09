"""Immutable CAD payload, usable by the mechanics worker without loading OCCT."""

import math
from dataclasses import dataclass, fields

import numpy as np


@dataclass(frozen=True)
class CadGeometry:
    brep_mm: str
    volume_m3: float
    unit_inertia_m2: tuple
    vertices_m: tuple
    triangles: tuple
    operations: tuple
    occt_version: str
    build123d_version: str

    def __post_init__(self):
        if not isinstance(self.brep_mm, str) or not 1 <= len(self.brep_mm) <= 2_000_000:
            raise ValueError("Missing or oversized BREP (2 MB maximum).")
        if not self.brep_mm.isascii() or "CASCADE Topology" not in self.brep_mm[:200]:
            raise ValueError("An ASCII OCCT BREP is required.")
        if (
            not isinstance(self.volume_m3, (int, float))
            or not math.isfinite(self.volume_m3)
            or self.volume_m3 <= 0
        ):
            raise ValueError("CAD volume must be positive.")
        if (
            not isinstance(self.unit_inertia_m2, tuple)
            or len(self.unit_inertia_m2) != 9
            or not all(
                type(v) in (int, float) and math.isfinite(v)
                for v in self.unit_inertia_m2
            )
        ):
            raise ValueError("A 3×3 CAD inertia tensor is required.")
        inertia = np.array(self.unit_inertia_m2, dtype=float).reshape(3, 3)
        if not np.isfinite(inertia).all() or not np.allclose(
            inertia, inertia.T, rtol=1e-10, atol=1e-16
        ):
            raise ValueError("The CAD tensor must be finite and symmetric.")
        eigen = np.linalg.eigvalsh(inertia)
        if eigen[0] <= 0 or eigen[2] > sum(eigen[:2]) * (1 + 1e-10):
            raise ValueError("Nonphysical CAD inertia tensor.")
        if (
            not isinstance(self.vertices_m, tuple)
            or not 4 <= len(self.vertices_m) <= 30_000
        ):
            raise ValueError("CAD vertex limit exceeded.")
        for point in self.vertices_m:
            if (
                not isinstance(point, tuple)
                or len(point) != 3
                or not all(
                    isinstance(v, (int, float)) and math.isfinite(v) and abs(v) <= 1000
                    for v in point
                )
            ):
                raise ValueError("Invalid CAD vertex.")
        if (
            not isinstance(self.triangles, tuple)
            or not 4 <= len(self.triangles) <= 30_000
        ):
            raise ValueError("CAD triangle limit exceeded.")
        for triangle in self.triangles:
            if (
                not isinstance(triangle, tuple)
                or len(triangle) != 3
                or len(set(triangle)) != 3
                or not all(
                    type(i) is int and 0 <= i < len(self.vertices_m) for i in triangle
                )
            ):
                raise ValueError("Invalid CAD triangle.")
        if (
            not isinstance(self.operations, tuple)
            or len(self.operations) > 100
            or not all(isinstance(v, str) and len(v) <= 2048 for v in self.operations)
        ):
            raise ValueError("Invalid CAD history.")
        for value in (self.occt_version, self.build123d_version):
            if not isinstance(value, str) or not 1 <= len(value) <= 64:
                raise ValueError("Missing CAD version.")
        if int(self.occt_version.split(".")[0]) < 8:
            raise ValueError("OCCT 8 or later is required.")

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict) or set(data) != {f.name for f in fields(cls)}:
            raise ValueError("Unknown or missing CAD fields.")
        data = dict(data)
        for key in ("unit_inertia_m2", "operations"):
            data[key] = tuple(data[key])
        for key in ("vertices_m", "triangles"):
            data[key] = tuple(tuple(row) for row in data[key])
        return cls(**data)
