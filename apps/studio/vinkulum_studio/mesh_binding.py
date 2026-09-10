"""Explicit physical intent on captured, outward-oriented mesh surfaces."""

import uuid
from dataclasses import asdict, dataclass
from functools import cached_property, lru_cache

import numpy as np

from .meshing import MeshedSolid, _fields, fingerprint
from .model import finite_number


def triangle_shape(r, s, quadratic):
    bary = np.array((1 - r - s, r, s))
    grad = np.array(((-1.0, -1.0), (1.0, 0.0), (0.0, 1.0)))
    if not quadratic:
        return bary, grad
    values = [v * (2 * v - 1) for v in bary]
    rows = [(4 * bary[i] - 1) * grad[i] for i in range(3)]
    for a, b in ((0, 1), (1, 2), (2, 0)):
        values.append(4 * bary[a] * bary[b])
        rows.append(4 * (bary[a] * grad[b] + bary[b] * grad[a]))
    return np.array(values), np.array(rows)


def triangle_quadrature():
    # Seven-point degree-five symmetric rule. Weights integrate a triangle
    # of area 1/2, not area 1. Every weight is positive.
    yield (1 / 3, 1 / 3), 0.1125
    for a, b, w in (
        (0.059715871789770, 0.470142064105115, 0.066197076394253),
        (0.797426985353087, 0.101286507323456, 0.062969590272414),
    ):
        for bary in ((a, b, b), (b, a, b), (b, b, a)):
            yield bary[1:], w


@lru_cache(maxsize=20000)
def surface_integrals(points):
    """Nodal area weights and integrated outward area vectors, in m².

    N*normal is polynomial and the rule is exact in exact arithmetic. N*|normal|
    on a curved face uses quadrature; its area is not a CAD-area certificate.
    """
    p = np.array(points, dtype=float)
    n = len(points)
    if p.shape not in ((3, 3), (6, 3)) or not np.isfinite(p).all():
        raise ValueError("A finite linear/quadratic boundary triangle is required.")
    local = p - p[0]
    area = np.zeros(n)
    normal = np.zeros((n, 3))
    for (r, s), weight in triangle_quadrature():
        basis, gradient = triangle_shape(r, s, n == 6)
        tangent = local.T @ gradient
        vector = np.cross(tangent[:, 0], tangent[:, 1])
        magnitude = np.linalg.norm(vector)
        if not np.isfinite(magnitude) or magnitude <= 0:
            raise ValueError("Degenerate mesh boundary triangle.")
        area += weight * magnitude * basis
        normal += weight * basis[:, None] * vector
    if not np.isfinite(area).all() or not np.isfinite(normal).all() or area.sum() <= 0:
        raise ValueError("Nonfinite mesh surface integral.")
    return tuple(map(float, area)), tuple(tuple(map(float, row)) for row in normal)


@dataclass(frozen=True)
class MeshCondition:
    id: str
    name: str
    kind: str
    surfaces: tuple
    axes: tuple = ()
    values: tuple = ()

    def __post_init__(self):
        if not isinstance(self.id, str) or str(uuid.UUID(self.id)) != self.id:
            raise ValueError("A canonical boundary-condition UUID is required.")
        if not isinstance(self.name, str) or not 1 <= len(self.name.strip()) <= 128:
            raise ValueError("A boundary-condition name is required.")
        if self.kind not in ("support", "pressure", "total_force"):
            raise ValueError("Unknown boundary-condition type.")
        for key in ("surfaces", "axes", "values"):
            if not isinstance(getattr(self, key), tuple):
                raise TypeError("Immutable boundary-condition parameters are required.")
        if (
            not self.surfaces
            or len(set(self.surfaces)) != len(self.surfaces)
            or any(type(v) is not int or v < 1 for v in self.surfaces)
        ):
            raise ValueError("Select distinct captured mesh faces.")
        if self.kind == "support":
            if (
                not self.axes
                or len(set(self.axes)) != len(self.axes)
                or any(type(a) is not int or a not in (1, 2, 3) for a in self.axes)
                or self.values
            ):
                raise ValueError("A zero support requires world axes X, Y or Z.")
        elif (
            self.axes
            or len(self.values) != (1 if self.kind == "pressure" else 3)
            or not all(finite_number(v) for v in self.values)
        ):
            raise ValueError(
                "A pressure requires Pa; a total force requires three world components in N."
            )

    @classmethod
    def from_dict(cls, data):
        data = _fields(data, cls)
        for key in ("surfaces", "axes", "values"):
            data[key] = tuple(data[key])
        return cls(**data)


@dataclass(frozen=True)
class MeshBinding:
    solid: MeshedSolid
    conditions: tuple
    load_factor: float = 1.0

    def __post_init__(self):
        if (
            not isinstance(self.solid, MeshedSolid)
            or not isinstance(self.conditions, tuple)
            or len(self.conditions) > 100
            or not all(isinstance(c, MeshCondition) for c in self.conditions)
        ):
            raise ValueError(
                "A captured mesh and at most 100 boundary conditions are required."
            )
        if len({c.id for c in self.conditions}) != len(self.conditions):
            raise ValueError("Duplicate boundary-condition identity.")
        available = {s.id for s in self.solid.surfaces}
        if any(not set(c.surfaces) <= available for c in self.conditions):
            raise ValueError("A boundary condition references a missing captured face.")
        if not finite_number(self.load_factor):
            raise ValueError("Boundary load multiplier must be finite.")

    @cached_property
    def arrays(self):
        mesh = self.solid.mesh
        points = np.array(mesh.nodes)
        forces = np.zeros((len(points), 3))
        fixed = set()
        groups = {s.id: s for s in self.solid.surfaces}
        for condition in self.conditions:
            faces = [row for i in condition.surfaces for row in groups[i].triangles]
            if condition.kind == "support":
                fixed.update(
                    (i, a) for face in faces for i in face for a in condition.axes
                )
                continue
            local = np.zeros_like(forces)
            area = 0.0
            for face in faces:
                ids = np.array(face) - 1
                weights, vectors = surface_integrals(
                    tuple(tuple(p) for p in points[ids])
                )
                if condition.kind == "pressure":
                    local[ids] -= condition.values[0] * np.array(vectors)
                else:
                    local[ids] += np.array(weights)[:, None] * np.array(
                        condition.values
                    )
                    area += sum(weights)
            if condition.kind == "total_force":
                if not np.isfinite(area) or area <= 0:
                    raise ValueError("Selected faces have no finite positive area.")
                local /= area
            forces += self.load_factor * local
        if not np.isfinite(forces).all():
            raise ValueError("Boundary loads exceed the numeric range.")
        rows = tuple(
            (i + 1, *map(float, row))
            for i, row in enumerate(forces)
            if np.any(row != 0)
        )
        return tuple(sorted(fixed)), rows

    @property
    def sha256(self):
        return fingerprint(asdict(self))

    def verify_study(self, study):
        mesh = self.solid.mesh
        fixed, forces = self.arrays
        if (
            study.nodes != mesh.nodes
            or study.elements != mesh.elements
            or study.element_type != mesh.element_type
            or study.fixed_dofs != fixed
            or study.forces != forces
        ):
            raise ValueError(
                "Study mesh/supports/loads do not match their captured CAD mesh and physical intent."
            )

    def study(self, young_pa, poisson):
        from .calculix import StaticStudy

        fixed, forces = self.arrays
        mesh = self.solid.mesh
        return StaticStudy(
            mesh.nodes,
            mesh.elements,
            fixed,
            forces,
            young_pa,
            poisson,
            mesh.element_type,
            mesh_binding=self,
        )

    @classmethod
    def from_dict(cls, data):
        data = _fields(data, cls)
        data["solid"] = MeshedSolid.from_dict(data["solid"])
        data["conditions"] = tuple(
            MeshCondition.from_dict(c) for c in data["conditions"]
        )
        return cls(**data)


def mesh_measurements(solid):
    points = np.array(solid.mesh.nodes) - np.array(solid.request.body.position)
    area = 0.0
    volume = 0.0
    for surface in solid.surfaces:
        for row in surface.triangles:
            p = points[np.array(row) - 1]
            weights, normals = surface_integrals(tuple(map(tuple, p)))
            area += sum(weights)
            volume += np.sum(p * np.array(normals)) / 3
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError("Mesh has no finite positive enclosed volume.")
    reference = solid.request.body.cad.volume_m3
    return {
        "surface_area_quadrature_m2": float(area),
        "signed_mesh_volume_m3": float(volume),
        "cad_volume_m3": reference,
        "relative_volume_difference": float(abs(volume / reference - 1)),
    }
