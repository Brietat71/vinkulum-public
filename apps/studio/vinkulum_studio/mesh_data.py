"""Validated solid mesh independent of material, loads or supports."""

from dataclasses import InitVar, dataclass
from functools import cached_property

import numpy as np

from .finite_elements import NODE_COUNTS, TET_EDGES, face_nodes, integration_weights
from .model import finite_number

MAX_NODES = 6000
MAX_ELEMENTS = 5000


@dataclass(frozen=True)
class FiniteMesh:
    nodes: tuple
    elements: tuple
    element_type: str
    cancellation: InitVar[object] = None

    def __post_init__(self, cancellation):
        if (
            not isinstance(self.element_type, str)
            or self.element_type not in NODE_COUNTS
        ):
            raise ValueError("Unsupported solid element type.")
        for name in ("nodes", "elements"):
            rows = getattr(self, name)
            if not isinstance(rows, (tuple, list)) or not all(
                isinstance(r, (tuple, list)) for r in rows
            ):
                raise ValueError("Expected rows of mesh coordinates/connectivity.")
            object.__setattr__(self, name, tuple(tuple(r) for r in rows))
        count = NODE_COUNTS[self.element_type]
        if (
            not count <= len(self.nodes) <= MAX_NODES
            or not 1 <= len(self.elements) <= MAX_ELEMENTS
        ):
            raise ValueError("Static mesh budget exceeded or mesh empty.")
        if any(
            len(row) != 3 or not all(finite_number(v) for v in row)
            for row in self.nodes
        ):
            raise ValueError("Mesh nodes require finite world coordinates in metres.")
        if any(
            len(e) != count
            or len(set(e)) != count
            or any(type(i) is not int or not 1 <= i <= len(self.nodes) for i in e)
            for e in self.elements
        ):
            raise ValueError("Invalid solid-element connectivity.")
        if len({tuple(sorted(e)) for e in self.elements}) != len(self.elements):
            raise ValueError("Duplicate elements are not allowed.")
        if {i for e in self.elements for i in e} != set(range(1, len(self.nodes) + 1)):
            raise ValueError("Every node must belong to an element.")
        object.__setattr__(self, "_cancellation", cancellation)
        try:
            _ = self.boundary_faces  # Validate before exposing the mesh.
        finally:
            object.__delattr__(self, "_cancellation")

    @cached_property
    def boundary_faces(self):
        points = np.array(self.nodes, dtype=float)
        oriented = {}
        faces, neighbours = {}, [set() for _ in self.elements]
        if self.element_type == "C3D10":
            corner_ids = {i for e in self.elements for i in e[:4]}
            mid_ids = {i for e in self.elements for i in e[4:]}
            if corner_ids & mid_ids:
                raise ValueError(
                    "A quadratic node cannot be both a corner and an edge node."
                )
            edges, midpoints = {}, {}
            for element in self.elements:
                for mid, (a, b) in zip(element[4:], TET_EDGES):
                    key = tuple(sorted((element[a], element[b])))
                    if (
                        edges.setdefault(key, mid) != mid
                        or midpoints.setdefault(mid, key) != key
                    ):
                        raise ValueError(
                            "Nonconforming quadratic edge-node identities."
                        )
        for index, element in enumerate(self.elements):
            if self._cancellation is not None and self._cancellation():
                raise InterruptedError("Mesh validation cancelled.")
            p = points[np.array(element) - 1]
            integration_weights(
                self.element_type, tuple(tuple(float(v) for v in row) for row in p)
            )
            for face in face_nodes(self.element_type):
                key = tuple(sorted(element[i] for i in face))
                oriented.setdefault(key, tuple(element[i] for i in face))
                owners = faces.setdefault(key, [])
                owners.append(index)
                if len(owners) > 2:
                    raise ValueError("A mesh face has more than two owners.")
                if len(owners) == 2:
                    left, right = owners
                    neighbours[left].add(right)
                    neighbours[right].add(left)
        seen, pending = set(), [0]
        while pending:
            index = pending.pop()
            if index not in seen:
                seen.add(index)
                pending.extend(neighbours[index] - seen)
        if len(seen) != len(self.elements):
            raise ValueError(
                "Elements must form one connected mesh through shared faces."
            )
        return tuple(oriented[key] for key in sorted(faces) if len(faces[key]) == 1)
