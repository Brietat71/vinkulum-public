"""Planar line sketches and exact affine dimension constraints, independent of CAD/Qt.

Coordinates are display/geometry seeds in mm. Dimension literals are bounded
decimal strings: 0.1 + 0.2 = 0.3 is an exact relation in the constraint graph.
Each equation is a difference of two coordinates, or an absolute coordinate.
Connected components give exact consistency, rank and remaining translations.
This is not a nonlinear angle, radius or Euclidean-distance solver.
"""

import math
import re
from collections import deque
from dataclasses import asdict, dataclass, fields, replace
from decimal import Decimal, localcontext
from fractions import Fraction
from uuid import UUID, uuid4

KINDS = {
    "horizontal": "Horizontal",
    "vertical": "Vertical",
    "fixed": "Fixed point",
    "distance_x": "X distance",
    "distance_y": "Y distance",
}
DECIMAL = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d{1,2})?\Z")


def identifier(value):
    if not isinstance(value, str) or str(UUID(value)) != value:
        raise ValueError("A canonical sketch entity UUID is required.")


def label(value):
    if not isinstance(value, str) or not 1 <= len(value.strip()) <= 128:
        raise ValueError("Sketch names require 1–128 characters.")


def decimal_mm(value):
    if (
        not isinstance(value, str)
        or not 1 <= len(value) <= 64
        or not DECIMAL.fullmatch(value)
    ):
        raise ValueError(
            "Enter a bounded decimal dimension in mm, such as 12.5 or -2e-3."
        )
    result = Fraction(value)
    if abs(result) > 1_000_000:
        raise ValueError("Sketch dimensions must be within ±10⁶ mm.")
    return result


def members(data, cls):
    if not isinstance(data, dict) or set(data) != {f.name for f in fields(cls)}:
        raise ValueError(f"Unknown or missing {cls.__name__} fields.")
    return dict(data)


@dataclass(frozen=True)
class SketchPoint:
    id: str
    name: str
    xy_mm: tuple

    def __post_init__(self):
        identifier(self.id)
        label(self.name)
        if (
            not isinstance(self.xy_mm, tuple)
            or len(self.xy_mm) != 2
            or not all(
                type(v) in (float, int) and math.isfinite(v) and abs(v) <= 1e6
                for v in self.xy_mm
            )
        ):
            raise ValueError(
                "Sketch points require finite X/Y coordinates within ±10⁶ mm."
            )


@dataclass(frozen=True)
class SketchEdge:
    id: str
    name: str
    start: str
    end: str

    def __post_init__(self):
        identifier(self.id)
        label(self.name)
        identifier(self.start)
        identifier(self.end)
        if self.start == self.end:
            raise ValueError("A sketch segment requires two distinct point identities.")


@dataclass(frozen=True)
class SketchConstraint:
    id: str
    name: str
    kind: str
    entities: tuple
    values_mm: tuple = ()

    def __post_init__(self):
        identifier(self.id)
        label(self.name)
        if self.kind not in KINDS:
            raise ValueError("Unsupported sketch constraint kind.")
        count = 2 if self.kind.startswith("distance_") else 1
        if (
            not isinstance(self.entities, tuple)
            or len(self.entities) != count
            or len(set(self.entities)) != count
        ):
            raise ValueError("Invalid sketch constraint references.")
        for value in self.entities:
            identifier(value)
        count = (
            2 if self.kind == "fixed" else 1 if self.kind.startswith("distance_") else 0
        )
        if not isinstance(self.values_mm, tuple) or len(self.values_mm) != count:
            raise ValueError("Invalid number of sketch dimensions.")
        for value in self.values_mm:
            decimal_mm(value)


@dataclass(frozen=True)
class Sketch:
    id: str
    name: str
    points: tuple = ()
    edges: tuple = ()
    constraints: tuple = ()

    def __post_init__(self):
        identifier(self.id)
        label(self.name)
        for rows, cls, limit in (
            (self.points, SketchPoint, 64),
            (self.edges, SketchEdge, 64),
            (self.constraints, SketchConstraint, 128),
        ):
            if (
                not isinstance(rows, tuple)
                or len(rows) > limit
                or not all(isinstance(r, cls) for r in rows)
            ):
                raise ValueError(
                    f"A sketch admits at most {limit} {cls.__name__} objects."
                )
        keys = [self.id] + [
            r.id for rows in (self.points, self.edges, self.constraints) for r in rows
        ]
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate sketch entity identities.")
        points, edges = {p.id for p in self.points}, {e.id for e in self.edges}
        if any(e.start not in points or e.end not in points for e in self.edges):
            raise ValueError("A sketch segment references a missing point.")
        if len({frozenset((e.start, e.end)) for e in self.edges}) != len(self.edges):
            raise ValueError("Duplicate sketch segments between the same points.")
        for c in self.constraints:
            available = edges if c.kind in ("horizontal", "vertical") else points
            if not set(c.entities) <= available:
                raise ValueError(f"Constraint '{c.name}' references a missing entity.")

    @classmethod
    def from_dict(cls, data):
        data = members(data, cls)
        for key, item_type, tuple_fields in (
            ("points", SketchPoint, ("xy_mm",)),
            ("edges", SketchEdge, ()),
            ("constraints", SketchConstraint, ("entities", "values_mm")),
        ):
            result = []
            for row in data[key]:
                row = members(row, item_type)
                for field in tuple_fields:
                    row[field] = tuple(row[field])
                result.append(item_type(**row))
            data[key] = tuple(result)
        return cls(**data)

    def canonical(self):
        data = asdict(self)
        for key in ("points", "edges", "constraints"):
            data[key] = sorted(data[key], key=lambda row: row["id"])
        return data

    def polygon(self):
        """One simple closed line loop, oriented +Z; test rounded OCCT inputs exactly."""
        result = solve(self)
        if result.conflicts:
            raise ValueError("Resolve conflicting sketch dimensions before extrusion.")
        if len(self.points) < 3 or len(self.edges) != len(self.points):
            raise ValueError(
                "Extrusion requires one closed profile with at least three segments."
            )
        graph = {p.id: [] for p in self.points}
        for edge in self.edges:
            graph[edge.start].append(edge.end)
            graph[edge.end].append(edge.start)
        if any(len(neighbors) != 2 for neighbors in graph.values()):
            raise ValueError(
                "Close the profile; branches and unused points cannot be extruded."
            )
        first, previous, order = min(graph), None, []
        current = first
        while current not in order:
            order.append(current)
            following = next(p for p in sorted(graph[current]) if p != previous)
            previous, current = current, following
        if current != first or len(order) != len(graph):
            raise ValueError(
                "Extrusion requires one connected profile; nested or separate loops are not admitted."
            )
        points = [result.coordinates[p] for p in order]
        exact = [tuple(Fraction.from_float(float(v)) for v in p) for p in points]
        pairs = list(zip(exact, exact[1:] + exact[:1]))
        if any(
            sum((v - u) ** 2 for u, v in zip(a, b)) < Fraction(1, 1_000_000)
            for a, b in pairs
        ):
            raise ValueError("Every sketch segment must be at least 0.001 mm long.")
        for i, (a, b) in enumerate(pairs):
            for j in range(i + 1, len(pairs)):
                c, d = pairs[j]
                if j == i + 1 or (i == 0 and j == len(pairs) - 1):
                    # Adjacent segments may meet only at their common endpoint.
                    outer, shared, end = (a, b, d) if j == i + 1 else (b, a, c)
                    if (
                        orientation(outer, shared, end) == 0
                        and sum(
                            (u - v) * (w - v) for u, v, w in zip(outer, shared, end)
                        )
                        > 0
                    ):
                        raise ValueError("Adjacent sketch segments overlap.")
                elif intersects(a, b, c, d):
                    raise ValueError("The sketch profile crosses or touches itself.")
        twice_area = sum(a[0] * b[1] - a[1] * b[0] for a, b in pairs)
        if not twice_area:
            raise ValueError("The sketch profile has zero area.")
        return tuple(points if twice_area > 0 else reversed(points))


def orientation(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def intersects(a, b, c, d):
    def on(p, q, r):
        return orientation(p, q, r) == 0 and all(
            min(u, v) <= w <= max(u, v) for u, v, w in zip(p, q, r)
        )

    if on(a, b, c) or on(a, b, d) or on(c, d, a) or on(c, d, b):
        return True
    return (orientation(a, b, c) > 0) != (orientation(a, b, d) > 0) and (
        orientation(c, d, a) > 0
    ) != (orientation(c, d, b) > 0)


@dataclass(frozen=True)
class SketchSolution:
    coordinates: dict
    free_axes: dict
    degrees_of_freedom: int | None
    redundant_equations: int
    conflicts: tuple = ()
    conflict_residual_mm: Fraction = Fraction(0)
    max_rounding_error_mm: float = 0.0
    exact_coordinates: dict | None = None


def solve(sketch, *, drag=None):
    """Exact affine solve; free translations minimize squared distance to seeds.

    Drag chooses a free component's translation through the selected coordinate.
    Anchored components stay fixed. No constraint is relaxed to fit a conflict.
    """
    edges = {e.id: e for e in sketch.edges}
    points = {p.id: p for p in sketch.points}
    rows = ([], [])
    anchor = "origin"
    for c in sorted(sketch.constraints, key=lambda c: c.id):
        if c.kind in ("horizontal", "vertical"):
            e = edges[c.entities[0]]
            rows[int(c.kind == "horizontal")].append(
                (e.start, e.end, Fraction(0), c.id)
            )
        elif c.kind == "fixed":
            for axis in (0, 1):
                rows[axis].append(
                    (anchor, c.entities[0], decimal_mm(c.values_mm[axis]), c.id)
                )
        else:
            rows[int(c.kind == "distance_y")].append(
                (*c.entities, decimal_mm(c.values_mm[0]), c.id)
            )
    coordinates = {p: [None, None] for p in points}
    exact_coordinates = {p: [None, None] for p in points}
    free = {p: [False, False] for p in points}
    degrees, rounding = 0, Fraction(0)
    for axis, equations in enumerate(rows):
        graph = {p: [] for p in (*points, anchor)}
        for a, b, difference, source in equations:
            graph[a].append((b, difference, source))
            graph[b].append((a, -difference, source))
        visited = set()
        for root in [anchor, *sorted(points)]:
            if root in visited:
                continue
            potential, parents = {root: Fraction(0)}, {root: None}
            queue = deque((root,))
            while queue:
                a = queue.popleft()
                for b, difference, source in graph[a]:
                    value = potential[a] + difference
                    if b in potential:
                        if value != potential[b]:

                            def path(node, parents=parents):
                                chain = []
                                while parents[node] is not None:
                                    node, constraint = parents[node]
                                    chain.append(constraint)
                                return chain

                            pa, pb = path(a), path(b)
                            while pa and pb and pa[-1] == pb[-1]:
                                pa.pop()
                                pb.pop()
                            return SketchSolution(
                                {},
                                {},
                                None,
                                0,
                                tuple(sorted(set(pa + pb + [source]))),
                                value - potential[b],
                            )
                    else:
                        potential[b] = value
                        parents[b] = (a, source)
                        queue.append(b)
            visited.update(potential)
            nodes = [p for p in potential if p != anchor]
            if not nodes:
                continue
            movable = root != anchor
            if movable:
                degrees += 1
                if drag is not None and drag[0] in potential:
                    shift = (
                        Fraction.from_float(float(drag[1][axis])) - potential[drag[0]]
                    )
                else:
                    shift = sum(
                        Fraction.from_float(float(points[p].xy_mm[axis])) - potential[p]
                        for p in nodes
                    ) / len(nodes)
            else:
                shift = Fraction(0)
            for p in nodes:
                exact = potential[p] + shift
                if abs(exact) > 1_000_000:
                    raise ValueError("Solved sketch coordinates exceed ±10⁶ mm.")
                value = float(exact)
                rounding = max(rounding, abs(Fraction.from_float(value) - exact))
                coordinates[p][axis], free[p][axis] = value, movable
                exact_coordinates[p][axis] = exact
    # A difference equation suffers at most twice this coordinate transport error.
    if rounding > Fraction(1, 1_000_000_000):
        raise ValueError(
            "Sketch coordinate transport exceeds the 10⁻⁹ mm display/CAD budget."
        )
    return SketchSolution(
        {p: tuple(v) for p, v in coordinates.items()},
        {p: tuple(v) for p, v in free.items()},
        degrees,
        sum(map(len, rows)) - (2 * len(points) - degrees),
        max_rounding_error_mm=float(rounding),
        exact_coordinates={p: tuple(v) for p, v in exact_coordinates.items()},
    )


def dimension_literal(value):
    """Keep finite decimal relations exact; bound a newly chosen free dimension."""
    denominator = value.denominator
    powers = []
    for prime in (2, 5):
        power = 0
        while denominator % prime == 0:
            denominator //= prime
            power += 1
        powers.append(power)
    with localcontext() as context:
        context.prec = (
            max(40, len(str(abs(value.numerator))) + max(powers) + 2)
            if denominator == 1
            else 40
        )
        text = str((Decimal(value.numerator) / Decimal(value.denominator)).normalize())
    decimal_mm(text)
    return text


def settled(sketch, *, drag=None):
    solution = solve(sketch, drag=drag)
    if solution.conflicts:
        return sketch
    return replace(
        sketch,
        points=tuple(
            replace(p, xy_mm=solution.coordinates[p.id]) for p in sketch.points
        ),
    )


def profile(points, name="Profile"):
    vertices = tuple(
        SketchPoint(str(uuid4()), f"Point {i + 1}", tuple(p))
        for i, p in enumerate(points)
    )
    edges = tuple(
        SketchEdge(str(uuid4()), f"Segment {i + 1}", a.id, b.id)
        for i, (a, b) in enumerate(zip(vertices, vertices[1:] + vertices[:1]))
    )
    return Sketch(str(uuid4()), name, vertices, edges)


def bracket_profile():
    sketch = profile(
        ((0, 0), (80, 0), (80, 20), (30, 20), (30, 60), (0, 60)), "Dimensioned bracket"
    )
    constraints = []
    for i, e in enumerate(sketch.edges):
        kind = "horizontal" if i % 2 == 0 else "vertical"
        constraints.append(
            SketchConstraint(str(uuid4()), f"{KINDS[kind]} {i + 1}", kind, (e.id,))
        )
    p = sketch.points
    constraints.append(
        SketchConstraint(str(uuid4()), "Origin", "fixed", (p[0].id,), ("0", "0"))
    )
    for name, kind, other, value in (
        ("Width", "distance_x", 1, "80"),
        ("Height", "distance_y", 5, "60"),
        ("Web width", "distance_x", 3, "30"),
        ("Foot height", "distance_y", 2, "20"),
    ):
        constraints.append(
            SketchConstraint(str(uuid4()), name, kind, (p[0].id, p[other].id), (value,))
        )
    return replace(sketch, constraints=tuple(constraints))
