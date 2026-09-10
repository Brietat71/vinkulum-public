"""CalculiX solid conventions and exact local Jacobian positivity for C3D10.

Element positivity is not global mesh injectivity or a displacement-error bound.
The curved-tetrahedron check uses cubic Bernstein bounds and exact dyadic input
coordinates, with bounded subdivision. An unresolved bound is never accepted.
"""

import itertools
import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

import numpy as np

TET_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
TET_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3))
HEX_FACES = (
    (0, 1, 2, 3),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)
NODE_COUNTS = {"C3D4": 4, "C3D8": 8, "C3D10": 10}
POINT_COUNTS = {"C3D4": 1, "C3D8": 8, "C3D10": 4}
_DL = ((-1, -1, -1), (1, 0, 0), (0, 1, 0), (0, 0, 1))
_VERTICES = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
_TRIPLES = tuple(itertools.product(range(4), repeat=3))
_GROUPS = tuple(sorted({tuple(sorted(t)) for t in _TRIPLES}))
_GROUP_INDEX = {key: index for index, key in enumerate(_GROUPS)}
_GROUP_MULTIPLICITY = tuple(
    sum(tuple(sorted(t)) == group for t in _TRIPLES) for group in _GROUPS
)


def face_nodes(kind):
    if kind == "C3D8":
        return HEX_FACES
    if kind == "C3D4":
        return TET_FACES
    if kind != "C3D10":
        raise ValueError("Unsupported solid element type.")
    return tuple(
        face
        + tuple(4 + i for i, (a, b) in enumerate(TET_EDGES) if a in face and b in face)
        for face in TET_FACES
    )


def tetra_shape(local, quadratic=True):
    """Values and reference gradients in CalculiX/VTK node order."""
    r, s, t = local
    bary = np.array((1 - r - s - t, r, s, t))
    gradients = np.array(_DL, dtype=float)
    if not quadratic:
        return bary, gradients
    values = [v * (2 * v - 1) for v in bary]
    rows = [(4 * bary[i] - 1) * gradients[i] for i in range(4)]
    for a, b in TET_EDGES:
        values.append(4 * bary[a] * bary[b])
        rows.append(4 * (bary[a] * gradients[b] + bary[b] * gradients[a]))
    return np.array(values), np.array(rows)


def tetra_quadrature(kind):
    if kind == "C3D4":
        return (((0.25, 0.25, 0.25), 1 / 6),)
    if kind != "C3D10":
        raise ValueError("Tetrahedral quadrature requires C3D4 or C3D10.")
    # Degree-two tetrahedral rule, in CalculiX gauss3d5 table order.
    small = (5 - math.sqrt(5)) / 20
    large = (5 + 3 * math.sqrt(5)) / 20
    return tuple(
        (p, 1 / 24)
        for p in (
            (small, small, small),
            (large, small, small),
            (small, large, small),
            (small, small, large),
        )
    )


def _det_columns(a, b, c):
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1])
    )


def _integer_jacobians(nodes):
    ratios = [[float(x).as_integer_ratio() for x in row] for row in nodes]
    denominator = max(d for row in ratios for _, d in row)
    points = [[n * (denominator // d) for n, d in row] for row in ratios]
    points = [[x - y for x, y in zip(row, points[0])] for row in points]
    matrices = []
    for vertex in range(4):
        bary = tuple(int(i == vertex) for i in range(4))
        gradients = [tuple((4 * bary[i] - 1) * v for v in _DL[i]) for i in range(4)]
        gradients += [
            tuple(4 * (bary[a] * _DL[b][k] + bary[b] * _DL[a][k]) for k in range(3))
            for a, b in TET_EDGES
        ]
        # A matrix is stored by columns, all entries integral.
        matrices.append(
            tuple(
                tuple(
                    sum(points[n][axis] * gradients[n][k] for n in range(10))
                    for axis in range(3)
                )
                for k in range(3)
            )
        )
    return tuple(matrices), denominator


def _bernstein_sums(matrices):
    # det(sum lambda_i J_i) is cubic. Group its 64 column products by
    # barycentric multi-index. Dividing by the positive multinomial count
    # gives the degree-three Bernstein coefficient (no interpolation solve).
    sums = [0] * 20
    for a, b, c in _TRIPLES:
        sums[_GROUP_INDEX[tuple(sorted((a, b, c)))]] += _det_columns(
            matrices[a][0], matrices[b][1], matrices[c][2]
        )
    return sums


@dataclass(frozen=True)
class JacobianCertificate:
    status: str
    cells_examined: int
    lower_bound: Fraction | None = None
    witness: tuple | None = None


def certify_tetra10(nodes, *, max_depth=10, max_cells=255):
    """Bound det(dx/dxi) over the closed reference tetrahedron, exactly.

    Input floats are treated as exact binary rationals. Subdivision averages
    integer Jacobian matrices with a shared power-of-two scale. Nonpositive
    vertex values give witnesses; exhausted budgets yield 'unresolved'.
    """
    if (
        len(nodes) != 10
        or any(len(row) != 3 for row in nodes)
        or not np.isfinite(np.array(nodes, dtype=float)).all()
    ):
        raise ValueError("Ten finite three-coordinate nodes are required.")
    if (
        type(max_depth) is not int
        or not 0 <= max_depth <= 16
        or type(max_cells) is not int
        or not 1 <= max_cells <= 4095
    ):
        raise ValueError("Invalid bounded Jacobian subdivision budget.")
    matrices, denominator = _integer_jacobians(nodes)
    pending = [(matrices, _VERTICES, 0)]
    lower, examined = None, 0
    while pending:
        matrices, vertices, depth = pending.pop()
        examined += 1
        for matrix, point in zip(matrices, vertices):
            if _det_columns(*matrix) <= 0:
                return JacobianCertificate("nonpositive", examined, witness=point)
        coefficients = _bernstein_sums(matrices)
        if min(coefficients) > 0:
            bound = min(
                Fraction(value, count * denominator**3 * 8**depth)
                for value, count in zip(coefficients, _GROUP_MULTIPLICITY)
            )
            lower = bound if lower is None else min(lower, bound)
            continue
        if depth >= max_depth or examined + len(pending) + 2 > max_cells:
            return JacobianCertificate("unresolved", examined)
        a, b = max(
            itertools.combinations(range(4), 2),
            key=lambda edge: sum(
                (x - y) ** 2 for x, y in zip(vertices[edge[0]], vertices[edge[1]])
            ),
        )
        midpoint = tuple((x + y) / 2 for x, y in zip(vertices[a], vertices[b]))
        middle_matrix = tuple(
            tuple(x + y for x, y in zip(c, d)) for c, d in zip(matrices[a], matrices[b])
        )
        doubled = tuple(tuple(tuple(2 * x for x in col) for col in m) for m in matrices)
        for endpoint in (a, b):
            child, reference = list(doubled), list(vertices)
            child[endpoint], reference[endpoint] = middle_matrix, midpoint
            pending.append((tuple(child), tuple(reference), depth + 1))
    return JacobianCertificate("positive", examined, lower_bound=lower)


@lru_cache(maxsize=8192)
def integration_weights(kind, nodes):
    """Validate element geometry and return physical quadrature weights [m³]."""
    p = np.array(nodes, dtype=float)
    if (
        kind not in NODE_COUNTS
        or p.shape != (NODE_COUNTS[kind], 3)
        or not np.isfinite(p).all()
    ):
        raise ValueError("Invalid solid element coordinates or type.")
    p = p - p[0]
    edges = p[[1, 3, 4]] if kind == "C3D8" else p[1:4]
    scale = max(np.linalg.norm(edge) for edge in edges)
    jac = edges.T / scale if scale > 0 and math.isfinite(scale) else np.zeros((3, 3))
    if not np.isfinite(jac).all() or np.linalg.det(jac) <= 1e-12:
        raise ValueError("An element is inverted, degenerate or too poorly scaled.")
    if kind == "C3D8":
        expected = np.array(
            (
                (0, 0, 0),
                (1, 0, 0),
                (1, 1, 0),
                (0, 1, 0),
                (0, 0, 1),
                (1, 0, 1),
                (1, 1, 1),
                (0, 1, 1),
            )
        )
        local = np.linalg.solve(jac, (p / scale).T).T
        if not np.isfinite(local).all() or np.max(abs(local - expected)) > 1e-10:
            raise ValueError("C3D8 geometry must be affine.")
        weights = (float(np.linalg.det(edges.T)) / 8,) * 8
    else:
        if kind == "C3D10":
            certificate = certify_tetra10(nodes)
            if certificate.status != "positive":
                raise ValueError(
                    "C3D10 Jacobian "
                    + certificate.status
                    + ": no positive bound over the whole element."
                )
            try:
                lower_bound = float(certificate.lower_bound)
            except OverflowError as error:
                raise ValueError("Element volume exceeds the numeric range.") from error
            if lower_bound / scale / scale / scale <= 1e-12:
                raise ValueError("C3D10 Jacobian bound is too poorly scaled.")
        weights = tuple(
            float(np.linalg.det(p.T @ tetra_shape(q, kind == "C3D10")[1])) * w
            for q, w in tetra_quadrature(kind)
        )
    if any(not math.isfinite(w) or w <= 0 for w in weights):
        raise ValueError("Nonpositive or nonfinite element quadrature weight.")
    return weights
