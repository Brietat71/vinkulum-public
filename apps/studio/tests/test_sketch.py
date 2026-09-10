"""Independent affine rank/consistency references and exact profile admission."""

import json
import unittest
from dataclasses import asdict, replace
from fractions import Fraction
from uuid import uuid4

import numpy as np
from vinkulum_studio.sketch import (
    Sketch,
    SketchConstraint,
    SketchEdge,
    SketchPoint,
    bracket_profile,
    decimal_mm,
    dimension_literal,
    profile,
    settled,
    solve,
)


def constraint(kind, entities, values=(), name=None):
    return SketchConstraint(
        str(uuid4()), name or kind, kind, tuple(entities), tuple(values)
    )


class AffineSketch(unittest.TestCase):
    def test_decimal_cycle_is_exact_and_conflict_identifies_its_equations(self):
        literal = "0." + "123456789" * 6
        self.assertEqual(
            decimal_mm(dimension_literal(decimal_mm(literal))), decimal_mm(literal)
        )
        p = tuple(SketchPoint(str(uuid4()), str(i), (i, i)) for i in range(3))
        c = (
            constraint("distance_x", (p[0].id, p[1].id), ("0.1",)),
            constraint("distance_x", (p[1].id, p[2].id), ("0.2",)),
            constraint("distance_x", (p[0].id, p[2].id), ("0.3",)),
        )
        sketch = Sketch(str(uuid4()), "Decimal triangle", p, (), c)
        result = solve(sketch)
        self.assertFalse(result.conflicts)
        self.assertEqual(result.degrees_of_freedom, 4)  # common X, three independent Y
        self.assertEqual(result.redundant_equations, 1)
        self.assertLessEqual(result.max_rounding_error_mm, 1e-9)
        changed = replace(c[2], values_mm=("0.30000000000000000001",))
        result = solve(replace(sketch, constraints=(c[0], c[1], changed)))
        self.assertEqual(set(result.conflicts), {r.id for r in c})
        self.assertEqual(abs(result.conflict_residual_mm), Fraction(1, 10**20))
        self.assertIsNone(result.degrees_of_freedom)
        self.assertFalse(result.coordinates)

    def test_remaining_motion_matches_independent_matrix_rank(self):
        rng = np.random.default_rng(642)
        for count in range(2, 13):
            with self.subTest(points=count):
                p = tuple(
                    SketchPoint(
                        str(uuid4()), str(i), tuple(map(float, rng.integers(-8, 8, 2)))
                    )
                    for i in range(count)
                )
                target = rng.integers(-20, 20, (count, 2))
                constraints, matrix = [], []
                for axis in (0, 1):
                    for i in range(count - 1):
                        if rng.random() < 0.8:
                            j = i + 1
                            constraints.append(
                                constraint(
                                    "distance_x" if axis == 0 else "distance_y",
                                    (p[i].id, p[j].id),
                                    (str(target[j, axis] - target[i, axis]),),
                                )
                            )
                            row = np.zeros(2 * count)
                            row[2 * i + axis], row[2 * j + axis] = -1, 1
                            matrix.append(row)
                constraints.append(
                    constraint("fixed", (p[0].id,), tuple(map(str, target[0])))
                )
                matrix += [np.eye(2 * count)[0], np.eye(2 * count)[1]]
                sketch = Sketch(
                    str(uuid4()), "Rank reference", p, (), tuple(constraints)
                )
                result = solve(sketch)
                rank = np.linalg.matrix_rank(matrix, tol=1e-10)
                self.assertEqual(result.degrees_of_freedom, 2 * count - rank)
                self.assertEqual(result.redundant_equations, len(matrix) - rank)
                # Independent least-squares projection onto these consistent
                # integer equations. This reference does not use graph traversal.
                A = np.asarray(matrix)
                b = A @ target.flatten()
                seed = np.array([v for row in p for v in row.xy_mm])
                projected = seed + np.linalg.lstsq(A, b - A @ seed, rcond=None)[0]
                actual = np.array([v for row in p for v in result.coordinates[row.id]])
                np.testing.assert_allclose(actual, projected, atol=2e-12, rtol=0)

    def test_bracket_dimension_edit_reorder_and_drag_preserve_identity_and_intent(self):
        sketch = bracket_profile()
        before = asdict(sketch)
        self.assertEqual(solve(sketch).degrees_of_freedom, 0)
        width = next(c for c in sketch.constraints if c.name == "Width")
        changed = settled(
            replace(
                sketch,
                constraints=tuple(
                    replace(c, values_mm=("100",)) if c.id == width.id else c
                    for c in sketch.constraints
                ),
            )
        )
        result = solve(changed)
        expected = ((0, 0), (100, 0), (100, 20), (30, 20), (30, 60), (0, 60))
        self.assertEqual(
            tuple(result.coordinates[p.id] for p in sketch.points), expected
        )
        shuffled = replace(
            changed,
            points=tuple(reversed(changed.points)),
            edges=tuple(reversed(changed.edges)),
            constraints=tuple(reversed(changed.constraints)),
        )
        self.assertEqual(solve(shuffled).coordinates, result.coordinates)
        self.assertEqual(shuffled.canonical(), changed.canonical())
        self.assertEqual(
            Sketch.from_dict(json.loads(json.dumps(asdict(shuffled)))), shuffled
        )
        fixed = settled(changed, drag=(sketch.points[2].id, (200, 200)))
        self.assertEqual(fixed, changed)
        free = replace(
            changed,
            constraints=tuple(c for c in changed.constraints if c.kind != "fixed"),
        )
        self.assertEqual(solve(free).degrees_of_freedom, 2)
        moved = settled(free, drag=(free.points[0].id, (7, -4)))
        np.testing.assert_allclose(
            [p.xy_mm for p in moved.points],
            np.array(expected) + (7, -4),
            atol=0,
            rtol=0,
        )
        self.assertEqual(asdict(sketch), before)

    def test_polygon_admission_rejects_crossings_touches_branches_and_multiple_loops(
        self,
    ):
        good = bracket_profile().polygon()
        area = (
            sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(good, good[1:] + good[:1]))
            / 2
        )
        self.assertEqual(area, 80 * 20 + 30 * 40)
        self.assertEqual(len(good), 6)
        for points in (
            ((0, 0), (10, 10), (0, 10), (10, 0)),
            ((0, 0), (10, 0), (5, 0), (10, 10), (0, 10)),
            ((0, 0), (10, 0), (10, 10), (5, 0), (0, 10)),
            ((0, 0), (1, 0), (2, 0)),
            ((0, 0), (0.0001, 0), (1, 1)),
        ):
            with self.subTest(points=points), self.assertRaises(ValueError):
                profile(points).polygon()
        a = profile(((0, 0), (1, 0), (0, 1)))
        b = profile(((2, 2), (3, 2), (2, 3)))
        with self.assertRaisesRegex(ValueError, "connected"):
            replace(a, points=a.points + b.points, edges=a.edges + b.edges).polygon()
        with self.assertRaisesRegex(ValueError, "closed"):
            replace(a, edges=a.edges[:-1]).polygon()
        # Near-collinear geometry at a large offset uses the exact binary64
        # coordinates delivered to CAD, not a cancellation-prone orientation.
        p = profile(
            (
                (999000.0, 999000.0),
                (999001.0, 999000.0),
                (999002.0, 999000.00001),
                (999000.0, 999001.0),
            )
        )
        self.assertEqual(len(p.polygon()), 4)

    def test_decimal_and_reference_readers_are_bounded_and_reject_unknown_data(self):
        for value in (
            "nan",
            "inf",
            "1e9999",
            "9" * 100,
            "__import__('os')",
            "1/3",
            "1 mm",
            2.0,
            "1000001",
        ):
            with self.subTest(value=value), self.assertRaises(ValueError):
                decimal_mm(value)
        sketch = bracket_profile()
        with self.assertRaisesRegex(ValueError, "missing point"):
            replace(sketch, points=sketch.points[:-1])
        with self.assertRaisesRegex(ValueError, "missing entity"):
            replace(sketch, edges=sketch.edges[:-1])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            replace(sketch, points=sketch.points + (sketch.points[0],))
        with self.assertRaises(ValueError):
            replace(sketch.points[0], xy_mm=(True, 0))
        payload = asdict(sketch)
        payload["unrecognized"] = "data"
        with self.assertRaises(ValueError):
            Sketch.from_dict(payload)
        with self.assertRaises(ValueError):
            SketchEdge(str(uuid4()), "Loop", sketch.points[0].id, sketch.points[0].id)


if __name__ == "__main__":
    unittest.main()
