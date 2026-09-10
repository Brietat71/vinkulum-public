"""Independent geometry and sign references for CAD-05 / UI-06 direction cues."""

import math
import unittest
from dataclasses import asdict
from types import SimpleNamespace
from uuid import uuid4

import numpy as np
from vinkulum_studio.boundary_display import (
    SurfaceFrame,
    boundary_frames,
    boundary_glyphs,
    effective_values,
    unit_direction,
)
from vinkulum_studio.mesh_binding import MeshCondition


def condition(kind, values=(), axes=(), surfaces=(1,)):
    return MeshCondition(str(uuid4()), kind, kind, surfaces, axes, values)


class BoundaryDirections(unittest.TestCase):
    def test_curved_map_normal_matches_analytic_surface_and_world_frame(self):
        # x=r, y=s, z=3/4*r*s. At the barycentre the position is
        # (1/3,1/3,1/12), and the outward normal is proportional to (-1/4,-1/4,1).
        nodes = np.array(
            (
                (0, 0, 0),
                (1, 0, 0),
                (0, 1, 0),
                (0.5, 0, 0),
                (0.5, 0.5, 3 / 16),
                (0, 0.5, 0),
            )
        )
        normal = np.array((-1 / 4, -1 / 4, 1)) / math.sqrt(9 / 8)
        point = np.array((1 / 3, 1 / 3, 1 / 12))
        R = np.array(((0, -1, 0), (1, 0, 0), (0, 0, 1)))
        for scale, rotation, shift in (
            (1, np.eye(3), np.zeros(3)),
            (1, R, np.array((1e5, -2e5, 3e5))),
            (1e-200, R, np.zeros(3)),
            (1e200, R, np.zeros(3)),
        ):
            with self.subTest(scale=scale, shift=shift):
                world = scale * nodes @ rotation.T + shift
                solid = SimpleNamespace(
                    mesh=SimpleNamespace(nodes=tuple(map(tuple, world))),
                    surfaces=(SimpleNamespace(id=1, triangles=((1, 2, 3, 4, 5, 6),)),),
                )
                frame = boundary_frames(solid)[1][0]
                np.testing.assert_allclose(
                    frame.normal, rotation @ normal, rtol=1e-14, atol=1e-14
                )
                np.testing.assert_allclose(
                    frame.point, scale * rotation @ point + shift, rtol=1e-14, atol=0
                )

    def test_normalization_keeps_tiny_and_large_nonzero_force_directions(self):
        expected = np.array((1, -2, 3)) / math.sqrt(14)
        for scale in (1e-300, 1, 1e300):
            np.testing.assert_allclose(
                unit_direction(scale * np.array((1, -2, 3))), expected, atol=1e-15
            )
        self.assertIsNone(unit_direction((0, 0, 0)))
        self.assertIsNone(unit_direction((float("inf"), 1, 0)))

    def test_pressure_force_and_supports_keep_distinct_sign_and_axis_conventions(self):
        frame = SurfaceFrame((1, 2, 3), (0, 0, 1))
        conditions = (
            condition("pressure", (20.0,)),
            condition("total_force", (3.0, 4.0, 0.0)),
            condition("support", axes=(1, 3)),
        )
        before = tuple(asdict(c) for c in conditions)
        for factor in (2.0, -2.0, 0.0):
            with self.subTest(factor=factor):
                glyphs, invalid = boundary_glyphs({1: (frame,)}, conditions, factor)
                self.assertFalse(invalid)
                kinds = {
                    kind: [g for g in glyphs if g.kind == kind]
                    for kind in ("pressure", "total_force", "support")
                }
                self.assertEqual(
                    [g.direction for g in kinds["support"]], [(1, 0, 0), (0, 0, 1)]
                )
                if factor:
                    sign = math.copysign(1, factor)
                    np.testing.assert_array_equal(
                        kinds["pressure"][0].direction, (0, 0, -sign)
                    )
                    np.testing.assert_allclose(
                        kinds["total_force"][0].direction, (sign * 0.6, sign * 0.8, 0)
                    )
                else:
                    self.assertFalse(kinds["pressure"] or kinds["total_force"])
                self.assertTrue(
                    all(g.point == frame.point and g.surface == 1 for g in glyphs)
                )
        self.assertEqual(tuple(asdict(c) for c in conditions), before)
        self.assertEqual(effective_values(conditions[0], -2), (-40.0,))

    def test_invalid_values_hide_only_load_symbols_and_spatial_sampling_is_bounded(
        self,
    ):
        frames = {
            1: tuple(SurfaceFrame((i / 100, 0, 0), (0, 0, 1)) for i in range(100))
        }
        pressure = condition("pressure", (1e300,))
        support = condition("support", axes=(2,))
        for factor in (None, float("nan"), 1e300):
            glyphs, invalid = boundary_glyphs(frames, (pressure, support), factor)
            self.assertEqual(invalid, (pressure.id,))
            self.assertEqual([g.kind for g in glyphs], ["support"])
        first, _ = boundary_glyphs(frames, (pressure,), 1, per_condition=8)
        second, _ = boundary_glyphs(frames, (pressure,), 1, per_condition=8)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 8)
        self.assertEqual(len({g.point for g in first}), 8)
        self.assertGreater(
            max(g.point[0] for g in first) - min(g.point[0] for g in first), 0.9
        )


if __name__ == "__main__":
    unittest.main()
