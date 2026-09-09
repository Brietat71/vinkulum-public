"""Independent qualification of explicit joint frames and temporal loads."""

import math
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from vinkulum import Noyau

I = np.eye(3).reshape(9).tolist()
Z = ("constante", [0.0])


def body():
    n = Noyau(g=[0.0, 0.0, 0.0])
    n.corps("body", 2.0, (2 * np.eye(3)).reshape(9).tolist(), [0.0, 0.0, 0.0])
    return n


class TemporalLoads(unittest.TestCase):
    def test_constant_matches_existing_force(self):
        old, new = body(), body()
        old.effort(0, [1.0, -2.0, 3.0], [0.1, 0.2, -0.3])
        new.effort_temporel(
            "load",
            0,
            [("constante", [x]) for x in (1.0, -2.0, 3.0)],
            [("constante", [x]) for x in (0.1, 0.2, -0.3)],
        )
        a, b = old.simule(0.2, 0.005), new.simule(0.2, 0.005)
        for left, right in zip(a, b):
            for i in range(5):
                np.testing.assert_allclose(left[i], right[i], atol=1e-13, rtol=1e-13)

    def test_linear_force_against_polynomial_solution(self):
        n = body()
        n.effort_temporel("ramp", 0, [("lineaire", [2.0, 6.0]), Z, Z], [Z, Z, Z])
        frames = n.simule(1.0, 0.005)
        errors = [
            abs(row[1][0][0] - (0.5 * row[0] ** 2 + 0.5 * row[0] ** 3))
            for row in frames
        ]
        self.assertLess(max(errors), 3e-5)
        self.assertAlmostEqual(frames[-1][3][0][0], 2.5, delta=3e-5)

    def test_torque_against_axis_rotation(self):
        n = body()
        n.effort_temporel("torque", 0, [Z, Z, Z], [Z, Z, ("constante", [0.4])])
        final = n.simule(1.0, 0.005)[-1]
        self.assertAlmostEqual(
            math.atan2(final[2][0][3], final[2][0][0]), 0.1, delta=1e-7
        )
        self.assertAlmostEqual(final[4][0][2], 0.2, delta=1e-7)

    def test_offset_force_tangent_and_refinement(self):
        def model():
            n = Noyau(g=[0.0, 0.0, 0.0])
            c, s = math.cos(0.3), math.sin(0.3)
            n.corps(
                "offset",
                2.0,
                (2 * np.eye(3)).reshape(9).tolist(),
                [0.0, 0.0, 0.0],
                rot=[c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0],
                w=[0.0, 0.0, 0.2],
            )
            n.effort_temporel(
                "offset",
                0,
                [("lineaire", [1.0, 2.0]), Z, Z],
                [Z, Z, Z],
                point=[0.0, 1.0, 0.0],
            )
            return n

        n = model()
        blocks, _ = n.audit_jacobien(0.02)
        self.assertLess(max(error for _, error, _ in blocks), 2e-5)
        reference = solve_ivp(
            lambda t, y: (y[1], -(1 + 2 * t) * math.cos(y[0]) / 2),
            (0, 0.5),
            [0.3, 0.2],
            method="DOP853",
            rtol=1e-12,
            atol=1e-14,
            dense_output=True,
        )
        errors = []
        for h in (0.02, 0.01, 0.005):
            frames = model().simule(0.5, h)
            errors.append(
                max(
                    abs(
                        math.atan2(row[2][0][3], row[2][0][0])
                        - reference.sol(row[0])[0]
                    )
                    for row in frames
                )
            )
        self.assertLess(errors[-1], 2e-5)
        self.assertGreater(min(errors[i] / errors[i + 1] for i in range(2)), 3.5)

    def test_conservation_and_unsupported_static_are_refused(self):
        for kwargs in ({"moment": True}, {"energie": True}):
            n = body()
            n.effort_temporel("load", 0, [("lineaire", [0.0, 1.0]), Z, Z], [Z, Z, Z])
            with self.assertRaises(ValueError):
                n.simule(0.1, 0.01, **kwargs)
        with self.assertRaisesRegex(ValueError, "statique"):
            n.statique()

    def test_table_values_endpoints_and_rejections(self):
        n = body()
        for law in (
            ("table", []),
            ("table", [0.0, 1.0, 0.0, 2.0]),
            ("constante", [float("nan")]),
        ):
            with self.assertRaises(ValueError):
                n.effort_temporel("invalid", 0, [law, Z, Z], [Z, Z, Z])
        # Constant table includes both clamped endpoints and interpolation.
        n.effort_temporel(
            "table", 0, [("table", [0.1, 2.0, 0.2, 2.0]), Z, Z], [Z, Z, Z]
        )
        final = n.simule(0.3, 0.005)[-1]
        self.assertAlmostEqual(final[1][0][0], 0.045, delta=1e-9)


class ExplicitFrames(unittest.TestCase):
    def model(self, explicit):
        n = Noyau()
        n.corps(
            "pendulum",
            0.2,
            (1e-8 * np.eye(3)).reshape(9).tolist(),
            [0.5, 0.0, -math.sqrt(0.75)],
        )
        if explicit:
            n.liaison_reperes(
                "joint",
                None,
                0,
                [0.0, 0.0, 0.0],
                I,
                [-0.5, 0.0, math.sqrt(0.75)],
                I,
                bloque_t=[0, 1, 2],
                bloque_r=[],
            )
        else:
            n.liaison("joint", None, 0, bloque_t=[0, 1, 2], bloque_r=[])
        return n

    def test_two_frames_reproduce_existing_joint(self):
        a, b = self.model(False).simule(0.2, 0.005), self.model(True).simule(0.2, 0.005)
        for left, right in zip(a, b):
            np.testing.assert_allclose(left[1], right[1], atol=1e-13, rtol=1e-13)

    def test_explicit_frames_do_not_reanchor_or_accept_invalid_rotations(self):
        n = body()
        with self.assertRaises(ValueError):
            n.liaison_reperes(
                "bad", None, 0, [0.0, 0.0, 0.0], I, [0.0, 0.0, 0.0], [0.0] * 9
            )
        n.liaison_reperes("explicit", None, 0, [0.0, 0.0, 0.0], I, [1.0, 0.0, 0.0], I)
        # The independent B point is one metre from A: it must survive creation.
        # Assembly can correct the pose, but creation itself must not do so.
        initial = n.etat()
        self.assertEqual(initial[1][0], [0.0, 0.0, 0.0])
        n.assemble()
        np.testing.assert_allclose(n.etat()[1][0], [-1.0, 0.0, 0.0], atol=1e-9)


if __name__ == "__main__":
    unittest.main()
