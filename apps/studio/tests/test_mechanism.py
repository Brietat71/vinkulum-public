import hashlib
import io
import json
import math
import tempfile
import unittest
import zipfile
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from vinkulum_studio.document import (
    Body,
    Law,
    Load,
    Project,
    joint_at,
    new_id,
    pendulum,
)
from vinkulum_studio.mechanism import MechanicalResult, build_solver, simulate_project


class MechanismAdapter(unittest.TestCase):
    def simulate(self, project):
        with tempfile.TemporaryDirectory() as directory:
            metadata = simulate_project(project, "test", directory)
            return MechanicalResult.read(metadata, "test", project, directory)

    def test_pendulum_adapter_and_immutable_arrays(self):
        from vinkulum_studio.model import Parameters
        from vinkulum_studio.worker import simulate

        params = Parameters(duration=0.2)
        project = pendulum(params)
        result = self.simulate(project)
        reference = simulate(params, "legacy")["samples"]
        np.testing.assert_allclose(
            result.position[:, 0, 0], [r[1] for r in reference], atol=1e-12
        )
        np.testing.assert_allclose(
            result.position[:, 0, 2], [r[2] for r in reference], atol=1e-12
        )
        with self.assertRaises(ValueError):
            result.position.setflags(write=True)
        with self.assertRaises(ValueError):
            result.position[0, 0, 0] = 1

    def test_commanded_slider_initial_velocity_and_motion(self):
        b = Body(new_id(), "Slider")
        p = Project(new_id(), bodies=(b,), gravity=(0.0, 0.0, 0.0), duration=0.2)
        j = joint_at(p, "glissiere", None, b.id, axis=(1.0, 0.0, 0.0))
        p = p.replace_object(replace(j, motion=Law("lineaire", (0.0, 0.3))))
        result = self.simulate(p)
        np.testing.assert_allclose(
            result.position[:, 0, 0], 0.3 * result.time, atol=1e-10
        )
        self.assertAlmostEqual(result.velocity[0, 0, 0], 0.3, delta=1e-10)
        np.testing.assert_allclose(
            result.joint_coordinate[:, 0], 0.3 * result.time, atol=1e-10
        )
        self.assertFalse(p.diagnostics())

    def test_temporal_load_and_export(self):
        b = Body(new_id(), "Body", mass=2.0)
        p = Project(new_id(), bodies=(b,), gravity=(0.0, 0.0, 0.0), duration=0.2)
        p = p.replace_object(
            Load(
                new_id(),
                "Ramp",
                b.id,
                force=(Law("lineaire", (0.0, 6.0)), Law(), Law()),
            )
        )
        result = self.simulate(p)
        np.testing.assert_allclose(
            result.position[:, 0, 0], 0.5 * result.time**3, atol=3e-5
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "result.csv"
            result.export_csv(path)
            metadata = json.loads(path.read_text().splitlines()[0][2:])
            self.assertEqual(metadata["project"]["bodies"][0]["id"], b.id)
            self.assertIn("m/s", path.read_text().splitlines()[1])

    def test_invalid_anchores_refused_before_build(self):
        p = pendulum()
        p = p.replace_object(replace(p.bodies[0], position=(0.0, 0.0, -2.0)))
        with self.assertRaises(ValueError):
            build_solver(p)

    def test_result_identity_hash_and_archive_shape_guards(self):
        p = replace(pendulum(), duration=0.1)
        with tempfile.TemporaryDirectory() as directory:
            data = simulate_project(p, "correct", directory)
            with self.assertRaises(ValueError):
                MechanicalResult.read(data, "wrong", p, directory)
            archive = Path(directory) / "trajectory.npz"
            content = archive.read_bytes()
            archive.write_bytes(content[:-20])
            with self.assertRaises(ValueError):
                MechanicalResult.read(data, "correct", p, directory)
            # A matching hash is insufficient: reject an allocation-bomb header
            # before numpy.load sees it, even though the ZIP payload is tiny.
            with zipfile.ZipFile(io.BytesIO(content)) as original:
                members = {name: original.read(name) for name in original.namelist()}
            header = io.BytesIO()
            np.lib.format.write_array_header_1_0(
                header,
                {"descr": "<f8", "fortran_order": False, "shape": (10**12,)},
            )
            members["time.npy"] = header.getvalue()
            with zipfile.ZipFile(archive, "w") as output:
                for name, payload in members.items():
                    output.writestr(name, payload)
            data["archive_sha256"] = hashlib.sha256(archive.read_bytes()).hexdigest()
            with self.assertRaisesRegex(ValueError, "declared"):
                MechanicalResult.read(data, "correct", p, directory)

    def test_double_pendulum_against_independent_lagrange_equations(self):
        from vinkulum_studio.examples3d import double_pendulum

        p = replace(double_pendulum(), duration=0.5)
        # Homogeneous rods of length 1, mass 1, square section .06 m.
        # Angles are absolute, measured from the downward vertical.
        J = (1 + 0.06**2) / 12
        A = 1.25 + J
        B = 0.5
        D = 0.25 + J
        g = 9.80665

        def rhs(t, y):
            a, b, wa, wb = y
            delta = a - b
            M = np.array([[A, B * math.cos(delta)], [B * math.cos(delta), D]])
            force = np.array(
                [
                    -B * math.sin(delta) * wb**2 - 1.5 * g * math.sin(a),
                    B * math.sin(delta) * wa**2 - 0.5 * g * math.sin(b),
                ]
            )
            acceleration = np.linalg.solve(M, force)
            return wa, wb, *acceleration

        ref = solve_ivp(
            rhs,
            (0, 0.5),
            [math.radians(30), math.radians(-15), 0.0, 0.0],
            method="DOP853",
            rtol=1e-12,
            atol=1e-14,
            dense_output=True,
        )
        errors = []
        for h in (0.01, 0.005, 0.0025):
            result = self.simulate(replace(p, step=h))
            angles = np.arctan2(-result.rotation[:, :, 2], result.rotation[:, :, 8])
            errors.append(float(np.max(np.abs(angles - ref.sol(result.time)[:2].T))))
        print("\nDouble pendulum maximum angle errors:", errors, flush=True)
        self.assertLess(errors[-1], 4e-4)
        self.assertGreater(min(errors[i] / errors[i + 1] for i in range(2)), 3.5)

    def test_slider_crank_against_geometric_closure(self):
        from vinkulum_studio.examples3d import slider_crank

        result = self.simulate(slider_crank())
        angle = 0.5 + 2 * result.time
        reference = 0.3 * np.cos(angle) + np.sqrt(0.8**2 - (0.3 * np.sin(angle)) ** 2)
        np.testing.assert_allclose(
            result.position[:, 2, 0], reference, atol=1e-9, rtol=0.0
        )


if __name__ == "__main__":
    unittest.main()
