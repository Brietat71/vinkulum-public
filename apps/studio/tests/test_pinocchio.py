"""Independent references for the optional Pinocchio document adapter."""

import math
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from importlib.util import find_spec
from pathlib import Path

import numpy as np
from vinkulum_studio.articulated import tree_links
from vinkulum_studio.document import (
    Body,
    Joint,
    Law,
    Load,
    Project,
    joint_at,
    new_id,
    save_project,
)
from vinkulum_studio.examples3d import double_pendulum, ry, rz, slider_crank
from vinkulum_studio.pinocchio_backend import ArticulatedModel, run_analysis


def two_link():
    lengths, masses, width = (0.8, 0.6), (1.2, 0.7), 0.03
    bodies = (
        Body(
            new_id(),
            "Upper",
            dimensions=(width, width, lengths[0]),
            mass=masses[0],
            position=(0, 0, -lengths[0] / 2),
        ),
        Body(
            new_id(),
            "Lower",
            dimensions=(width, width, lengths[1]),
            mass=masses[1],
            position=(0, 0, -lengths[0] - lengths[1] / 2),
        ),
    )
    project = Project(new_id(), bodies=bodies, gravity=(0, 0, -9.81))
    joints = (
        joint_at(project, "pivot", None, bodies[0].id, axis=(0, 1, 0)),
        joint_at(
            project,
            "pivot",
            bodies[0].id,
            bodies[1].id,
            point=(0, 0, -lengths[0]),
            axis=(0, 1, 0),
        ),
    )
    return replace(project, joints=joints)


def lagrange(q, v, a):
    l1, l2, m1, m2, w, g = 0.8, 0.6, 1.2, 0.7, 0.03, 9.81
    c1, c2 = l1 / 2, l2 / 2
    i1, i2 = m1 * (l1 * l1 + w * w) / 12, m2 * (l2 * l2 + w * w) / 12
    k = m2 * l1 * c2
    s, c = math.sin(q[1]), math.cos(q[1])
    B = m1 * c1 + m2 * l1
    mass = np.array(
        [
            [
                i1 + i2 + m1 * c1 * c1 + m2 * (l1 * l1 + c2 * c2) + 2 * k * c,
                i2 + m2 * c2 * c2 + k * c,
            ],
            [i2 + m2 * c2 * c2 + k * c, i2 + m2 * c2 * c2],
        ]
    )
    h = np.array([-k * s * (2 * v[0] * v[1] + v[1] ** 2), k * s * v[0] ** 2])
    h += g * np.array(
        [B * math.sin(q[0]) + m2 * c2 * math.sin(sum(q)), m2 * c2 * math.sin(sum(q))]
    )
    dq = np.column_stack(
        (
            g
            * np.array(
                [
                    B * math.cos(q[0]) + m2 * c2 * math.cos(sum(q)),
                    m2 * c2 * math.cos(sum(q)),
                ]
            ),
            -k * s * np.array([[2, 1], [1, 0]]) @ a
            + np.array([-k * c * (2 * v[0] * v[1] + v[1] ** 2), k * c * v[0] ** 2])
            + g * m2 * c2 * math.cos(sum(q)),
        )
    )
    dv = np.array([[-2 * k * s * v[1], -2 * k * s * sum(v)], [2 * k * s * v[0], 0]])
    return mass, h, dq, dv


class TreeAdmission(unittest.TestCase):
    def test_loops_detached_unsupported_motion_and_bad_anchors_are_rejected(self):
        original = two_link()
        self.assertEqual(len(tree_links(original)), 2)
        for project in (
            replace(original, joints=original.joints[:1]),
            replace(
                original,
                joints=(*original.joints, replace(original.joints[0], id=new_id())),
            ),
            replace(
                original,
                joints=(replace(original.joints[0], kind="rotule"), original.joints[1]),
            ),
            replace(
                original,
                joints=(replace(original.joints[0], motion=Law()), original.joints[1]),
            ),
            replace(
                original,
                joints=(
                    replace(original.joints[0], pa=(0.01, 0, 0)),
                    original.joints[1],
                ),
            ),
            slider_crank(),
        ):
            with self.subTest(joints=project.joints), self.assertRaises(ValueError):
                tree_links(project)


@unittest.skipUnless(
    find_spec("pinocchio"), "Install the separate qualified Pinocchio environment"
)
class PinocchioReferences(unittest.TestCase):
    def test_ground_on_side_b_and_independent_branches(self):
        project = two_link()
        first, second = project.bodies
        root = project.joints[0]
        reversed_root = replace(
            root, a=first.id, b=None, pa=root.pb, pb=root.pa, ra=root.rb, rb=root.ra
        )
        other_root = joint_at(
            project, "pivot", None, second.id, point=(0, 0, -0.8), axis=(0, 1, 0)
        )
        model = ArticulatedModel(replace(project, joints=(other_root, reversed_root)))
        q = [0.3, -0.5]
        result = model.analyse(q=q)
        bodies = {body.id: body for body in project.bodies}
        expected_mass, expected_bias = [], []
        for i, link in enumerate(model.links):
            body = bodies[link.child]
            radius = body.dimensions[2] / 2
            expected_mass.append(body.inertia()[4] + body.mass * radius**2)
            expected_bias.append(body.mass * 9.81 * radius * math.sin(q[i]))
        np.testing.assert_allclose(
            result["mass_matrix"], np.diag(expected_mass), atol=1e-14
        )
        np.testing.assert_allclose(result["intrinsic_bias"], expected_bias, atol=1e-13)
        first_values = next(body for body in result["bodies"] if body["id"] == first.id)
        index = next(i for i, link in enumerate(model.links) if link.child == first.id)
        np.testing.assert_allclose(
            np.array(first_values["jacobian"])[3:, index], [0, -1, 0], atol=1e-14
        )

    def test_one_link_absolute_pose_full_inertia_and_scalar_matrix_shapes(self):
        angle, mass, length = 0.35, 1.2, 0.5
        S = np.array(rz(0.7)).reshape(3, 3) @ np.array(ry(0.2)).reshape(3, 3)
        inertia = S @ np.diag([0.03, 0.04, 0.05]) @ S.T
        R = np.array(ry(angle)).reshape(3, 3)
        body = Body(
            new_id(),
            "Offset arm",
            mass=mass,
            position=tuple(R @ np.array([0, 0, -length])),
            orientation=tuple(R.flat),
            inertia_mode="explicit",
            explicit_inertia=tuple(inertia.flat),
        )
        frame = (0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
        joint = Joint(
            new_id(),
            "Hinge",
            "pivot",
            None,
            body.id,
            pb=(0, 0, length),
            ra=frame,
            rb=frame,
        )
        load = Load(
            new_id(),
            "Wrench",
            body.id,
            point=(0.15, 0.2, -0.1),
            force=tuple(Law(values=(x,)) for x in (3, 4, 5)),
            moment=tuple(Law(values=(x,)) for x in (7, 11, 13)),
        )
        project = Project(
            new_id(),
            bodies=(body,),
            joints=(joint,),
            loads=(load,),
            gravity=(0, 0, -9.81),
        )
        r = ArticulatedModel(project).analyse(velocity=[0.6], acceleration=[0.8])
        self.assertAlmostEqual(r["state"]["q"][0], angle, places=14)
        M = inertia[1, 1] + mass * length**2
        gravity = mass * 9.81 * length * math.sin(angle)
        arm = R @ np.array([0.15, 0.2, -length - 0.1])
        external = np.cross(arm, [3, 4, 5])[1] + 11
        np.testing.assert_allclose(r["mass_matrix"], [[M]], rtol=1e-13, atol=1e-14)
        np.testing.assert_allclose(r["external_effort"], [external], atol=1e-13)
        np.testing.assert_allclose(
            r["inverse_effort"], [M * 0.8 + gravity - external], atol=1e-13
        )
        np.testing.assert_allclose(
            r["forward_acceleration"], [(external - gravity) / M], atol=1e-12
        )
        self.assertEqual(np.array(r["bodies"][0]["jacobian"]).shape, (6, 1))
        np.testing.assert_allclose(
            r["intrinsic_inverse_derivatives"]["q"],
            [[mass * 9.81 * length * math.cos(angle)]],
            atol=1e-13,
        )

    def test_two_link_lagrange_operators_derivatives_and_kinematics(self):
        project = two_link()
        model = ArticulatedModel(project)
        rng = np.random.default_rng(408)
        for _ in range(64):
            q, v, a = (
                rng.uniform(-2.5, 2.5, 2),
                rng.uniform(-2, 2, 2),
                rng.uniform(-4, 4, 2),
            )
            M, h, dq, dv = lagrange(q, v, a)
            tau = M @ a + h
            r = model.analyse(
                q=q.tolist(),
                velocity=v.tolist(),
                acceleration=a.tolist(),
                effort=tau.tolist(),
            )
            np.testing.assert_allclose(r["mass_matrix"], M, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(r["intrinsic_bias"], h, rtol=1e-12, atol=1e-13)
            np.testing.assert_allclose(r["inverse_effort"], tau, rtol=1e-12, atol=1e-13)
            np.testing.assert_allclose(
                r["forward_acceleration"], a, rtol=1e-11, atol=1e-12
            )
            np.testing.assert_allclose(
                r["intrinsic_inverse_derivatives"]["q"], dq, rtol=1e-12, atol=1e-13
            )
            np.testing.assert_allclose(
                r["intrinsic_inverse_derivatives"]["velocity"],
                dv,
                rtol=1e-12,
                atol=1e-13,
            )
            np.testing.assert_allclose(
                r["intrinsic_inverse_derivatives"]["acceleration"],
                M,
                rtol=1e-12,
                atol=1e-14,
            )
            self.assertAlmostEqual(r["kinetic_energy_J"], 0.5 * v @ M @ v, places=12)
            bodies = {body["id"]: body for body in r["bodies"]}
            angle = q.sum()
            pos = (
                -0.8 * math.sin(q[0]) - 0.3 * math.sin(angle),
                0,
                -0.8 * math.cos(q[0]) - 0.3 * math.cos(angle),
            )
            J = np.array(
                [
                    [
                        -0.8 * math.cos(q[0]) - 0.3 * math.cos(angle),
                        -0.3 * math.cos(angle),
                    ],
                    [0, 0],
                    [
                        0.8 * math.sin(q[0]) + 0.3 * math.sin(angle),
                        0.3 * math.sin(angle),
                    ],
                    [0, 0],
                    [1, 1],
                    [0, 0],
                ]
            )
            lower = bodies[project.bodies[1].id]
            np.testing.assert_allclose(lower["position_m"], pos, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(
                lower["body_to_world"],
                np.array(ry(angle)).reshape(3, 3),
                rtol=1e-12,
                atol=1e-14,
            )
            np.testing.assert_allclose(lower["jacobian"], J, rtol=1e-12, atol=1e-14)
            np.testing.assert_allclose(
                lower["velocity_world"], J @ v, rtol=1e-12, atol=1e-14
            )

    def test_reversed_edges_preserve_declared_sign_and_document_order_independence(
        self,
    ):
        project = two_link()
        q, v, a, effort = (
            np.array([0.4, -0.7]),
            np.array([0.2, 0.6]),
            np.array([-0.2, 0.8]),
            np.array([0.9, -0.3]),
        )
        base = ArticulatedModel(project).analyse(
            q=q.tolist(),
            velocity=v.tolist(),
            acceleration=a.tolist(),
            effort=effort.tolist(),
        )
        joint = project.joints[1]
        reversed_joint = replace(
            joint,
            a=joint.b,
            b=joint.a,
            pa=joint.pb,
            pb=joint.pa,
            ra=joint.rb,
            rb=joint.ra,
        )
        reversed_project = replace(
            project,
            bodies=project.bodies[::-1],
            joints=(reversed_joint, project.joints[0]),
        )
        S = np.diag([1, -1])
        r = ArticulatedModel(reversed_project).analyse(
            q=(S @ q).tolist(),
            velocity=(S @ v).tolist(),
            acceleration=(S @ a).tolist(),
            effort=(S @ effort).tolist(),
        )
        np.testing.assert_allclose(
            r["mass_matrix"], S @ base["mass_matrix"] @ S, rtol=1e-13, atol=1e-14
        )
        for key in ("inverse_effort", "forward_acceleration", "intrinsic_bias"):
            np.testing.assert_allclose(r[key], S @ base[key], rtol=1e-13, atol=1e-13)
        original = {b["id"]: b for b in base["bodies"]}
        for b in r["bodies"]:
            np.testing.assert_allclose(
                b["position_m"], original[b["id"]]["position_m"], atol=1e-14
            )
            np.testing.assert_allclose(
                b["jacobian"], np.array(original[b["id"]]["jacobian"]) @ S, atol=1e-14
            )

    def test_world_rotation_translation_and_full_body_inertia_reframing(self):
        project = two_link()
        q, v = [0.5, -0.2], [0.7, 0.3]
        original = ArticulatedModel(project).analyse(q=q, velocity=v)
        R = np.array(rz(0.7)).reshape(3, 3) @ np.array(ry(-0.4)).reshape(3, 3)
        t = np.array([1000, -2000, 3000])
        joints = tuple(
            replace(
                j,
                pa=tuple(R @ j.pa + t),
                ra=tuple((R @ np.array(j.ra).reshape(3, 3)).flat),
            )
            if j.a is None
            else j
            for j in project.joints
        )
        moved = replace(
            project,
            bodies=tuple(
                replace(
                    b,
                    position=tuple(R @ b.position + t),
                    orientation=tuple((R @ np.array(b.orientation).reshape(3, 3)).flat),
                )
                for b in project.bodies
            ),
            joints=joints,
            gravity=tuple(R @ project.gravity),
        )
        r = ArticulatedModel(moved).analyse(q=q, velocity=v)
        np.testing.assert_allclose(
            r["mass_matrix"], original["mass_matrix"], rtol=1e-12, atol=1e-13
        )
        np.testing.assert_allclose(
            r["inverse_effort"], original["inverse_effort"], rtol=1e-12, atol=1e-13
        )
        for first, second in zip(original["bodies"], r["bodies"]):
            np.testing.assert_allclose(
                second["position_m"], R @ first["position_m"] + t, rtol=0, atol=2e-12
            )
            J = np.array(first["jacobian"])
            np.testing.assert_allclose(
                second["jacobian"],
                np.vstack((R @ J[:3], R @ J[3:])),
                rtol=1e-12,
                atol=1e-13,
            )
        self.assertAlmostEqual(
            r["potential_energy_J"] - original["potential_energy_J"],
            -sum(b.mass for b in project.bodies) * np.array(moved.gravity) @ t,
            places=8,
        )
        S = np.array(rz(0.6)).reshape(3, 3) @ np.array(ry(0.4)).reshape(3, 3)
        reframed_bodies = tuple(
            replace(
                b,
                orientation=tuple(S.flat),
                inertia_mode="explicit",
                explicit_inertia=tuple(
                    (S.T @ np.array(b.inertia()).reshape(3, 3) @ S).flat
                ),
            )
            for b in project.bodies
        )
        reframed_joints = tuple(
            replace(
                j,
                pa=j.pa if j.a is None else tuple(S.T @ j.pa),
                ra=j.ra
                if j.a is None
                else tuple((S.T @ np.array(j.ra).reshape(3, 3)).flat),
                pb=tuple(S.T @ j.pb),
                rb=tuple((S.T @ np.array(j.rb).reshape(3, 3)).flat),
            )
            for j in project.joints
        )
        reframe = ArticulatedModel(
            replace(project, bodies=reframed_bodies, joints=reframed_joints)
        ).analyse(q=q, velocity=v)
        np.testing.assert_allclose(
            reframe["mass_matrix"], original["mass_matrix"], rtol=1e-12, atol=1e-13
        )
        np.testing.assert_allclose(
            reframe["inverse_effort"],
            original["inverse_effort"],
            rtol=1e-12,
            atol=1e-13,
        )

    def test_prismatic_initial_coordinate_and_world_point_load(self):
        body = Body(new_id(), "Slider", mass=2, position=(0, 0, 0.4))
        joint = Joint(new_id(), "Guide", "glissiere", None, body.id)
        force = Load(
            new_id(),
            "World load",
            body.id,
            point=(0.1, 0, 0),
            force=(Law(), Law(), Law("lineaire", (2, 4))),
            moment=(Law(values=(1,)), Law(values=(2,)), Law(values=(3,))),
        )
        project = Project(
            new_id(),
            bodies=(body,),
            joints=(joint,),
            loads=(force,),
            gravity=(0, 0, -9.81),
        )
        r = ArticulatedModel(project).analyse(acceleration=[1.5], effort=[4], time_s=1)
        self.assertEqual(r["state"]["q"], [0.4])
        np.testing.assert_allclose(r["mass_matrix"], [[2]], atol=1e-14)
        np.testing.assert_allclose(r["external_effort"], [6], atol=1e-14)
        np.testing.assert_allclose(r["inverse_effort"], [16.62], atol=1e-13)
        np.testing.assert_allclose(r["forward_acceleration"], [-4.81], atol=1e-13)
        np.testing.assert_allclose(
            r["bodies"][0]["position_m"], [0, 0, 0.4], atol=1e-14
        )
        self.assertEqual(r["coordinate_order"][0]["effort_unit"], "N")

    def test_world_load_virtual_work_and_intrinsic_derivatives_remain_explicit(self):
        project = two_link()
        body = project.bodies[1]
        load = Load(
            new_id(),
            "Offset wrench",
            body.id,
            point=(0.1, 0.2, -0.3),
            force=(Law(values=(2,)), Law(values=(3,)), Law(values=(5,))),
            moment=(Law(values=(7,)), Law(values=(11,)), Law(values=(13,))),
        )
        project = replace(project, loads=(load,))
        model = ArticulatedModel(project)
        q, v, a = [0.2, 0.7], [0.4, -0.1], [0.6, 0.8]
        result = model.analyse(q=q, velocity=v, acceleration=a)
        base = ArticulatedModel(replace(project, loads=())).analyse(
            q=q, velocity=v, acceleration=a
        )
        # Differentiate an independently evaluated body point position, rather
        # than reusing the adapter's wrench/Jacobian assembly.
        F = np.array([2, 3, 5])
        M = np.array([7, 11, 13])

        def point(x):
            angle = x[0] + x[1]
            centre = np.array(
                [
                    -0.8 * math.sin(x[0]) - 0.3 * math.sin(angle),
                    0,
                    -0.8 * math.cos(x[0]) - 0.3 * math.cos(angle),
                ]
            )
            return centre + np.array(ry(angle)).reshape(3, 3) @ load.point

        expected = []
        for i in range(2):
            step = np.eye(2)[i] * 1e-5
            expected.append(
                F @ (point(np.array(q) + step) - point(np.array(q) - step)) / (2e-5)
                + M[1]
            )
        np.testing.assert_allclose(
            result["external_effort"], expected, rtol=1e-10, atol=1e-10
        )
        np.testing.assert_allclose(
            result["inverse_effort"],
            np.array(base["inverse_effort"]) - expected,
            rtol=1e-10,
            atol=1e-10,
        )
        self.assertEqual(
            result["intrinsic_inverse_derivatives"],
            base["intrinsic_inverse_derivatives"],
        )
        self.assertIn(
            "Excludes applied world loads", result["conventions"]["derivatives"]
        )

    def test_mixed_revolute_prismatic_lagrange_reference(self):
        m1, m2, l1 = 1.2, 0.7, 0.8
        first = Body(
            new_id(),
            "Arm",
            mass=m1,
            dimensions=(0.03, 0.03, l1),
            position=(0, 0, -l1 / 2),
        )
        second = Body(
            new_id(),
            "Radial slider",
            shape="sphere",
            dimensions=(0.05,),
            mass=m2,
            position=(0, 0, -1.1),
        )
        project = Project(new_id(), bodies=(first, second), gravity=(0, 0, -9.81))
        hinge = joint_at(project, "pivot", None, first.id, axis=(0, 1, 0))
        frame = (1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)
        slide = Joint(
            new_id(),
            "Radial guide",
            "glissiere",
            first.id,
            second.id,
            pa=(0, 0, l1 / 2),
            ra=frame,
            rb=frame,
        )
        model = ArticulatedModel(replace(project, joints=(slide, hinge)))
        self.assertAlmostEqual(model.initial_q[1], 1.1)
        q, v, a = [0.4, 1.3], [0.6, -0.3], [0.2, 0.7]
        M = np.diag(
            [
                first.inertia()[4]
                + m1 * (l1 / 2) ** 2
                + second.inertia()[4]
                + m2 * q[1] ** 2,
                m2,
            ]
        )
        h = np.array(
            [
                2 * m2 * q[1] * v[0] * v[1]
                + 9.81 * (m1 * l1 / 2 + m2 * q[1]) * math.sin(q[0]),
                -m2 * q[1] * v[0] ** 2 - m2 * 9.81 * math.cos(q[0]),
            ]
        )
        r = model.analyse(q=q, velocity=v, acceleration=a, effort=(M @ a + h).tolist())
        np.testing.assert_allclose(r["mass_matrix"], M, atol=1e-14)
        np.testing.assert_allclose(r["intrinsic_bias"], h, atol=1e-13)
        np.testing.assert_allclose(r["forward_acceleration"], a, atol=1e-13)

    def test_snapshot_files_reject_overwrite_and_invalid_state(self):
        project = double_pendulum()
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "operators"
            report = run_analysis(project, directory, {"velocity": [0.1, 0.2]})
            self.assertEqual(report["scientific_status"], "NotAssessed")
            self.assertTrue((directory / "project.vinkulum.json").is_file())
            before = (directory / "result.json").read_bytes()
            with self.assertRaises(FileExistsError):
                run_analysis(project, directory)
            self.assertEqual((directory / "result.json").read_bytes(), before)
            model = ArticulatedModel(project)
            for state in (
                {"q": [1]},
                {"q": [True, 0]},
                {"velocity": [0, float("nan")]},
                {"time_s": -1},
                {"time_s": True},
            ):
                with self.subTest(state=state), self.assertRaises(ValueError):
                    model.analyse(**state)

    def test_cli_runs_outside_source_and_keeps_graphics_out_of_worker_imports(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            project = root / "captured mechanism.json"
            save_project(project, two_link())
            output = root / "result with spaces"
            child = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "vinkulum_studio.pinocchio_backend",
                    str(project),
                    "--output",
                    str(output),
                ],
                cwd=root,
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            self.assertEqual(child.returncode, 0, child.stderr)
            self.assertTrue((output / "result.json").is_file())
            imports = subprocess.run(
                [
                    sys.executable,
                    "-c",
                    "import sys; from vinkulum_studio.pinocchio_backend import ArticulatedModel; from vinkulum_studio.examples3d import double_pendulum; ArticulatedModel(double_pendulum()).analyse(); assert not any(name in sys.modules for name in ('PySide6', 'vtk', 'OCP', 'vinkulum'))",
                ],
                cwd=root,
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            self.assertEqual(imports.returncode, 0, imports.stderr)


if __name__ == "__main__":
    unittest.main()
