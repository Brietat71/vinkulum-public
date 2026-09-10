"""Applied-load derivatives from rotations and virtual work, independently of Pinocchio."""

import math
import json
import runpy
import tempfile
from pathlib import Path
import unittest
from dataclasses import replace
from importlib.util import find_spec

import numpy as np

from vinkulum_studio.document import Body, Law, Load, Project, joint_at, new_id
from vinkulum_studio.pinocchio_backend import ArticulatedModel, run_analysis
from vinkulum_studio.articulated_result import load_operators
from vinkulum_studio.examples3d import ry, rz

REFERENCE = runpy.run_path(
    Path(__file__).parent / "fixtures/applied_loads/reference.py"
)
spatial_pair = REFERENCE["spatial_pair"]
rotation_reference = REFERENCE["rotation_reference"]


@unittest.skipUnless(
    find_spec("pinocchio"), "Install the separate qualified Pinocchio environment"
)
class AppliedLoadReferences(unittest.TestCase):
    def test_issue_six_force_and_nonsymmetric_free_moment(self):
        q, length, force, moment = (0.3, -0.4), 0.7, 5.0, 13.0
        x, y = q
        force_effort = (
            -force
            * length
            * np.array([math.sin(x) * math.cos(y), math.cos(x) * math.sin(y)])
        )
        force_derivative = (
            force
            * length
            * np.array(
                [
                    [-math.cos(x) * math.cos(y), math.sin(x) * math.sin(y)],
                    [math.sin(x) * math.sin(y), -math.cos(x) * math.cos(y)],
                ]
            )
        )
        moment_effort = np.array([0, moment * math.sin(x)])
        moment_derivative = np.array([[0, 0], [moment * math.cos(x), 0]])
        for F, M, effort, derivative in (
            ((0, 0, force), (0, 0, 0), force_effort, force_derivative),
            ((0, 0, 0), (0, 0, moment), moment_effort, moment_derivative),
            (
                (0, 0, force),
                (0, 0, moment),
                force_effort + moment_effort,
                force_derivative + moment_derivative,
            ),
        ):
            with self.subTest(force=F, moment=M):
                project = spatial_pair(
                    centre=(0, 0, 0.2), point=(0, 0, length - 0.2), force=F, moment=M
                )
                result = ArticulatedModel(project).analyse(
                    q=q, velocity=(0.2, 0.4), acceleration=(0.7, -0.3)
                )
                np.testing.assert_allclose(
                    result["external_effort"], effort, rtol=1e-13, atol=2e-14
                )
                np.testing.assert_allclose(
                    result["external_effort_derivatives"]["q"],
                    derivative,
                    rtol=1e-13,
                    atol=2e-14,
                )
                np.testing.assert_allclose(
                    result["loaded_inverse_derivatives"]["q"],
                    np.array(result["intrinsic_inverse_derivatives"]["q"]) - derivative,
                    rtol=1e-13,
                    atol=2e-14,
                )
                for key in ("velocity", "acceleration"):
                    np.testing.assert_array_equal(
                        result["external_effort_derivatives"][key], np.zeros((2, 2))
                    )
                    self.assertEqual(
                        result["loaded_inverse_derivatives"][key],
                        result["intrinsic_inverse_derivatives"][key],
                    )
        self.assertGreater(np.linalg.norm(moment_derivative - moment_derivative.T), 10)

    def test_spatial_and_mixed_pairs_against_explicit_rotations(self):
        rng = np.random.default_rng(62106)
        for mixed in (False, True):
            project = spatial_pair(mixed=mixed)
            model = ArticulatedModel(project)
            load = project.loads[0]
            arm = np.array(project.bodies[1].position) + load.point
            F, M = np.array([2, -3, 5]), np.array([7, 11, 13])
            for q in rng.uniform(-0.8, 0.8, (16, 2)):
                with self.subTest(mixed=mixed, q=q):
                    _, _, _, _, _, effort, derivative = rotation_reference(
                        q, arm, F, M, mixed=mixed
                    )
                    r = model.analyse(
                        q=q.tolist(), velocity=(0.2, -0.1), acceleration=(0.4, 0.3)
                    )
                    np.testing.assert_allclose(
                        r["external_effort"], effort, rtol=1e-13, atol=3e-14
                    )
                    np.testing.assert_allclose(
                        r["external_effort_derivatives"]["q"],
                        derivative,
                        rtol=1e-13,
                        atol=3e-14,
                    )

    def test_world_body_reframing_and_reversed_edges_keep_units_and_signs(self):
        q, v, a = np.array([0.3, -0.4]), np.array([0.2, 0.6]), np.array([0.7, -0.3])
        R = np.array(rz(0.7)).reshape(3, 3) @ np.array(ry(-0.4)).reshape(3, 3)
        S = np.array(ry(0.6)).reshape(3, 3) @ np.array(rz(-0.2)).reshape(3, 3)
        for mixed in (False, True):
            project = spatial_pair(mixed=mixed)
            baseline = ArticulatedModel(project).analyse(
                q=q.tolist(), velocity=v.tolist(), acceleration=a.tolist()
            )
            variants = [
                (project, np.eye(2)),
                (
                    REFERENCE["move_world"](project, R, np.array([1000, -2000, 3000])),
                    np.eye(2),
                ),
                (REFERENCE["reframe_bodies"](project, S), np.eye(2)),
            ]
            for reversed_index in (0, 1):
                joints = list(project.joints)
                j = joints[reversed_index]
                joints[reversed_index] = replace(
                    j, a=j.b, b=j.a, pa=j.pb, pb=j.pa, ra=j.rb, rb=j.ra
                )
                signs = np.eye(2)
                signs[reversed_index, reversed_index] = -1
                variants.append(
                    (
                        replace(
                            project,
                            joints=tuple(joints[::-1]),
                            bodies=project.bodies[::-1],
                        ),
                        signs,
                    )
                )
            for changed, signs in variants:
                with self.subTest(mixed=mixed, signs=signs, project=changed.id):
                    result = ArticulatedModel(changed).analyse(
                        q=(signs @ q).tolist(),
                        velocity=(signs @ v).tolist(),
                        acceleration=(signs @ a).tolist(),
                    )
                    for key in ("external_effort", "inverse_effort"):
                        np.testing.assert_allclose(
                            result[key], signs @ baseline[key], rtol=2e-12, atol=2e-12
                        )
                    for channel in (
                        "intrinsic_inverse_derivatives",
                        "external_effort_derivatives",
                        "loaded_inverse_derivatives",
                    ):
                        for variable in ("q", "velocity", "acceleration"):
                            np.testing.assert_allclose(
                                result[channel][variable],
                                signs @ baseline[channel][variable] @ signs,
                                rtol=2e-12,
                                atol=2e-12,
                            )
                    self.assertEqual(
                        [x["joint_id"] for x in result["coordinate_order"]],
                        [x["joint_id"] for x in baseline["coordinate_order"]],
                    )
                    self.assertEqual(
                        result["coordinate_order"][1]["coordinate_unit"],
                        "m" if mixed else "rad",
                    )
                    self.assertEqual(
                        result["coordinate_order"][1]["effort_unit"],
                        "N" if mixed else "N m",
                    )

    def test_branches_multiple_loads_and_fixed_load_time(self):
        p = spatial_pair()
        other = Body(new_id(), "Separate branch", position=(1, 2, 3))
        p = replace(p, bodies=(*p.bodies, other))
        p = replace(
            p,
            joints=(
                *p.joints,
                joint_at(p, "pivot", None, other.id, point=(1, 2, 3), axis=(0, 0, 1)),
            ),
        )
        load = p.loads[0]
        extra = replace(
            load,
            id=new_id(),
            force=tuple(Law("lineaire", (x, y)) for x, y in ((1, 2), (3, -1), (5, 4))),
            moment=(Law(), Law(), Law("lineaire", (7, -2))),
        )
        p = replace(p, loads=(load, extra))
        model = ArticulatedModel(p)
        active = [i for i, link in enumerate(model.links) if link.child != other.id]
        inactive = next(
            i for i, link in enumerate(model.links) if link.child == other.id
        )
        q = np.zeros(3)
        q[active] = [0.3, -0.4]
        q[inactive] = 0.7
        arm = np.array(p.bodies[1].position) + load.point
        for t in (0, 0.5, 2):
            F = np.sum(
                [[law.value(t) for law in item.force] for item in p.loads], axis=0
            )
            M = np.sum(
                [[law.value(t) for law in item.moment] for item in p.loads], axis=0
            )
            *_, expected_effort, expected = rotation_reference(q[active], arm, F, M)
            result = model.analyse(q=q.tolist(), time_s=t)
            actual = np.array(result["external_effort_derivatives"]["q"])
            np.testing.assert_allclose(
                actual[np.ix_(active, active)], expected, rtol=1e-13, atol=5e-14
            )
            np.testing.assert_allclose(
                np.array(result["external_effort"])[active],
                expected_effort,
                rtol=1e-13,
                atol=5e-14,
            )
            np.testing.assert_array_equal(actual[inactive, :], np.zeros(3))
            np.testing.assert_array_equal(actual[:, inactive], np.zeros(3))
            no_loads = ArticulatedModel(replace(p, loads=())).analyse(
                q=q.tolist(), time_s=t
            )
            self.assertEqual(
                result["intrinsic_inverse_derivatives"],
                no_loads["intrinsic_inverse_derivatives"],
            )
            np.testing.assert_array_equal(
                no_loads["external_effort_derivatives"]["q"], np.zeros((3, 3))
            )

    def test_single_coordinate_shapes_and_prismatic_zero_derivative(self):
        for kind, expected in (
            ("pivot", -3 * math.sin(0.3) - 5 * math.cos(0.3)),
            ("glissiere", 0),
        ):
            body = Body(new_id(), "Single", position=(0, 0, 0.7))
            p = Project(new_id(), bodies=(body,))
            j = joint_at(p, kind, None, body.id, axis=(0, 1, 0))
            load = Load(
                new_id(),
                "Point force",
                body.id,
                point=(0, 0, 0.3),
                force=(Law(values=(3,)), Law(), Law(values=(5,))),
            )
            p = replace(p, joints=(j,), loads=(load,))
            result = ArticulatedModel(p).analyse(q=[0.3])
            np.testing.assert_allclose(
                result["external_effort_derivatives"]["q"],
                [[expected]],
                rtol=1e-13,
                atol=2e-14,
            )

    def test_reader_rejects_loaded_channel_corruption_and_accepts_legacy_archive(self):
        old = (
            Path(__file__).resolve().parents[3]
            / "examples/studio/articulated/double-pendulum"
        )
        legacy = load_operators(old)
        self.assertEqual(legacy.report["schema"], 1)
        self.assertNotIn("external_effort_derivatives", legacy.report)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "result"
            report = run_analysis(
                spatial_pair(mixed=True), directory, {"q": [0.3, -0.4]}
            )
            self.assertEqual(load_operators(directory).report, report)
            edits = (
                lambda r: r.pop("external_effort_derivatives"),
                lambda r: r["external_effort_derivatives"]["q"][0].__setitem__(1, True),
                lambda r: r["external_effort_derivatives"]["q"].pop(),
                lambda r: r["external_effort_derivatives"]["velocity"][0].__setitem__(
                    0, 1
                ),
                lambda r: r["loaded_inverse_derivatives"]["q"][0].__setitem__(0, 99),
                lambda r: r["loaded_inverse_derivatives"]["acceleration"][
                    0
                ].__setitem__(0, 99),
                lambda r: r["conventions"].__setitem__(
                    "external_derivatives", "configuration-dependent follower force"
                ),
                lambda r: r.update(schema=1, adapter_version="0.1.0"),
            )
            for edit in edits:
                altered = json.loads(json.dumps(report))
                edit(altered)
                (directory / "result.json").write_text(json.dumps(altered))
                with (
                    self.subTest(edit=edit),
                    self.assertRaises((ValueError, TypeError)),
                ):
                    load_operators(directory)
            # Even coordinated edits preserving loaded = intrinsic - external
            # must disagree with the captured Jacobians and fixed world loads.
            altered = json.loads(json.dumps(report))
            altered["external_effort_derivatives"]["q"][0][1] += 1
            altered["loaded_inverse_derivatives"]["q"][0][1] -= 1
            (directory / "result.json").write_text(json.dumps(altered))
            with self.assertRaisesRegex(ValueError, "applied world load derivative"):
                load_operators(directory)


if __name__ == "__main__":
    unittest.main()
