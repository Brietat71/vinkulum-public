import json
import tempfile
import unittest
from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import numpy as np
from vinkulum_studio.document import (
    Body,
    History,
    Law,
    Load,
    Project,
    import_g0,
    joint_at,
    load_project,
    new_id,
    pendulum,
    save_project,
)
from vinkulum_studio.model import Parameters, save_parameters


class MechanicalDocument(unittest.TestCase):
    def test_laws_and_invalid_tables(self):
        self.assertEqual(Law("lineaire", (2.0, 3.0)).value(4), 14)
        table = Law("table", (0.0, 2.0, 1.0, 4.0, 2.0, 0.0))
        self.assertEqual([table.value(t) for t in (-1, 0.5, 1.5, 3)], [2, 3, 2, 0])
        for values in ((), (0.0,), (0.0, 1.0, 0.0, 2.0), (float("nan"), 1.0)):
            with self.subTest(values=values), self.assertRaises(ValueError):
                Law("table", values)

    def test_homogeneous_inertias_and_physical_custom_tensor(self):
        b = Body(new_id(), "Box", dimensions=(2.0, 3.0, 4.0), mass=12.0)
        np.testing.assert_allclose(
            np.diag(np.array(b.inertia()).reshape(3, 3)), (25, 20, 13)
        )
        s = replace(b, shape="sphere", dimensions=(2.0,), mass=5.0)
        np.testing.assert_allclose(s.inertia(), np.diag([8.0, 8.0, 8.0]).flat)
        c = replace(b, shape="cylinder", dimensions=(2.0, 3.0))
        np.testing.assert_allclose(c.inertia(), np.diag([21.0, 21.0, 24.0]).flat)
        with self.assertRaises(ValueError):
            replace(
                b,
                inertia_mode="explicit",
                explicit_inertia=(10.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
            )

    def test_inputs_are_deeply_immutable(self):
        position = [1.0, 2.0, 3.0]
        b = Body(new_id(), "Bodies", position=position)
        position[0] = 100
        self.assertEqual(b.position, (1.0, 2.0, 3.0))
        with self.assertRaises(FrozenInstanceError):
            b.name = "modifié"

    def test_rename_reorder_and_delete_preserve_references(self):
        p = pendulum()
        b = p.bodies[0]
        second = Body(new_id(), "Second")
        p = p.replace_object(second).replace_object(replace(b, name="Renommé"))
        p = replace(p, bodies=tuple(reversed(p.bodies)))
        self.assertFalse(p.diagnostics())
        self.assertEqual(p.joints[0].b, b.id)
        deleted = p.remove(b.id)
        self.assertEqual(deleted.joints, p.joints)
        self.assertIn("deleted", deleted.diagnostics()[0].message)

    def test_move_invalidates_joint_undo_restores_with_new_revision(self):
        p = pendulum()
        history = History(p)
        history.commit(
            p.replace_object(replace(p.bodies[0], position=(0.0, 0.0, -2.0)))
        )
        self.assertTrue(history.current.diagnostics())
        self.assertFalse(p.diagnostics())
        history.undo()
        self.assertEqual(history.current.bodies, p.bodies)
        self.assertFalse(history.current.diagnostics())
        self.assertEqual(history.current.revision, 2)
        history.redo()
        self.assertTrue(history.current.diagnostics())
        self.assertEqual(history.current.revision, 3)

    def test_joint_frames_and_allowed_free_translation(self):
        body = Body(new_id(), "Slider")
        p = Project(new_id(), bodies=(body,))
        j = joint_at(p, "glissiere", None, body.id, axis=(1.0, 0.0, 0.0))
        p = p.replace_object(j)
        self.assertFalse(p.diagnostics())
        moved = p.replace_object(replace(body, position=(0.3, 0.0, 0.0)))
        self.assertFalse(moved.diagnostics())
        self.assertTrue(
            moved.replace_object(replace(body, position=(0.0, 0.3, 0.0))).diagnostics()
        )
        commanded = moved.replace_object(replace(j, motion=Law("lineaire", (0.0, 0.1))))
        self.assertTrue(commanded.diagnostics())

    def test_json_round_trip_and_unknown_types(self):
        p = pendulum()
        p = p.replace_object(
            Load(
                new_id(),
                "Force",
                p.bodies[0].id,
                force=(Law("lineaire", (0.0, 1.0)), Law(), Law()),
            )
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "project.json"
            save_project(path, p)
            self.assertEqual(load_project(path), p)
            data = json.loads(path.read_text())
            data["project"]["bodies"][0]["shape"] = "unavailable-type"
            path.write_text(json.dumps(data))
            with self.assertRaises(ValueError):
                load_project(path)
            self.assertIn("unavailable-type", path.read_text())

    def test_g0_import_preserves_original_inertia_and_parameters(self):
        parameters = Parameters(length=2.0, angle_deg=-45.0, mass=0.3)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "old.json"
            save_parameters(path, parameters)
            before = path.read_bytes()
            p = import_g0(path)
            self.assertEqual(path.read_bytes(), before)
            np.testing.assert_allclose(p.bodies[0].inertia(), np.diag([1e-8] * 3).flat)
            self.assertAlmostEqual(np.linalg.norm(p.bodies[0].position), 2.0)
            self.assertEqual(p.bodies[0].mass, 0.3)
            self.assertFalse(p.diagnostics())

    def test_budget_scales_with_bodies(self):
        with self.assertRaises(ValueError):
            Project(
                new_id(),
                bodies=tuple(Body(new_id(), str(i)) for i in range(10)),
                duration=2.0,
                step=0.0001,
            )


if __name__ == "__main__":
    unittest.main()
