"""OCCT extrusion against two analytic cuboids, then persistent upstream edits."""

import importlib.util
import json
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from vinkulum_studio.cad_history import CadFeature, CadRecipe, SketchExtrusion
from vinkulum_studio.document import Body, Project, load_project, new_id, save_project
from vinkulum_studio.sketch import bracket_profile, profile, settled


@unittest.skipUnless(
    importlib.util.find_spec("build123d"), "Install OCCT 8 CAD dependencies"
)
class SketchCAD(unittest.TestCase):
    def execute(self, **request):
        from vinkulum_studio.cad import execute

        return Body.from_dict(execute(request)["body"])

    def test_bracket_mass_com_and_full_inertia_match_two_cuboids(self):
        sketch = bracket_profile()
        body = self.execute(
            operation="extrude_sketch",
            profile=asdict(sketch),
            dimensions_mm=(10,),
            density=7800,
        )
        pieces = (
            (np.array((0.08, 0.02, 0.01)), np.array((0.04, 0.01, 0.005))),
            (np.array((0.03, 0.04, 0.01)), np.array((0.015, 0.04, 0.005))),
        )
        masses = [7800 * np.prod(size) for size, _ in pieces]
        mass = sum(masses)
        centre = sum(m * point for m, (_, point) in zip(masses, pieces)) / mass
        inertia = np.zeros((3, 3))
        for m, (size, point) in zip(masses, pieces):
            d = point - centre
            inertia += m * (
                np.diag((sum(size**2) - size**2) / 12)
                + (d @ d) * np.eye(3)
                - np.outer(d, d)
            )
        self.assertAlmostEqual(body.mass / mass, 1, places=12)
        np.testing.assert_allclose(body.position, centre, rtol=0, atol=1e-13)
        np.testing.assert_allclose(
            np.array(body.inertia()).reshape(3, 3), inertia, rtol=1e-11, atol=1e-15
        )
        self.assertIsInstance(body.cad.recipe.features[0], SketchExtrusion)
        self.assertEqual(body.cad.recipe.features[0].profile, sketch)
        self.assertGreater(len(body.cad.triangles), 4)

    def test_upstream_dimension_preserves_graph_design_frame_and_schema_contract(self):
        from vinkulum_studio.editor import euler_matrix

        body = self.execute(
            operation="extrude_sketch",
            profile=asdict(bracket_profile()),
            dimensions_mm=(10,),
            position_mm=(50, -70, 90),
        )
        feature = body.cad.recipe.features[0]
        changed = settled(
            replace(
                feature.profile,
                constraints=tuple(
                    replace(c, values_mm=("100",)) if c.name == "Width" else c
                    for c in feature.profile.constraints
                ),
            )
        )
        recipe = body.cad.recipe.replace_feature(replace(feature, profile=changed))
        R = np.array(euler_matrix((20, -40, 30))).reshape(3, 3)
        body = replace(body, orientation=tuple(R.flat), position=(3, -4, 5))
        origin = np.array(body.position) + R @ body.cad.recipe.origin_in_body_m
        result = self.execute(
            operation="regenerate", a=asdict(body), recipe=asdict(recipe), density=7800
        )
        np.testing.assert_allclose(
            np.array(result.position) + R @ result.cad.recipe.origin_in_body_m,
            origin,
            rtol=0,
            atol=1e-13,
        )
        self.assertAlmostEqual(
            result.mass, 7800 * (100 * 20 + 30 * 40) * 10 * 1e-9, places=12
        )
        self.assertEqual(result.id, body.id)
        self.assertEqual(result.cad.recipe.features[0].id, feature.id)
        self.assertEqual(
            {p.id for p in changed.points}, {p.id for p in feature.profile.points}
        )
        shuffled = replace(
            changed,
            points=tuple(reversed(changed.points)),
            edges=tuple(reversed(changed.edges)),
            constraints=tuple(reversed(changed.constraints)),
        )
        self.assertEqual(
            recipe.fingerprint(),
            recipe.replace_feature(replace(feature, profile=shuffled)).fingerprint(),
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sketch.vinkulum.json"
            project = Project(new_id(), bodies=(result,))
            save_project(path, project)
            self.assertEqual(json.loads(path.read_text())["schema_version"], 4)
            self.assertEqual(load_project(path), project)
            data = json.loads(path.read_text())
            data["schema_version"] = 3
            path.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError, "schema 4"):
                load_project(path)
        with self.assertRaisesRegex(ValueError, "persistent profile"):
            CadFeature(
                new_id(), "Incomplete extrusion", "extrude_sketch", dimensions_mm=(10,)
            )
        self.assertEqual(CadRecipe.from_dict(asdict(recipe)), recipe)

    def test_invalid_profile_and_conflicting_dimensions_are_refused_before_extrusion(
        self,
    ):
        sketch = bracket_profile()
        c = next(c for c in sketch.constraints if c.name == "Width")
        conflict = replace(c, id=new_id(), values_mm=("90",))
        bad = replace(sketch, constraints=sketch.constraints + (conflict,))
        for value in (bad, profile(((0, 0), (10, 10), (0, 10), (10, 0)))):
            with self.assertRaises(ValueError):
                self.execute(
                    operation="extrude_sketch",
                    profile=asdict(value),
                    dimensions_mm=(10,),
                )


if __name__ == "__main__":
    unittest.main()
