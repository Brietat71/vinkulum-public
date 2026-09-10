"""Parametric regeneration, design frames and persistent solid-feature identities."""

import importlib.util
import json
import math
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
from vinkulum_studio.cad_history import CadFeature, CadRecipe
from vinkulum_studio.document import (
    Body,
    History,
    Law,
    Load,
    Project,
    joint_at,
    load_project,
    new_id,
    replace_cad_body,
    save_project,
)


class FeatureGraph(unittest.TestCase):
    def test_identity_ordering_and_invalid_dependencies(self):
        base = CadFeature(new_id(), "Stock", "box", dimensions_mm=(100, 60, 20))
        tool = CadFeature(new_id(), "Tool", "cylinder", dimensions_mm=(10, 40))
        cut = CadFeature(new_id(), "Hole", "cut", inputs=(base.id, tool.id))
        recipe = CadRecipe((cut, tool, base), cut.id)
        self.assertEqual(recipe.ordered(), (base, tool, cut))
        self.assertEqual(
            recipe.fingerprint(),
            replace(recipe, features=(base, cut, tool)).fingerprint(),
        )
        renamed = recipe.replace_feature(replace(base, name="Renamed stock"))
        self.assertEqual(renamed.ordered()[-1].inputs, cut.inputs)
        self.assertEqual(CadRecipe.from_dict(asdict(renamed)), renamed)
        for features in (
            (base, tool),
            (base, base, cut),
            (base, tool, replace(cut, inputs=(base.id, new_id()))),
            (
                base,
                tool,
                cut,
                CadFeature(new_id(), "Unused", "sphere", dimensions_mm=(1,)),
            ),
            (base, tool, replace(cut, inputs=(cut.id, tool.id))),
        ):
            with self.subTest(features=features), self.assertRaises(ValueError):
                CadRecipe(features, cut.id)
        with self.assertRaises(ValueError):
            replace(base, dimensions_mm=(True, 10, 10))
        with self.assertRaises(ValueError):
            replace(base, orientation=(1, 0, 0, 0, 1, 0, 0, 0, -1))
        with self.assertRaises(ValueError):
            replace(cut, position_mm=(1, 0, 0))


@unittest.skipUnless(
    importlib.util.find_spec("build123d"), "Install the OCCT 8 CAD environment"
)
class ParametricSolids(unittest.TestCase):
    def execute(self, **request):
        from vinkulum_studio.cad import execute

        return Body.from_dict(execute(request)["body"])

    def box(self, dimensions=(100, 60, 20), **kw):
        return self.execute(
            operation="box", dimensions_mm=dimensions, density=7800, **kw
        )

    def regenerate(self, body, recipe, density=7800):
        return self.execute(
            operation="regenerate",
            a=asdict(body),
            recipe=asdict(recipe),
            density=density,
        )

    def test_upstream_edit_propagates_to_cut_and_fillet_with_stable_ids(self):
        a = self.box()
        b = self.execute(
            operation="cylinder", dimensions_mm=(10, 40), position_mm=(20, 0, 0)
        )
        cut = self.execute(operation="cut", a=asdict(a), b=asdict(b), density=7800)
        recipe = cut.cad.recipe
        base = recipe.features[0]
        changed = recipe.replace_feature(
            replace(base, name="Wider stock", dimensions_mm=(120, 60, 20))
        )
        result = self.regenerate(cut, changed)
        expected = 0.12 * 0.06 * 0.02 - math.pi * 0.01**2 * 0.02
        self.assertAlmostEqual(result.cad.volume_m3 / expected, 1, places=10)
        self.assertEqual(result.id, a.id)
        self.assertEqual(
            tuple(f.id for f in result.cad.recipe.features),
            tuple(f.id for f in recipe.features),
        )
        self.assertEqual(result.cad.recipe.features[-2].source_body_id, b.id)
        # Newer edits of tool body B do not alter the captured cutter.
        moved = replace(b, position=(1, 0, 0))
        self.assertNotEqual(asdict(moved), asdict(b))
        repeated = self.regenerate(result, result.cad.recipe)
        self.assertAlmostEqual(repeated.mass / result.mass, 1, places=11)
        rounded = self.execute(
            operation="fillet", a=asdict(result), radius_mm=1, density=7800
        )
        narrower = rounded.cad.recipe.replace_feature(
            replace(rounded.cad.recipe.features[0], dimensions_mm=(110, 60, 20))
        )
        updated = self.regenerate(rounded, narrower)
        self.assertLess(updated.mass, rounded.mass)
        self.assertEqual(updated.cad.recipe.root, rounded.cad.recipe.root)
        reordered = self.regenerate(
            updated,
            replace(
                updated.cad.recipe,
                features=tuple(reversed(updated.cad.recipe.features)),
            ),
        )
        self.assertAlmostEqual(reordered.mass / updated.mass, 1, places=11)

    def test_extrusion_design_plane_and_exact_inertia_after_height_edit(self):
        from vinkulum_studio.editor import euler_matrix

        a = self.execute(
            operation="extrude_rectangle",
            dimensions_mm=(20, 30, 40),
            position_mm=(100, 200, 300),
            density=7800,
        )
        R = np.array(euler_matrix((24, 31, -17))).reshape(3, 3)
        a = replace(a, orientation=tuple(R.flat), position=(1000, -500, 700))
        origin = np.array(a.position) + R @ a.cad.recipe.origin_in_body_m
        base = a.cad.recipe.features[0]
        result = self.regenerate(
            a, a.cad.recipe.replace_feature(replace(base, dimensions_mm=(20, 30, 80)))
        )
        np.testing.assert_allclose(
            result.position, origin + R @ np.array([0, 0, 0.04]), rtol=0, atol=2e-13
        )
        np.testing.assert_allclose(
            np.array(result.position) + R @ result.cad.recipe.origin_in_body_m,
            origin,
            rtol=0,
            atol=2e-13,
        )
        self.assertAlmostEqual(result.mass, 7800 * 0.02 * 0.03 * 0.08, places=12)
        expected = (
            np.diag(
                [
                    (0.03**2 + 0.08**2) / 12,
                    (0.02**2 + 0.08**2) / 12,
                    (0.02**2 + 0.03**2) / 12,
                ]
            )
            * result.mass
        )
        np.testing.assert_allclose(
            np.array(result.inertia()).reshape(3, 3), expected, rtol=1e-10, atol=1e-15
        )
        bad_origin = replace(result.cad.recipe, origin_in_body_m=(0, 0, 0))
        with self.assertRaisesRegex(ValueError, "design origin"):
            self.regenerate(result, bad_origin)

    def test_tool_placement_after_asymmetric_cut_and_rotated_body_frame(self):
        from vinkulum_studio.editor import euler_matrix

        a = self.box()
        first = self.execute(
            operation="cylinder", dimensions_mm=(5, 40), position_mm=(25, 0, 0)
        )
        a = self.execute(operation="cut", a=asdict(a), b=asdict(first))
        R = np.array(euler_matrix((20, 30, 40))).reshape(3, 3)
        a = replace(a, orientation=tuple(R.flat), position=(10, -12, 2))
        design_origin = np.array(a.position) + R @ a.cad.recipe.origin_in_body_m
        second = self.execute(operation="cylinder", dimensions_mm=(4, 40))
        second = replace(
            second,
            orientation=tuple(R.flat),
            position=tuple(design_origin + R @ np.array([-0.025, 0, 0])),
        )
        a = self.execute(operation="cut", a=asdict(a), b=asdict(second))
        expected = 0.1 * 0.06 * 0.02 - math.pi * (0.005**2 + 0.004**2) * 0.02
        self.assertAlmostEqual(a.cad.volume_m3 / expected, 1, places=9)
        np.testing.assert_allclose(
            np.array(a.position) + R @ a.cad.recipe.origin_in_body_m,
            design_origin,
            rtol=0,
            atol=2e-13,
        )
        self.assertAlmostEqual(
            self.regenerate(a, a.cad.recipe).mass / a.mass, 1, places=10
        )

    def test_failed_feature_and_attachment_preservation_with_undo(self):
        a = self.box()
        p = Project(new_id(), "Mounted stock", bodies=(a,))
        joint = joint_at(p, "encastrement", None, a.id, point=(0, 0, 0))
        load = Load(
            new_id(),
            "Force",
            a.id,
            point=(0.01, 0.02, 0),
            force=(Law(values=(1.0,)), Law(), Law()),
        )
        p = p.replace_object(joint).replace_object(load)
        base = a.cad.recipe.features[0]
        result = self.regenerate(
            a, a.cad.recipe.replace_feature(replace(base, position_mm=(10, 0, 0)))
        )
        changed = replace_cad_body(p, result)
        np.testing.assert_allclose(
            np.array(a.position) + p.joints[0].pb,
            np.array(result.position) + changed.joints[0].pb,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            np.array(a.position) + p.loads[0].point,
            np.array(result.position) + changed.loads[0].point,
            atol=1e-12,
        )
        self.assertEqual(changed.diagnostics(), ())
        history = History(p)
        history.commit(changed)
        history.undo()
        self.assertEqual(history.current.bodies, p.bodies)
        bad = CadFeature(
            new_id(),
            "Impossible rounding",
            "fillet",
            inputs=(base.id,),
            dimensions_mm=(1000,),
        )
        with self.assertRaisesRegex(ValueError, "Impossible rounding"):
            self.regenerate(a, replace(a.cad.recipe, features=(base, bad), root=bad.id))
        self.assertEqual(p.bodies, (a,))

    def test_schema_round_trip_legacy_snapshot_and_downgrade_rejection(self):
        a = self.box()
        p = Project(new_id(), "History", bodies=(a,))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "project.json"
            save_project(path, p)
            self.assertEqual(json.loads(path.read_text())["schema_version"], 3)
            self.assertEqual(load_project(path), p)
            bad = json.loads(path.read_text())
            bad["schema_version"] = 2
            path.write_text(json.dumps(bad))
            with self.assertRaisesRegex(ValueError, "schema 3"):
                load_project(path)
            legacy = replace(a, cad=replace(a.cad, recipe=None))
            save_project(path, replace(p, bodies=(legacy,)))
            self.assertEqual(json.loads(path.read_text())["schema_version"], 2)
            self.assertNotIn(
                "recipe", json.loads(path.read_text())["project"]["bodies"][0]["cad"]
            )
            self.assertEqual(load_project(path).bodies[0], legacy)
            rounded = self.execute(operation="fillet", a=asdict(legacy), radius_mm=1)
            self.assertEqual(rounded.cad.recipe.features[0].kind, "snapshot")
            self.assertAlmostEqual(
                self.regenerate(rounded, rounded.cad.recipe).mass / rounded.mass,
                1,
                places=10,
            )
