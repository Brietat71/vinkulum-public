"""An original STEP input built independently of the OCCT stack under test."""

import importlib.util
import json
import runpy
import tempfile
import unittest
from pathlib import Path

import numpy as np

FIXTURE = Path(__file__).with_name("fixtures") / "independent_step"


class FixtureSources(unittest.TestCase):
    def test_original_sources_reproduce_the_retained_part_and_reference(self):
        render = runpy.run_path(FIXTURE / "generate.py")["render"]
        reference = runpy.run_path(FIXTURE / "reference.py")["reference"]
        self.assertEqual(render(), (FIXTURE / "l-bracket.step").read_text())
        self.assertEqual(
            reference(), json.loads((FIXTURE / "reference.json").read_text())
        )


@unittest.skipUnless(
    importlib.util.find_spec("build123d"), "Install the OCCT 8 CAD extra"
)
class IndependentStep(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from vinkulum_studio.cad import execute
        from vinkulum_studio.document import Body

        cls.expected = json.loads((FIXTURE / "reference.json").read_text())
        cls.body = Body.from_dict(
            execute(
                {
                    "operation": "import_step",
                    "path": str(FIXTURE / "l-bracket.step"),
                    "density": cls.expected["density_kg_m3"],
                }
            )["body"]
        )

    def compare_reference(self, body):
        expected = self.expected
        epsilon, length = (
            expected["relative_scale_tolerance"],
            expected["characteristic_length_m"],
        )
        np.testing.assert_allclose(
            body.cad.volume_m3,
            expected["volume_m3"],
            rtol=0,
            atol=epsilon * expected["volume_m3"],
        )
        np.testing.assert_allclose(
            body.mass, expected["mass_kg"], rtol=0, atol=epsilon * expected["mass_kg"]
        )
        np.testing.assert_allclose(
            body.position, expected["centre_m"], rtol=0, atol=epsilon * length
        )
        np.testing.assert_allclose(
            np.array(body.inertia()).reshape(3, 3),
            expected["inertia_kg_m2"],
            rtol=0,
            atol=epsilon * expected["mass_kg"] * length**2,
        )

    def test_imported_single_solid_matches_box_integrals_in_si(self):
        from vinkulum_studio.cad import read_brep

        solid = read_brep(self.body.cad.brep_mm)
        self.assertEqual(len(solid.solids()), 1)
        self.assertTrue(solid.is_valid)
        self.compare_reference(self.body)
        tensor = np.array(self.expected["inertia_kg_m2"])
        self.assertTrue(
            all(abs(tensor[i, j]) > 1e-7 for i, j in ((0, 1), (0, 2), (1, 2)))
        )

    def test_wrong_step_length_unit_fails_the_same_physical_comparison(self):
        from vinkulum_studio.cad import execute
        from vinkulum_studio.document import Body

        text = (FIXTURE / "l-bracket.step").read_text()
        self.assertEqual(text.count("SI_UNIT(.MILLI.,.METRE.)"), 1)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wrong-unit.step"
            path.write_text(
                text.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)")
            )
            wrong = Body.from_dict(
                execute(
                    {
                        "operation": "import_step",
                        "path": str(path),
                        "density": self.expected["density_kg_m3"],
                    }
                )["body"]
            )
        with self.assertRaises(AssertionError):
            self.compare_reference(wrong)
        # Real reader conversion: cm instead of mm multiplies m by 10^3 and I by 10^5.
        np.testing.assert_allclose(wrong.mass / self.body.mass, 1000, rtol=1e-8)
        np.testing.assert_allclose(
            wrong.inertia(), np.array(self.body.inertia()) * 100000, rtol=1e-8, atol=0
        )


if __name__ == "__main__":
    unittest.main()
