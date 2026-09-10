"""Curved metre STEP and a closed offset void, independent of the CAD exporter."""

import importlib.util
import json
import os
import re
import runpy
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

FIXTURE = Path(__file__).with_name("fixtures") / "analytic_step_corpus"
REFERENCE = json.loads((FIXTURE / "reference.json").read_text())
CAD = importlib.util.find_spec("build123d") is not None


def compare(body, expected):
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


class AnalyticSources(unittest.TestCase):
    def test_standard_library_sources_reproduce_both_parts_and_references(self):
        generator = runpy.run_path(FIXTURE / "generate.py")
        for name, function in (
            ("eccentric-bore-metres", "eccentric_bore"),
            ("offset-closed-void-mm", "offset_void"),
        ):
            self.assertEqual(
                generator[function](), (FIXTURE / f"{name}.step").read_text()
            )
        self.assertEqual(
            runpy.run_path(FIXTURE / "reference.py")["references"](), REFERENCE
        )


@unittest.skipUnless(CAD, "Install the OCCT 8 CAD extra")
class AnalyticStepCorpus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from vinkulum_studio.cad import execute
        from vinkulum_studio.document import Body

        cls.bodies = {
            name: Body.from_dict(
                execute(
                    {
                        "operation": "import_step",
                        "path": str(FIXTURE / f"{name}.step"),
                        "density": expected["density_kg_m3"],
                    }
                )["body"]
            )
            for name, expected in REFERENCE.items()
        }

    def test_eccentric_curved_bore_in_metres_matches_cylinder_integrals(self):
        from vinkulum_studio.cad import read_brep

        name = "eccentric-bore-metres"
        body = self.bodies[name]
        shape = read_brep(body.cad.brep_mm)
        self.assertTrue(shape.is_valid)
        self.assertEqual(len(shape.solids()), 1)
        self.assertEqual(len(shape.shells()), 1)
        self.assertEqual([f.geom_type.name for f in shape.faces()].count("CYLINDER"), 8)
        self.assertEqual([f.geom_type.name for f in shape.faces()].count("PLANE"), 2)
        compare(body, REFERENCE[name])
        tensor = np.array(REFERENCE[name]["inertia_kg_m2"])
        self.assertTrue(
            all(abs(tensor[i, j]) > 1e-7 for i, j in ((0, 1), (0, 2), (1, 2)))
        )

    def test_closed_offset_void_matches_signed_box_integrals(self):
        from vinkulum_studio.cad import read_brep

        name = "offset-closed-void-mm"
        body = self.bodies[name]
        shape = read_brep(body.cad.brep_mm)
        self.assertTrue(shape.is_valid)
        self.assertEqual(len(shape.solids()), 1)
        self.assertEqual(len(shape.shells()), 2)
        self.assertEqual(len(shape.faces()), 12)
        compare(body, REFERENCE[name])
        tensor = np.array(REFERENCE[name]["inertia_kg_m2"])
        self.assertTrue(
            all(abs(tensor[i, j]) > 1e-7 for i, j in ((0, 1), (0, 2), (1, 2)))
        )

    def import_variant(self, text, density):
        from vinkulum_studio.cad import execute
        from vinkulum_studio.document import Body

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "variant.step"
            path.write_text(text)
            return Body.from_dict(
                execute(
                    {"operation": "import_step", "path": str(path), "density": density}
                )["body"]
            )

    def test_metre_to_millimetre_misdeclaration_changes_real_reader_scale(self):
        name = "eccentric-bore-metres"
        text = (FIXTURE / f"{name}.step").read_text()
        self.assertEqual(text.count("SI_UNIT($,.METRE.)"), 1)
        wrong = self.import_variant(
            text.replace("SI_UNIT($,.METRE.)", "SI_UNIT(.MILLI.,.METRE.)"),
            REFERENCE[name]["density_kg_m3"],
        )
        with self.assertRaises(AssertionError):
            compare(wrong, REFERENCE[name])
        np.testing.assert_allclose(wrong.mass / self.bodies[name].mass, 1e-9, rtol=1e-7)
        np.testing.assert_allclose(
            wrong.inertia(),
            np.array(self.bodies[name].inertia()) * 1e-15,
            rtol=1e-7,
            atol=0,
        )

    def test_omitting_the_void_is_detected_by_the_same_physical_reference(self):
        from vinkulum_studio.cad import read_brep

        name = "offset-closed-void-mm"
        text, count = re.subn(
            r"BREP_WITH_VOIDS\('Offset closed cavity',(#\d+),\(#\d+\)\)",
            r"MANIFOLD_SOLID_BREP('Filled box',\1)",
            (FIXTURE / f"{name}.step").read_text(),
        )
        self.assertEqual(count, 1)
        filled = self.import_variant(text, REFERENCE[name]["density_kg_m3"])
        self.assertEqual(len(read_brep(filled.cad.brep_mm).shells()), 1)
        np.testing.assert_allclose(
            filled.cad.volume_m3, 60 * 44 * 30 * 1e-9, rtol=1e-12
        )
        with self.assertRaises(AssertionError):
            compare(filled, REFERENCE[name])

    def test_two_disjoint_parts_are_rejected_instead_of_becoming_one_body(self):
        text = runpy.run_path(FIXTURE / "generate.py")["offset_void"](extra_solid=True)
        with self.assertRaisesRegex(ValueError, "one valid, closed, connected solid"):
            self.import_variant(
                text, REFERENCE["offset-closed-void-mm"]["density_kg_m3"]
            )


@unittest.skipUnless(
    CAD and os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires desktop and OCCT 8"
)
class AnalyticWorkerImport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from PySide6.QtWidgets import QApplication

        cls.app = QApplication.instance() or QApplication([])

    def test_real_import_worker_preserves_previous_part_after_multipart_failure(self):
        from PySide6.QtTest import QTest
        from vinkulum_studio.cad_controller import CadController
        from vinkulum_studio.document import Project, new_id

        project = Project(new_id(), "Analytic imports")
        controller = CadController(project=project)
        self.addCleanup(controller.shutdown)
        errors = []
        controller.failed.connect(lambda kind, message: errors.append((kind, message)))

        def wait():
            deadline = time.monotonic() + 30
            while controller.process is not None and time.monotonic() < deadline:
                QTest.qWait(5)
            self.assertIsNone(controller.process, controller.last_log)

        for name, expected in REFERENCE.items():
            controller.start(
                {
                    "operation": "import_step",
                    "path": str(FIXTURE / f"{name}.step"),
                    "density": expected["density_kg_m3"],
                }
            )
            wait()
            self.assertFalse(errors, errors)
            admission = controller.last_admission
            compare(admission.body, expected)
            self.assertIs(admission.source_project, project)
            self.assertEqual(admission.project.bodies, (admission.body,))
        previous = controller.last_result
        accepted = controller.last_admission
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "two-parts.step"
            path.write_text(
                runpy.run_path(FIXTURE / "generate.py")["offset_void"](extra_solid=True)
            )
            controller.start(
                {"operation": "import_step", "path": str(path), "density": 2700}
            )
            wait()
        self.assertEqual(errors[-1][0], "operation_failed")
        self.assertIn("one valid, closed, connected solid", errors[-1][1])
        self.assertIs(controller.last_result, previous)
        self.assertIs(controller.last_admission, accepted)
        self.assertIsNone(controller.temporary)


if __name__ == "__main__":
    unittest.main()
