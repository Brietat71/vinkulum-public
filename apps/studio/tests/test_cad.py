"""Independent SI mass-property references and real OCCT 8 / desktop recipes."""

import importlib.util
import json
import math
import os
import tempfile
import time
import unittest
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
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

CAD = importlib.util.find_spec("build123d") is not None


@unittest.skipUnless(CAD, "Install the OCCT 8 CAD extra")
class CadKernelRecipe(unittest.TestCase):
    def create(self, kind="box", dimensions=(20, 30, 40), **options):
        from vinkulum_studio.cad import execute

        return Body.from_dict(
            execute(
                {
                    "operation": kind,
                    "dimensions_mm": dimensions,
                    "density": 7800,
                    **options,
                }
            )["body"]
        )

    def test_actual_occt_version_and_analytic_box_sphere_cylinder(self):
        import OCP

        self.assertGreaterEqual(int(OCP.__version__.split(".")[0]), 8)
        cases = (
            (
                "box",
                (20, 30, 40),
                0.02 * 0.03 * 0.04,
                np.diag(
                    [
                        (0.03**2 + 0.04**2) / 12,
                        (0.02**2 + 0.04**2) / 12,
                        (0.02**2 + 0.03**2) / 12,
                    ]
                ),
            ),
            ("sphere", (20,), 4 / 3 * math.pi * 0.02**3, np.eye(3) * (2 / 5 * 0.02**2)),
            (
                "cylinder",
                (20, 40),
                math.pi * 0.02**2 * 0.04,
                np.diag([(3 * 0.02**2 + 0.04**2) / 12] * 2 + [0.02**2 / 2]),
            ),
        )
        for kind, dims, volume, unit_inertia in cases:
            with self.subTest(kind=kind):
                body = self.create(kind, dims)
                self.assertAlmostEqual(body.cad.volume_m3 / volume, 1, places=10)
                self.assertAlmostEqual(body.mass / (7800 * volume), 1, places=10)
                np.testing.assert_allclose(
                    np.array(body.inertia()).reshape(3, 3),
                    unit_inertia * body.mass,
                    rtol=1e-9,
                    atol=1e-15,
                )
                np.testing.assert_allclose(body.position, (0, 0, 0), atol=1e-12)
                self.assertGreater(len(body.cad.triangles), 3)

    def test_extrusion_center_and_mass_scaling(self):
        body = self.create(
            "extrude_rectangle", (20, 30, 40), position_mm=(100, 200, 300)
        )
        np.testing.assert_allclose(body.position, (0.1, 0.2, 0.32), atol=1e-12)
        np.testing.assert_allclose(
            replace(body, mass=body.mass * 3).inertia(),
            np.array(body.inertia()) * 3,
            rtol=1e-14,
        )
        np.testing.assert_allclose(
            np.min(body.cad.vertices_m, axis=0), (-0.01, -0.015, -0.02), atol=1e-10
        )

    def test_boolean_volume_and_failed_operation(self):
        from vinkulum_studio.cad import execute

        a = self.create("box", (100, 60, 20))
        b = self.create("cylinder", (10, 40))
        request = {"operation": "cut", "a": asdict(a), "b": asdict(b), "density": 7800}
        result = Body.from_dict(execute(request)["body"])
        expected = 0.1 * 0.06 * 0.02 - math.pi * 0.01**2 * 0.02
        self.assertAlmostEqual(result.cad.volume_m3 / expected, 1, places=10)
        self.assertEqual(a.id, result.id)
        self.assertEqual(
            a.cad.operations, self.create("box", (100, 60, 20)).cad.operations
        )
        with self.assertRaises(ValueError):
            execute({**request, "b": asdict(a)})
        # A union of disconnected pieces cannot silently become one mechanical body.
        with self.assertRaises(ValueError):
            execute(
                {
                    **request,
                    "operation": "fuse",
                    "b": asdict(replace(b, position=(10, 0, 0))),
                }
            )

    def test_step_round_trip_and_unit_conversion(self):
        import build123d as bd
        from vinkulum_studio.cad import execute

        a = self.create(position_mm=(100, 200, 300))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "part.step"
            path.write_text(
                execute({"operation": "export_step", "a": asdict(a)})["step"]
            )
            b = Body.from_dict(
                execute(
                    {"operation": "import_step", "path": str(path), "density": 7800}
                )["body"]
            )
            np.testing.assert_allclose(a.position, b.position, atol=1e-11)
            np.testing.assert_allclose(a.inertia(), b.inertia(), rtol=1e-9, atol=1e-14)
            self.assertAlmostEqual(a.mass / b.mass, 1, places=10)
            # STEP in metres must be converted by the reader to its internal mm.
            bd.export_step(bd.Box(0.02, 0.03, 0.04), path, unit=bd.Unit.M)
            c = Body.from_dict(
                execute(
                    {"operation": "import_step", "path": str(path), "density": 7800}
                )["body"]
            )
            np.testing.assert_allclose(c.dimensions, (0.02, 0.03, 0.04), atol=1e-10)

    def test_fillet_and_rotated_inertia(self):
        import build123d as bd
        from vinkulum_studio.cad import capture_shape, execute
        from vinkulum_studio.editor import euler_matrix

        a = self.create()
        rounded = Body.from_dict(
            execute(
                {"operation": "fillet", "a": asdict(a), "radius_mm": 1, "density": 7800}
            )["body"]
        )
        self.assertLess(rounded.mass, a.mass)
        R = np.array(euler_matrix((25, 30, 40))).reshape(3, 3)
        rotated = capture_shape(
            bd.Box(20, 30, 40).moved(bd.Location((100, 200, 300), (25, 30, 40))),
            name="Rotated",
            density=7800,
        )
        # build123d uses intrinsic XYZ here; derive the rotation from its matrix.
        transform = bd.Location((100, 200, 300), (25, 30, 40)).wrapped.Transformation()
        R = np.array(
            [[transform.Value(i + 1, j + 1) for j in range(3)] for i in range(3)]
        )
        np.testing.assert_allclose(
            np.array(rotated.inertia()).reshape(3, 3),
            R @ np.array(a.inertia()).reshape(3, 3) @ R.T,
            atol=1e-13,
        )
        np.testing.assert_allclose(rotated.position, (0.1, 0.2, 0.3), atol=1e-12)

    def test_rebase_attachments_preserves_world_and_history(self):
        from vinkulum_studio.cad import execute

        a = self.create("box", (100, 60, 20))
        tool = self.create("box", (60, 100, 100), position_mm=(50, 0, 0))
        p = Project(new_id(), "CAD", bodies=(a, tool), gravity=(0, 0, 0))
        joint = joint_at(p, "pivot", None, a.id, point=(-0.05, 0, 0))
        load = Load(
            new_id(),
            "Load",
            a.id,
            point=(0.01, 0.02, 0),
            force=(Law("constante", (1.0,)), Law(), Law()),
        )
        p = p.replace_object(joint).replace_object(load)
        b = Body.from_dict(
            execute(
                {"operation": "cut", "a": asdict(a), "b": asdict(tool), "density": 7800}
            )["body"]
        )
        updated = replace_cad_body(p, b)
        self.assertGreater(np.linalg.norm(np.array(b.position) - a.position), 0.001)
        np.testing.assert_allclose(
            np.array(a.position) + joint.pb,
            np.array(b.position) + updated.joints[0].pb,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            np.array(a.position) + load.point,
            np.array(b.position) + updated.loads[0].point,
            atol=1e-12,
        )
        self.assertEqual(updated.diagnostics(), ())
        history = History(p)
        history.commit(updated)
        history.undo()
        self.assertEqual(history.current.bodies, p.bodies)

    def test_project_brep_persistence_and_invalid_mesh(self):
        body = self.create()
        p = Project(new_id(), "CAD", bodies=(body,))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cad.json"
            save_project(path, p)
            self.assertEqual(json.loads(path.read_text())["schema_version"], 3)
            self.assertEqual(load_project(path), p)
            invalid = asdict(body)
            invalid["cad"]["triangles"] = [[0, 1, 999999]] * 4
            with self.assertRaises(ValueError):
                Body.from_dict(invalid)


@unittest.skipUnless(
    CAD and os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires CAD and desktop"
)
class CadDesktopRecipe(unittest.TestCase):
    def test_real_cad_process_then_mechanical_worker_and_render(self):
        from PySide6.QtTest import QTest
        from PySide6.QtWidgets import QApplication
        from vinkulum_studio.bundle_check import check_scene_image
        from vinkulum_studio.cad_dialog import CadDialog
        from vinkulum_studio.editor import EditorWindow

        self.app = QApplication.instance() or QApplication([])
        window = EditorWindow()
        window._discard_allowed = lambda: True
        self.addCleanup(window.close)
        window.show()
        window.new_project()
        dialog = CadDialog(window.project, None, window)
        self.addCleanup(dialog.reject)
        dialog.show()
        dialog.start()
        deadline = time.monotonic() + 30
        while dialog.process is not None and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertIsNotNone(dialog.result_data, dialog.status.text())
        body = Body.from_dict(dialog.result_data["body"])
        window._commit(replace_cad_body(window.project, body), fit=True)
        window.duration.setText("0.05")
        window.step.setText("0.005")
        window.run()
        deadline = time.monotonic() + 30
        while window.controller.process is not None and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertIsNotNone(window.result, window.status.text())
        self.assertEqual(window.result.project.bodies[0].cad, body.cad)
        with tempfile.TemporaryDirectory() as directory:
            image = Path(directory) / "cad.png"
            window.viewport.screenshot(image)
            self.assertEqual(check_scene_image(image)["nonblack_probes"], 100)

    def test_failed_cad_process_keeps_document(self):
        from PySide6.QtTest import QTest
        from PySide6.QtWidgets import QApplication
        from vinkulum_studio.cad_dialog import CadDialog

        self.app = QApplication.instance() or QApplication([])
        p = Project(new_id(), "Empty")
        dialog = CadDialog(p, None)
        self.addCleanup(dialog.reject)
        dialog.operation.setCurrentIndex(dialog.operation.findData("fillet"))
        dialog.show()
        dialog.start()
        deadline = time.monotonic() + 30
        while dialog.process is not None and time.monotonic() < deadline:
            QTest.qWait(10)
        self.assertIsNone(dialog.result_data)
        self.assertEqual(dialog.project, p)
        self.assertIn("Select", dialog.status.text())


if __name__ == "__main__":
    unittest.main()
