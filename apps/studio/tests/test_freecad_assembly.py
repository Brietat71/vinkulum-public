"""Actual two-body engine runs from a retained native FreeCAD Assembly capture."""

import importlib.util
import json
import math
import os
import runpy
import subprocess
import sys
import tempfile
import unittest
import uuid
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
ARCHIVE = ROOT / "docs/bancs/freecad-assembly-2026/record.zip"
WORKER = ROOT / "apps/freecad/assembly_worker.py"
CAD = importlib.util.find_spec("build123d") is not None


@unittest.skipUnless(CAD, "Install the OCCT 8 CAD extra")
class FreecadAssemblyWorker(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        # Native parent placement is translated (125,-75,210) mm and yawed 17°.
        with zipfile.ZipFile(ARCHIVE) as archive:
            raw = archive.read("record/trial-2/request.json")
            self.request = json.loads(raw)
            (self.directory / "request.json").write_bytes(raw)
            for body in self.request["bodies"]:
                name = body["step_file"]
                (self.directory / name).write_bytes(
                    archive.read("record/trial-2/" + name)
                )

    def launch(self):
        environment = os.environ.copy()
        for name in ("PYTHONHOME", "PYTHONPATH"):
            environment.pop(name, None)
        return subprocess.run(
            [sys.executable, str(WORKER), str(self.directory)],
            cwd=self.directory,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=60,
        )

    def test_transported_native_assembly_matches_independent_double_pendulum(self):
        process = self.launch()
        self.assertEqual(process.returncode, 0, process.stdout)
        result = json.loads((self.directory / "result/playback.json").read_text())
        self.assertEqual(result["manifest"]["threads"], 2)
        self.assertEqual(result["manifest"]["scientific_status"], "NotAssessed")
        names = [body["name"] for body in self.request["bodies"]]
        self.assertEqual(
            result["body_ids"], [b["body_id"] for b in self.request["bodies"]]
        )
        for name, mass in (("Upper", 1.296), ("Lower", 0.972)):
            self.assertAlmostEqual(
                result["properties_si"][names.index(name)]["mass_kg"], mass, places=10
            )
        yaw = math.radians(17)
        rotation = np.array(
            (
                (math.cos(yaw), -math.sin(yaw), 0),
                (math.sin(yaw), math.cos(yaw), 0),
                (0, 0, 1),
            )
        )
        positions = (np.array(result["position_m"]) - (0.125, -0.075, 0.210)) @ rotation
        upper = positions[:, names.index("Upper")]
        lower = positions[:, names.index("Lower")] - 2 * upper
        actual = np.column_stack(
            [np.arctan2(p[:, 0], -p[:, 2]) for p in (upper, lower)]
        )
        expected = runpy.run_path(ROOT / "apps/freecad/verify_assembly.py")[
            "reference"
        ](np.array(result["time_s"]), 32)
        np.testing.assert_allclose(actual, expected, atol=2e-3, rtol=0)

    def assert_rejected_before_simulation(self, diagnostic):
        (self.directory / "request.json").write_text(json.dumps(self.request))
        process = self.launch()
        self.assertNotEqual(process.returncode, 0)
        self.assertIn(diagnostic, process.stdout)
        for name in (
            "project.vinkulum.json",
            "result/trajectory.npz",
            "result/playback.json",
        ):
            self.assertFalse((self.directory / name).exists(), name)

    def test_one_inconsistent_body_mass_rejects_the_whole_mechanism(self):
        self.request["bodies"][1]["properties_si"]["mass_kg"] *= 2
        self.assert_rejected_before_simulation("Not equal to tolerance")

    def test_joint_to_uncaptured_component_rejects_the_mechanism(self):
        self.request["joints"][0]["body_b"] = str(uuid.uuid4())
        self.assert_rejected_before_simulation("A joint references an uncaptured body")


if __name__ == "__main__":
    unittest.main()
