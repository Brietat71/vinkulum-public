"""Real external bridge worker against a retained official-FreeCAD STEP capture."""

import importlib.util
import json
import os
import runpy
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
ARCHIVE = ROOT / "docs/bancs/freecad-bridge-2026/record.zip"
WORKER = ROOT / "apps/freecad/worker.py"
CAD = importlib.util.find_spec("build123d") is not None


@unittest.skipUnless(CAD, "Install the OCCT 8 CAD extra")
class FreecadBridgeWorker(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        with zipfile.ZipFile(ARCHIVE) as archive:
            for name in ("request.json", "part.step"):
                (self.directory / name).write_bytes(
                    archive.read(f"record/length-800-step-0.005/{name}")
                )

    def launch(self):
        environment = os.environ.copy()
        for name in ("PYTHONHOME", "PYTHONPATH"):
            environment.pop(name, None)
        return subprocess.run(
            [sys.executable, str(WORKER), str(self.directory)],
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=60,
        )

    def test_captured_freecad_part_runs_and_matches_independent_pendulum(self):
        result = self.launch()
        self.assertEqual(result.returncode, 0, result.stdout)
        data = json.loads((self.directory / "result/playback.json").read_text())
        self.assertAlmostEqual(data["properties_si"]["mass_kg"], 3.744, places=10)
        self.assertEqual(data["manifest"]["threads"], 2)
        self.assertEqual(data["manifest"]["scientific_status"], "NotAssessed")
        times = np.array(data["time_s"])
        positions = np.array(data["position_m"])
        expected = runpy.run_path(ROOT / "apps/freecad/verify.py")["reference"](
            times, 0.8, 32
        )
        np.testing.assert_allclose(
            np.arctan2(positions[:, 0], -positions[:, 2]),
            expected,
            atol=1.3e-4,
            rtol=0,
        )

    def test_disagreement_with_captured_mass_stops_before_simulation(self):
        path = self.directory / "request.json"
        request = json.loads(path.read_text())
        request["properties_si"]["mass_kg"] *= 2
        path.write_text(json.dumps(request))
        result = self.launch()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Not equal to tolerance", result.stdout)
        self.assertFalse((self.directory / "project.vinkulum.json").exists())
        self.assertFalse((self.directory / "result/trajectory.npz").exists())
        self.assertFalse((self.directory / "result/playback.json").exists())


if __name__ == "__main__":
    unittest.main()
