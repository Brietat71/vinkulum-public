"""Retained worker measurements remain readable without loading the CAD engine."""

import json
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


@unittest.skipUnless(sys.platform == "linux", "The process experiment targets Linux")
class CadServiceArchive(unittest.TestCase):
    def test_retained_geometry_reopens_and_changed_response_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with zipfile.ZipFile(
                ROOT / "docs/bancs/cad-service-2026/captured-results.zip"
            ) as archive:
                archive.extractall(root)
            command = [
                sys.executable,
                str(ROOT / "ci/qualify_cad_service.py"),
                str(root),
                "--verify",
            ]
            result = subprocess.run(
                command, cwd=root, capture_output=True, text=True, timeout=30
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn(
                "Verified 45 retained CAD results without importing OCCT", result.stdout
            )
            path = root / "resident/001-output.json"
            response = json.loads(path.read_text())
            response["body"]["mass"] *= 1.01
            path.write_text(json.dumps(response))
            changed = subprocess.run(
                command, cwd=root, capture_output=True, text=True, timeout=30
            )
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn("resident/001-output.json", changed.stderr)


if __name__ == "__main__":
    unittest.main()
