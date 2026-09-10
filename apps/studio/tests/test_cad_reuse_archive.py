"""Retained GUI geometry is checked independently and rejects changed captures."""

import json
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


class CadReuseArchive(unittest.TestCase):
    def test_retained_geometry_reopens_and_changed_mass_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with zipfile.ZipFile(
                ROOT / "docs/bancs/cad-reuse-dev6/captured-gui.zip"
            ) as archive:
                archive.extractall(root)
            command = [sys.executable, str(ROOT / "ci/verify_cad_reuse.py"), str(root)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("Verified 86 captured GUI CAD results", result.stdout)
            path = root / "cached-100/0-box.json"
            body = json.loads(path.read_text())
            body["mass"] *= 1.01
            path.write_text(json.dumps(body))
            changed = subprocess.run(
                command, capture_output=True, text=True, timeout=30
            )
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn("cached-100/0-box.json", changed.stderr)


if __name__ == "__main__":
    unittest.main()
