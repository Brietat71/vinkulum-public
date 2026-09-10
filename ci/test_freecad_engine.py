"""Installer refusal must preserve existing environments and avoid partial writes."""

import os
import platform
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("freecad_engine.py")


class InstallationRefusal(unittest.TestCase):
    def test_existing_environment_is_unchanged(self):
        with tempfile.TemporaryDirectory() as folder:
            destination = Path(folder) / "working-engine"
            destination.mkdir()
            marker = destination / "pyvenv.cfg"
            marker.write_bytes(b"existing engine configuration\n")
            before = marker.stat().st_mtime_ns
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(destination)],
                capture_output=True,
                check=False,
                text=True,
                timeout=10,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Choose a new directory", result.stderr)
            self.assertEqual(list(destination.iterdir()), [marker])
            self.assertEqual(marker.read_bytes(), b"existing engine configuration\n")
            self.assertEqual(marker.stat().st_mtime_ns, before)

    @unittest.skipUnless(
        sys.platform == "linux" and platform.machine() == "x86_64", "Linux recipe"
    )
    def test_missing_prerequisite_does_not_create_destination(self):
        with tempfile.TemporaryDirectory() as folder:
            destination = Path(folder) / "new-engine"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), str(destination)],
                env={**os.environ, "PATH": ""},
                check=False,
                capture_output=True,
                text=True,
                timeout=10,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Install uv", result.stderr)
            self.assertFalse(destination.exists())


if __name__ == "__main__":
    unittest.main()
