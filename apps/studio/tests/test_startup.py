"""The public launcher must reach an exposed, English, rendered window."""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from vinkulum_studio import __version__


@unittest.skipUnless(
    os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires native desktop OpenGL"
)
class NormalStartup(unittest.TestCase):
    def test_module_launcher_creates_the_normal_window(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "vinkulum_studio",
                    "--startup-check",
                    str(output),
                ],
                cwd=directory,
                env={**os.environ, "XDG_CONFIG_HOME": str(output / "config")},
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            report = json.loads((output / "startup-check.json").read_text())
            self.assertEqual(report["status"], "passed")
            self.assertEqual(report["locale"], "en_GB")
            self.assertEqual(report["app_version"], __version__)
            self.assertIn("Vinkulum Studio", report["title"])
            self.assertFalse(report["frozen"])
