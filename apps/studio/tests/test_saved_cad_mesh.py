"""The published CAD/CalculiX capture opens in a fresh engine-free process."""

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path


class SavedCADMesh(unittest.TestCase):
    def test_capture_hashes_and_fresh_reopening_without_gui_or_engines(self):
        root = Path(__file__).resolve().parents[3]
        archive = root / "docs/bancs/studio-cad-meshing-060/mesh-and-static-example.zip"
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory)
            with zipfile.ZipFile(archive) as stream:
                stream.extractall(destination)
            manifest = json.loads((destination / "files.sha256.json").read_text())
            self.assertEqual(
                set(manifest),
                {
                    str(p.relative_to(destination))
                    for p in destination.rglob("*")
                    if p.is_file() and p.name != "files.sha256.json"
                },
            )
            for name, digest in manifest.items():
                self.assertEqual(
                    hashlib.sha256((destination / name).read_bytes()).hexdigest(),
                    digest,
                )
            code = """import sys, subprocess
from pathlib import Path
def forbidden(*args, **kwargs):
    raise AssertionError('Reopening must not execute an external process')
subprocess.run = forbidden
from vinkulum_studio.meshing import load_mesh
from vinkulum_studio.calculix import load_static_result
root=Path(sys.argv[1])
before={p: (p.read_bytes(), p.stat().st_mtime_ns) for p in root.rglob('*') if p.is_file()}
solid,_=load_mesh(root/'mesh')
study,result,_=load_static_result(root/'calculation')
assert study.mesh_binding.solid == solid
assert result['schema'] == 3
assert result['input_transport'] == study.input_transport
assert result['strain_energy_J'] > 0
assert before == {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in before}
for name in ('PySide6','vtkmodules','OCP','build123d','gmsh'):
    assert name not in sys.modules, name
print('Captured CAD mesh, physical intent and solver output rechecked without engines or GUI.')
"""
            result = subprocess.run(
                [sys.executable, "-c", code, str(destination)],
                cwd=destination,
                env={**os.environ, "PATH": ""},
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("rechecked without engines", result.stdout)


if __name__ == "__main__":
    unittest.main()
