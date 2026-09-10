"""The distributed scientific example remains readable without any solver/UI."""

import subprocess
import sys
import unittest
from pathlib import Path


class CapturedOperators(unittest.TestCase):
    def test_saved_example_in_a_fresh_engine_free_process(self):
        root = Path(__file__).resolve().parents[3]
        example = root / "examples/studio/articulated/double-pendulum"
        source = """
import sys
from pathlib import Path
from vinkulum_studio.articulated_result import load_operators
p = Path(sys.argv[1])
before = {f.name: (f.read_bytes(), f.stat().st_mtime_ns) for f in p.glob('*.json')}
result = load_operators(p)
assert result.state['q'] == [0.2, -0.3]
assert result.report['nq'] == result.report['nv'] == 2
assert result.report['scientific_status'] == 'NotAssessed'
assert len(result.display_project.bodies) == 2
assert not any(name.split('.')[0] in {'pinocchio', 'PySide6', 'vtk', 'vtkmodules', 'OCP', 'vinkulum'} for name in sys.modules)
assert before == {f.name: (f.read_bytes(), f.stat().st_mtime_ns) for f in p.glob('*.json')}
"""
        completed = subprocess.run(
            [sys.executable, "-c", source, str(example)],
            cwd=example,
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
