"""Real external processes must not load the frozen desktop's environment."""

import json
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication
from test_meshing import request
from vinkulum_studio.examples3d import double_pendulum
from vinkulum_studio.meshing_controller import MeshingController
from vinkulum_studio.pinocchio_controller import PinocchioController
from vinkulum_studio.static_controller import StaticController
from vinkulum_studio.static_window import tension_example


@unittest.skipUnless(os.environ.get("VINKULUM_3D_TESTS") == "1", "Requires Qt")
class ExternalEngineEnvironment(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_actual_engine_processes_restore_caller_libraries_without_mutating_gui(
        self,
    ):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for frozen, original in (
                (True, "caller-libs"),
                (True, ""),
                (False, "caller-libs"),
            ):
                for kind in ("gmsh", "ccx", "pinocchio"):
                    with self.subTest(frozen=frozen, original=original, engine=kind):
                        folder = root / f"{frozen}-{original}-{kind}"
                        folder.mkdir()
                        record = folder / "environment.json"
                        engine = folder / "external-engine"
                        keys = (
                            "LD_LIBRARY_PATH",
                            "PYTHONPATH",
                            "PYTHONHOME",
                            "VIRTUAL_ENV",
                            "VINKULUM_TEST_MARKER",
                        )
                        engine.write_text(
                            f"#!{sys.executable}\nimport json,os\nfrom pathlib import Path\n"
                            f"Path({str(record)!r}).write_text(json.dumps({{k:os.environ.get(k) for k in {keys!r}}}))\n"
                            "raise SystemExit(2)\n"
                        )
                        engine.chmod(0o755)
                        environment = {
                            "LD_LIBRARY_PATH": "/not-a-real-desktop-library-directory",
                            "LD_LIBRARY_PATH_ORIG": str(root / original)
                            if original
                            else "",
                            "PYTHONPATH": "/not-an-engine-module-directory",
                            "PYTHONHOME": "/not-an-engine-python-home",
                            "VIRTUAL_ENV": "/not-an-engine-venv",
                            "VINKULUM_TEST_MARKER": "inherited",
                        }
                        with (
                            patch.dict(os.environ, environment),
                            patch.object(sys, "frozen", frozen, create=True),
                        ):
                            before = dict(os.environ)
                            if kind == "gmsh":
                                controller = MeshingController()
                                controller.start(
                                    request(), folder / "run", executable=engine
                                )
                            elif kind == "ccx":
                                controller = StaticController()
                                controller.start(
                                    tension_example(), folder / "run", executable=engine
                                )
                            else:
                                controller = PinocchioController()
                                controller.start(
                                    double_pendulum(),
                                    {},
                                    folder / "run",
                                    interpreter=engine,
                                )
                            try:
                                deadline = time.monotonic() + 10
                                while (
                                    controller.busy
                                    if kind == "gmsh"
                                    else controller.process is not None
                                ) and time.monotonic() < deadline:
                                    QTest.qWait(10)
                                self.assertIsNone(controller.process)
                                self.assertTrue(
                                    record.is_file(),
                                    "External interpreter could not start independently",
                                )
                                actual = json.loads(record.read_text())
                                expected_path = (
                                    (environment["LD_LIBRARY_PATH_ORIG"] or None)
                                    if frozen
                                    else environment["LD_LIBRARY_PATH"]
                                )
                                self.assertEqual(
                                    actual,
                                    dict.fromkeys(keys[:1], expected_path)
                                    | dict.fromkeys(keys[1:4])
                                    | {keys[4]: "inherited"},
                                )
                                self.assertEqual(dict(os.environ), before)
                            finally:
                                controller.shutdown()


if __name__ == "__main__":
    unittest.main()
