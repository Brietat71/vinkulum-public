"""Real Linux processes: host death kills the guardian and all static children."""

import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def alive(pid):
    path = Path(f"/proc/{pid}/stat")
    return path.exists() and path.read_text().split()[2] != "Z"


@unittest.skipUnless(sys.platform == "linux", "Linux static guardian")
class GuardianTests(unittest.TestCase):
    def test_parent_death_kills_worker_and_grandchild(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            guard = directory / "static_guard.py"
            shutil.copyfile(ROOT / "apps/freecad/static_guard.py", guard)
            (directory / "static_worker.py").write_text(
                "import os,subprocess,sys,time\nfrom pathlib import Path\n"
                'p=subprocess.Popen([sys.executable,"-c","import time; time.sleep(60)"])\n'
                'Path("children").write_text(f"{os.getpid()} {p.pid}")\ntime.sleep(60)\n'
            )
            host = subprocess.Popen(
                [
                    sys.executable,
                    "-c",
                    (
                        "import os,subprocess,sys,time\nfrom pathlib import Path\n"
                        'p=subprocess.Popen([sys.executable,"static_guard.py",str(os.getpid())],start_new_session=True)\n'
                        'Path("guardian").write_text(str(p.pid))\ntime.sleep(60)\n'
                    ),
                ],
                cwd=directory,
            )
            group = None
            try:
                deadline = time.monotonic() + 10
                while not (directory / "children").exists():
                    self.assertLess(time.monotonic(), deadline)
                    time.sleep(0.02)
                group = int((directory / "guardian").read_text())
                children = [
                    int(p) for p in (directory / "children").read_text().split()
                ]
                self.assertTrue(all(alive(p) for p in [group, *children]))
                host.kill()
                host.wait(timeout=5)
                while any(alive(p) for p in [group, *children]):
                    self.assertLess(time.monotonic(), deadline)
                    time.sleep(0.02)
            finally:
                host.kill()
                host.wait(timeout=5)
                if group is None and (directory / "guardian").exists():
                    group = int((directory / "guardian").read_text())
                if group:
                    try:
                        os.killpg(group, signal.SIGKILL)
                    except ProcessLookupError:
                        pass

    def test_worker_exit_code_is_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            guard = directory / "static_guard.py"
            shutil.copyfile(ROOT / "apps/freecad/static_guard.py", guard)
            (directory / "static_worker.py").write_text("raise SystemExit(7)\n")
            result = subprocess.run(
                [sys.executable, str(guard), str(os.getpid())],
                start_new_session=True,
                timeout=10,
                check=False,
            )
            self.assertEqual(result.returncode, 7)


if __name__ == "__main__":
    unittest.main()
