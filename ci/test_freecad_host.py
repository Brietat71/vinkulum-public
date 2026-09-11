"""Linux process-ownership regressions; no FreeCAD SDK or network required."""

import ctypes
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from freecad_host import destination_lock, require_idle, run_recorded


@unittest.skipUnless(sys.platform == "linux", "Linux host installer")
class HostOwnershipTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="vinkulum-host-test-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def lock_probe(self):
        return subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "from pathlib import Path; import sys; "
                    "from freecad_host import destination_lock; "
                    "ctx=destination_lock(Path(sys.argv[1])); ctx.__enter__(); ctx.__exit__(None,None,None)"
                ),
                str(self.root),
            ],
            cwd=Path(__file__).parent,
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )

    def test_competing_process_refused_and_exception_releases_lock(self):
        with (
            self.assertRaisesRegex(RuntimeError, "owner failure"),
            destination_lock(self.root),
        ):
            result = self.lock_probe()
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Another installer owns", result.stderr)
            raise RuntimeError("owner failure")
        self.assertEqual(self.lock_probe().returncode, 0)
        self.assertTrue((self.root / ".build.lock").is_file())

    def test_lock_symlink_refused(self):
        target = self.root / "untouched"
        target.write_text("original")
        (self.root / ".build.lock").symlink_to(target)
        with self.assertRaises(OSError), destination_lock(self.root):
            self.fail("symlink accepted")
        self.assertEqual(target.read_text(), "original")

    def test_unrecorded_launch_refused(self):
        with self.assertRaisesRegex(ValueError, "no process record"):
            require_idle({"steps": [{"status": "starting"}]})

    def test_command_success_and_failure(self):
        for code in (0, 7):
            with self.subTest(code=code):
                step = {}
                snapshots = []
                with (self.root / "log").open("w") as log:
                    command = [sys.executable, "-c", f"raise SystemExit({code})"]
                    if code:
                        with self.assertRaises(subprocess.CalledProcessError) as error:
                            run_recorded(
                                command,
                                self.root,
                                os.environ.copy(),
                                log,
                                step,
                                lambda snapshots=snapshots, step=step: snapshots.append(
                                    dict(step)
                                ),
                            )
                        self.assertEqual(error.exception.returncode, code)
                        self.assertNotEqual(step["status"], "passed")
                    else:
                        run_recorded(
                            command,
                            self.root,
                            os.environ.copy(),
                            log,
                            step,
                            lambda snapshots=snapshots, step=step: snapshots.append(
                                dict(step)
                            ),
                        )
                        self.assertEqual(step["status"], "passed")
                self.assertEqual(snapshots[0]["status"], "starting")
                self.assertIn("process_session", snapshots[1])
                require_idle({"steps": [step]})

    def test_orphan_in_another_group_still_prevents_resume(self):
        libc = ctypes.CDLL(None, use_errno=True)
        previous = ctypes.c_int()
        self.assertEqual(libc.prctl(37, ctypes.byref(previous), 0, 0, 0), 0)
        self.assertEqual(libc.prctl(36, 1, 0, 0, 0), 0)
        self.addCleanup(lambda: libc.prctl(36, previous.value, 0, 0, 0))
        # The command survives its session leader and changes group, as Ninja does.
        child = (
            "import json,os,sys; os.setpgrp(); "
            "print(json.dumps({'pid':os.getpid(),'group':os.getpgrp(),'session':os.getsid(0)}),flush=True); "
            "sys.stdin.read()"
        )
        parent = (
            "import subprocess,sys; "
            "p=subprocess.Popen([sys.executable,'-c',sys.argv[1]],stdout=subprocess.PIPE,text=True); "
            "print(p.stdout.readline(),end='',flush=True)"
        )
        owner = subprocess.Popen(
            [sys.executable, "-c", parent, child],
            start_new_session=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )
        pid = None
        try:
            # communicate closes stdin, so read the one bounded readiness record separately.
            import select

            self.assertTrue(select.select([owner.stdout], [], [], 10)[0])
            info = json.loads(owner.stdout.readline())
            pid = info["pid"]
            owner.wait(timeout=10)
            self.assertEqual(info["session"], owner.pid)
            self.assertNotEqual(info["group"], owner.pid)
            with self.assertRaisesRegex(ValueError, "still active"):
                require_idle(
                    {"steps": [{"status": "running", "process_session": owner.pid}]}
                )
            owner.stdin.close()
            os.waitpid(pid, 0)
            pid = None
            require_idle(
                {"steps": [{"status": "running", "process_session": owner.pid}]}
            )
        finally:
            if owner.poll() is None:
                owner.kill()
                owner.wait(timeout=10)
            if pid is not None:
                try:
                    os.kill(pid, 9)
                    os.waitpid(pid, 0)
                except ProcessLookupError:
                    pass
            for stream in (owner.stdin, owner.stdout):
                if stream is not None:
                    stream.close()


if __name__ == "__main__":
    unittest.main()
