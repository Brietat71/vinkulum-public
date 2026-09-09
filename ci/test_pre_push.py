"""Exercise real local pushes: foreign Git repositories and CI failure propagation."""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def check_hook(hook):
    env = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    with tempfile.TemporaryDirectory(prefix="vinkulum hook ") as directory:
        root = Path(directory)
        checkout, linked, foreign, remote = (
            root / name
            for name in ("checkout", "linked checkout", "dependency", "remote")
        )
        env["VINKULUM_TEST_FOREIGN"] = str(foreign)

        def git(where, *args, success=True):
            result = subprocess.run(
                ["git", "-C", str(where), *args],
                env=env,
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            if success and result.returncode:
                raise AssertionError(f"git {args}: {result.stderr}")
            return result

        git(root, "init", "--bare", str(remote))
        git(root, "init", "-b", "main", str(checkout))
        git(root, "init", str(foreign))
        git(
            foreign,
            "config",
            "remote.origin.url",
            "https://example.invalid/mathlib.git",
        )
        git(checkout, "config", "user.name", "Vinkulum hook test")
        git(checkout, "config", "user.email", "hook-test@example.invalid")
        git(checkout, "config", "commit.gpgsign", "false")
        git(checkout, "config", "core.hooksPath", "ci/hooks")
        hooks = checkout / "ci/hooks"
        hooks.mkdir(parents=True)
        shutil.copyfile(hook, hooks / "pre-push")
        (hooks / "pre-push").chmod(0o755)
        # A foreign repository query represents Lake's dependency lookup.
        ci = checkout / "ci/local.sh"
        ci.write_text("""#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
test "$(git -C "$VINKULUM_TEST_FOREIGN" config remote.origin.url)" = 'https://example.invalid/mathlib.git'
test "${VINKULUM_TEST_FAIL:-0}" = 0
""")
        ci.chmod(0o755)
        git(checkout, "add", "ci")
        git(checkout, "commit", "-m", "Test fixture")
        git(checkout, "remote", "add", "origin", str(remote))
        git(checkout, "worktree", "add", "-b", "linked", str(linked))
        for work in (checkout, linked):
            name = "normal" if work == checkout else "worktree"
            git(work, "push", "origin", f"HEAD:refs/heads/{name}")
            env["VINKULUM_TEST_FAIL"] = "1"
            failed = git(
                work,
                "push",
                "origin",
                f"HEAD:refs/heads/{name}-rejected",
                success=False,
            )
            assert failed.returncode != 0, "CI failure did not reject the push"
            assert (
                git(
                    remote,
                    "show-ref",
                    "--verify",
                    f"refs/heads/{name}-rejected",
                    success=False,
                ).returncode
                != 0
            )
            env.pop("VINKULUM_TEST_FAIL")
    print(
        "Pre-push: normal and linked checkouts isolate dependencies and reject CI failures."
    )


if __name__ == "__main__":
    check_hook(
        Path(sys.argv[1])
        if len(sys.argv) > 1
        else Path(__file__).parent / "hooks/pre-push"
    )
