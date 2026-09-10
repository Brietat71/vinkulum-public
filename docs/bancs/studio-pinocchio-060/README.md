# Studio 0.6.0.dev1 — installed-package qualification

These checks ran locally on Linux x86-64 using an installed Studio wheel outside
the source tree, in a CAD/GUI environment with a separate installed Pinocchio
4.1.0 worker. The [manifest](qualification.json) records the wheel SHA-256,
the exact packaged Python source hashes, evidence hashes and test outcomes.
This qualifies the source workflow and Python wheel, **not a standalone archive**.

| Check | Result | Evidence |
|---|---|---|
| Full Studio desktop suite with OCCT 8, native Vinkulum and external ccx/Pinocchio | 73 passed, 10 numerical-engine cases skipped in the GUI environment | [Desktop log](studio-tests.log) |
| Pinocchio numerical suite in its separate installed environment | 11 passed, none skipped | [Operator log](pinocchio-tests.log) |
| Saved public example in a fresh process without importing any GUI or engine | 1 passed | [Archive log](saved-example-tests.log) |
| Real articulated window, captured two-link state and lower-body Jacobian | Passed | [Recipe report](recipe.json), [application capture](../../assets/studio-pinocchio.png) |

The desktop suite includes five Pinocchio GUI/controller tests: keyboard
activation, full-precision inputs and CSV, immutable captures, engine-free
reopening, unsupported models, invalid input preservation, close vetoes,
cancellation/timeout and rejection of inconsistent output. The numerical suite
includes independent Lagrange references, not just agreement between operators
from the same engine. Its shared admission test also runs in the desktop suite;
the rows above are separate executions, not a claim of 85 distinct test cases.

The saved example was added after the full desktop pass and checked separately;
no packaged runtime source changed afterward. Public CI includes it in subsequent
full passes. `qualification.json` fingerprints identify exactly the Python files
that were installed and tested, without claiming a clean Git build retroactively.

The public [Studio 0.5.0 Linux release](https://github.com/Brietat71/vinkulum-public/releases/tag/studio-v0.5.0-linux)
predates this Pinocchio workspace. Engine bundling, macOS qualification,
trajectory integration, general derivative verification and comparisons with
native-kernel observables remain separate work. Results stay `NotAssessed`.
