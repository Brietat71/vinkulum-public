# Contributing to Vinkulum

Start with the [contributor projects](docs/CONTRIBUTOR_PROJECTS.md): concrete
first deliverables for CAD, dynamics, Pinocchio, numerical verification, Rust,
Qt and documentation. English and French are both welcome. A useful contribution
can be a reproducible defect, an independent physical reference, a proof, clearer
documentation or a measured performance improvement.

## Choose a small first result

Open an [issue](https://github.com/Brietat71/vinkulum-public/issues/new/choose)
with the problem, a minimal example, your reference and how success can be checked.
Check existing work before starting a substantial change. You do not need to
understand the whole kernel to contribute a regression part or improve a workflow.

For bugs, include the version, platform, steps, observed and expected results.
For numerical problems, include units, tolerances, time steps and the source of
the expected result. Keep secrets and confidential models out of public reports.

## Develop and validate

Follow the [root README](README.md) for the kernel or the
[FreeCAD extension guide](apps/freecad/README.md) for the desktop interface.
The standalone Studio GUI is retired. Its engine adapters remain maintained;
new interface contributions and desktop deliveries target FreeCAD on Linux.

```sh
PY="$VIRTUAL_ENV/bin/python" bash ci/studio.sh
```

CAD changes need the [adapted OCCT 8 environment](docs/STUDIO_CAD.md).
For a new solver, begin with a conversion contract and a reference case, as in
the [Pinocchio integration plan](docs/PINOCCHIO_INTEGRATION.md).

Prepare Lean 4.19.0 and Mathlib as described in
[the proof workspace](preuves/README.md), then run `ci/local.sh` before a pull
request. Mechanical changes also need `ci/local.sh --bancs` and relevant
counterexamples. GitHub Actions calls the same CI with `--bancs` and checks an
installed wheel outside the repository. External campaigns can require separate
solver installations. Report which checks ran and which could not run.

Performance proposals need a reproducible workload, an initial profile and an
accuracy comparison. Report transfer costs, memory use and regressions alongside
speedups. Keep benchmark inputs and measurement procedures available for review.

## Submit a reviewable change

Describe the concrete problem and resulting behaviour, then the evidence used
to check it. State supported physics and remaining limitations. A benchmark result
or a passing test does not establish a general certification claim.

Original implementations rely on literature, public interfaces, input models
and observations. Do not copy a third-party solver's code to implement Vinkulum.
Identify the provenance and licences of any proposed third-party files; keep
adaptations explicit and preserve attribution.

Unless explicitly stated otherwise, intentionally submitted contributions are
offered under Apache-2.0, as provided in section 5 of the licence. Submit only
content you are authorised to distribute.
