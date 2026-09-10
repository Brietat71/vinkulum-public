# Installed FreeCAD static task — 2026-09-10

[Workflow, installation and scope](../../FREECAD_STATIC_TASK.md).

![Actual installed native static task](static-task.png)

`record.zip` retains the exact development package, source hashes and executed
sources, input/result FCStd files, STEP captures, raw mesh/CalculiX outputs, console
logs, the task report, engine inventory, guardian test output and the existing
Motion report from the **same byte-identical package**. Verify with
`sha256sum -c SHA256SUMS`. The package manifest identifies **0.1.0a3.dev1**. Its provenance deliberately
records an uncommitted composition candidate, not a new official extension release.
Local absolute paths record the experiment; they are not portable defaults.

Runtime: FreeCAD 1.1.3, Qt 6, Python 3.11, host OCCT 7.8.1; separate Python 3.14.7
engine with Vinkulum 0.20 and OCCT 8.0.1; Gmsh HXT 5.0.0-git-91b4154 with OCCT 8.0.1;
CalculiX 2.21. Linux x86_64, two engine threads.

Eleven installed-task checks pass: native menu without a workbench switch; face
selection and boundary buttons; native boundary Undo/Redo; pressure editing and
units; clicking the displayed native Close button; input save/reopen with source
face links; real asynchronous statics and native coloured results; saved-result
reopening and invalidation without an open task; refusing changed face geometry;
cancellation and document closure terminating a worker plus its spawned child;
and refusing the late result of a real solve after material inputs change.
The report groups some related assertions into a single check.

The accepted reference is the 120 × 20 × 15 mm bar in 2 MPa tension with E = 210 GPa
and Poisson ratio zero. Maximum nodal displacement error against the exact affine
solution is 4.762e-13 m. The GUI timer executes 151 ticks during that solve; this
shows event-loop activity, not a general latency or responsiveness bound. A second
real calculation is deliberately made stale and creates no result object.
The wider transfer/reference budgets remain documented in the
[base experiment](../freecad-static-2026/README.md).

The cancellation fixture is an executable that spawns a sleeping child. It checks
the actual Qt-created process group used by the task; it does not claim cancellation
at every internal Gmsh or CalculiX phase. Two separate standard-library tests use
real Linux processes to check that the guardian propagates worker exit status and
kills its own group when its parent disappears. They run in `ci/local.sh` without
FreeCAD, CAD or solver dependencies.

All eighteen existing native single-solid Motion checks and all twenty-six native
Assembly Motion checks also pass on the exact same installed package. All three
FreeCAD processes exit zero and all three native consoles contain no Python
traceback. The static task runner also rejects Python tracebacks from callbacks
that FreeCAD catches internally. The Assembly trajectories independently pass
`verify_assembly_analysis.py`; their requests, native arrays, report and reference
verification are retained under `assembly/` in the same archive.

This final composition retains both the native Assembly menu/example and the new
static task. It was requalified after merging the Assembly contribution and
versioning the common package. The earlier isolated static-task qualification is
preserved in repository history; this record covers the combined development build.

Intermediate qualification exposed an invalid unit arithmetic expression in the
test and an invalid integer conversion of a Qt button enum in the host. Both were
corrected. The visible-button test waits for the normal GUI event loop to display
the task before clicking it; it does not replace the click with a direct method call.
No numerical acceptance budget was relaxed.

This qualifies the stated bar workflow and lifecycle cases, not arbitrary geometry,
nonlinear analysis, general topology persistence, FEM convergence or a release on
other platforms. The original 0.1.0a2 public release files are unchanged.
