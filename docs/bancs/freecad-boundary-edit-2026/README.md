# Native FreeCAD boundary editing — 2026-09-10

[Static task guide](../../FREECAD_STATIC_TASK.md).

Development **0.1.0a3.dev4** fixes the pressure-editing interaction: selecting a
pressure loads its stored MPa value; updating and native Undo preserve the selected
object and refresh its value. Stable object names identify rows independently of
labels. Supports and empty selections cannot activate the pressure-update button.
The list context menu offers **Edit pressure…**, which focuses and selects the
pressure field. Applying a change still requires the explicit update button.

The installed FreeCAD qualification passes sixteen grouped checks, including the
previous first-run/static lifecycle cases and two new editing checks. It sends a
native mouse-context-menu event to the displayed list viewport, then clicks the
visible menu action with QtTest. This establishes native Qt event routing and
focus, not physical mouse hardware behavior on all Linux desktops.

The actual accepted CalculiX bar result has maximum nodal displacement error
4.762e-13 m against the affine reference. A second actual calculation rejects
results after the material changes. The GUI timer fires 149 times during the
accepted solve; this is not a general latency bound. The FreeCAD process exits
zero and its console contains no Python traceback.

`record.zip` retains the exact documented package, source hashes and executed
sources, FCStd inputs/results, STEP, mesh and raw CalculiX files, native screenshot,
console and [report](static-check.json). Verify `sha256sum -c SHA256SUMS` before
extracting. Both producer and packaged runtime/installation hashes were matched
against the final checkout. The manifest honestly records a dirty development
candidate based on d52a3d683f52b2b312f1ddde4d93bb19c3242c1d, not an official release.

Runtime: FreeCAD 1.1.3 / Qt 6 / Python 3.11 / host OCCT 7.8.1 on Linux x86-64;
separate Python 3.14.7 / Vinkulum 0.20 / OCCT 8.0.1; Gmsh HXT
5.0.0-git-91b4154 / OCCT 8.0.1 and CalculiX 2.21, two engine threads. Repeat using
`ci/freecad_static.py --recipe static-task` and the executable arguments documented
in the static task guide. Recorded absolute paths identify this experiment.

This does not establish general FEM accuracy, arbitrary topology persistence,
finished viewport polish or other-platform qualification. Assembly and single-solid
Motion checks were not rerun for this isolated static-panel change. No solver,
physical model, numerical tolerance or published release asset was changed.
