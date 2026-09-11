# FreeCAD engine without Qt or VTK — 2026-09-11

Adapter package **0.6.1.dev1** retains the `vinkulum_studio` import namespace used
by FreeCAD workers, but its required dependencies no longer include PySide6 or
VTK. The optional `legacy-desktop` extra preserves the dependency set needed to
reproduce archived controller tests. This is not a new Studio GUI release.

The Linux installer creates a new environment, imports the CAD/mechanics/FEM
adapters outside the checkout and refuses an environment in which PySide6,
shiboken6, vtk or vtkmodules can be found. It still checks the independent motion
reference and rejects a corrupted captured mass. Existing environments are never
modified by this recipe.

The retained installation uses Python 3.14.7, kernel 0.20.0 and OCCT binding
8.0.1.0. All 62 installed adapter Python modules match the reviewed source bytes.
A negative probe against the older Qt-equipped environment is refused before a
headless success report is written. No GUI libraries were removed from that
existing environment.

Actual FreeCAD 1.1.3 / Qt 6 on Linux uses the new external engine successfully:

- Single-solid motion extension: 18 grouped checks, including native document
  lifecycle and the independent pendulum reference.
- Native Assembly analysis: 26 grouped checks, including repeated captures,
  native joint changes and temporary playback cleanup.
- Static task: 18 grouped checks, including actual Gmsh/CalculiX execution,
  boundary controls, mouse/keyboard menus, cancellation and late-result refusal.
  Maximum displacement error for the existing affine bar is
  `4.761904763738946e-13 m`.

FreeCAD's own Qt remains in its host process. The new environment contains the
backend package's historical GUI source files but has no desktop dependencies;
this change does not claim to split every package or rename its public imports.
The separate Pinocchio CI environment remains separately qualified. The existing
CI also installs the optional desktop dependencies explicitly for retained
headless controller regressions; it does not build or launch a Studio desktop.

`record.zip` retains the installer and reports, source hashes, native macros and
scripts, installed extension ZIPs, native documents and raw solver captures.
Use `sha256sum -c SHA256SUMS` and inspect the scoped reports. These tests do not
fix the separately reviewed BREP save/reopen fingerprint issue, qualify other
platforms or establish new general FEM/trajectory accuracy guarantees.
