# Vinkulum Studio 0.6.0a1 — CAD and mechanism analysis

A local PySide6/VTK application for creating bodies and joints, manipulating
geometry, editing numerical properties, and defining motion laws and time-varying
loads. Vinkulum 0.19.0 runs in a separate process; Studio displays its positions,
orientations, velocities and joint coordinates.

Source version **0.6.0.dev4** introduced the [CAD-to-statics workspace](../../docs/STUDIO_CAD_MESHING.md):
generate an OCCT 8 / Gmsh mesh, select boundary faces in 3D, add supports,
pressure or total forces, save the captured study and open CalculiX. Meshing and
numerical admission run outside the GUI thread. The [input transport contract](../../docs/CALCULIX_NUMERIC_TRANSPORT.md)
revalidates the geometry actually sent to the solver's 20-character numeric fields.

The **0.6.0.dev3 source version** extended CalculiX studies to linear and quadratic
tetrahedra, including curved C3D10 geometry. The
[tetrahedral guide](../../docs/STUDIO_TETRAHEDRA.md) provides an independent
bending reference and the exact local Jacobian-positivity check. This is the
element contract now used by the CAD meshing workspace.

The **0.6.0.dev2 source version** adds [editable solid features](../../docs/STUDIO_CAD_HISTORY.md):
change an upstream dimension or placement, preview the regenerated part and apply
one undoable transaction. New CAD parts retain a feature graph; old BREP parts
remain readable as captured solid inputs. This workflow requires project schema 3.

The **0.6.0.dev1 source version** adds an **Articulated operators / Pinocchio**
workspace: edit a captured state, run a separate Pinocchio 4.1 worker, inspect
mass matrices, inverse/forward dynamics and Jacobians, export CSV, and reopen
results without the engine. See the [installation and operator guide](../../docs/PINOCCHIO_OPERATORS.md).
These workspaces are included in the **0.6.0a1 Linux preview**, along with
[examples to try in five minutes](packaging/EXAMPLES.md). New Gmsh, CalculiX
and Pinocchio computations require separately installed engines.

Studio 0.5.0 adds an experimental **Linear statics / CalculiX** workspace,
captured displacement/stress inspection, CSV export and checked reopening of
calculation folders. The solver remains a separately installed executable.
See the [FEM guide](../../docs/CALCULIX_INTEGRATION.md) for its supported domain.
Camera fitting also includes visible attachment points and ignores hidden
objects; annotation size is based on assembly span rather than world position.

Studio 0.4.2 fixed the normal launcher on PySide6 6.11.2 and qualified startup
through the actual application entry point. Version 0.4.1 introduced an English interface and preserved the precision of
numerical fields and compatibility with existing project files. The CAD workflow
introduced in 0.4.0 uses **OCCT 8.0.1** and an explicit adaptation of **build123d**.

## Install and run

From the repository root, with **Python 3.14**, Rust/Cargo, a C++17 compiler and
a system linker:

```sh
python3.14 -m venv .venv-studio
. .venv-studio/bin/activate
python -m pip install .
python -m pip install -e './apps/studio[test]'
python -m vinkulum_studio
```

Enable CAD in the same environment:

```sh
python ci/prepare_cad.py build/cad-sources
python -m pip install build/cad-sources/build123d-0.11.1 \
  build/cad-sources/ocpsvg-0.6.0 -e './apps/studio[cad,test]'
python -m vinkulum_studio
```

The CAD source preparation requires `patch`. Versions, adaptations and limits
are documented in [Studio CAD](../../docs/STUDIO_CAD.md).
A kernel **0.19.0** wheel matching your Python and platform can replace
`pip install .`. Linux wheels do not work on macOS.

On Ubuntu 24.04, install the desktop prerequisites with:

```sh
sudo apt-get install libegl1 libgl1 libgl1-mesa-dri \
  libfontconfig1 libfreetype6 libdbus-1-3 libglib2.0-0t64 \
  libx11-6 libx11-xcb1 libxkbcommon0 libxkbcommon-x11-0 \
  libxcb1 libxcb-cursor0 libxcb-icccm4 libxcb-keysyms1 libxcb-image0 \
  libxcb-randr0 libxcb-render0 libxcb-render-util0 libxcb-shape0 \
  libxcb-shm0 libxcb-sync1 libxcb-util1 libxcb-xfixes0 libxcb-xkb1
# Automated desktop tests also require:
sudo apt-get install xvfb xauth
```

An OpenGL driver is required. PySide6 **6.11.2** and VTK **9.7.0** have their own
licences; original Studio code is Apache-2.0. Check the platform and version in
[public releases](https://github.com/Brietat71/vinkulum-public/releases).
An older Apple Silicon DMG predates the current CAD/UI changes. It requires
macOS 14 and uses ad-hoc signing without Apple notarisation. Qualification of
the current CAD release on macOS ARM64 is still pending.

### Local development and Linux packaging

Restart Studio after a Python/UI edit. The editable installation uses the updated
sources and reuses the compiled native kernel. Rebuild the extension after Rust
changes; package an archive when preparing a delivery.

To build a Linux archive from an environment containing the current native
kernel, Studio CAD and PyInstaller:

```sh
python -m pip install 'pyinstaller==6.22.2'
# First build the external mesher with ci/build_mesher.py; see the CAD meshing guide.
VINKULUM_BUNDLE_GMSH=/absolute/path/to/gmsh-install/bin/gmsh \
  PY=$(command -v python) bash ci/linux_bundle.sh
```

This runs locally. It builds, archives and extracts the application, then tests
the extracted executable outside the repository: real CAD, STEP round-trip,
CAD feature regeneration, native simulation, OpenGL rendering, Gmsh/OCCT 8,
CalculiX, and the delivered CAD/FEM/Pinocchio captures. Without a display it
uses Xvfb. `dist/linux/` receives the verified archive, `SHA256SUMS`, build
provenance and the `check/` report. Set `VINKULUM_LINUX_OUT` for another destination.

Extract `Vinkulum-Studio-0.6.0a1-linux-x86_64.tar.gz`, then run
`./Vinkulum\ Studio/Vinkulum\ Studio`. Keep `_internal` beside the executable.
The qualified 0.4.0 baseline used Linux x86-64, glibc 2.39, X11 and Mesa OpenGL
on Ubuntu 24.04; each new package must pass its own extracted-binary checks.
Older glibc, Linux ARM64 and native Wayland are not qualified. The system still
provides graphics libraries and its OpenGL driver. See
[Linux installation](packaging/INSTALLATION-LINUX.txt).

## Workbench

Three workspaces organise the interface: **Model**, **Simulate** and **Inspect**.
The viewport occupies the main area. Browser, Inspector, diagnostics and results
panels can be resized, detached and toggled in **View → Panels**. The launcher
preserves the theme and layout. **Restore layout** restores the active workspace.

- **Ctrl/Cmd+K** searches commands; **S** opens tools for the selection. File
  and editing shortcuts follow the platform.
- **Run → Articulated operators · Pinocchio…** (**Ctrl+Shift+P**) opens a
  captured-state study. **Ctrl+Return** evaluates supported tree operators;
  this workspace does not integrate a trajectory.
- **Create → Edit CAD features…** (**Ctrl+Shift+H**) edits the selected part's
  feature parameters. Preview a complete regeneration before applying it to the model.
- The XY reference grid fades towards its edges and can be hidden. **Views**
  provides perspective and orthographic projections. Joint rings and axes are
  symbols, not added mechanical volume; attachment lines appear on selection.
- The Browser filters names, types and IDs. Its context menu can isolate or
  hide objects without changing the simulation.
- The Inspector separates X/Y/Z fields. Compact displayed numbers preserve the
  complete unedited values; focusing a field exposes its full precision.
- **Select**, **Move** and **Rotate** control the manipulator. The camera widget
  selects orientations directly in the viewport.
- **Inspect** retains up to eight runs within a 128 MiB array budget. Selecting
  a run displays its own captured, read-only model.
- Click a curve to select a sample; use the wheel to zoom, Shift-drag to pan,
  double-click to fit, and arrow keys to step. Time, scene and table stay aligned.
- **Compare with** overlays a dashed reference and describes differences in
  models and settings. Series keep native time grids; objects are associated
  by stable identity rather than list position.

The interface is an early workbench iteration. Its
[design criteria and references](../../docs/STUDIO_GUI_2026.md) do not imply
feature parity with a mature CAD/CAE suite.

## Build and simulate

1. Start a new project or load an example: G0 pendulum, double pendulum,
   slider-crank, or a body under a time-varying force.
2. Add a box, cylinder or sphere, or use **CAD** for exact solid modelling.
   Select a body, edit dimensions, mass and pose, then apply its properties.
   Homogeneous inertia is calculated; an explicit matrix can be supplied.
3. Move or rotate it with the manipulator. A completed operation creates one
   history entry; undo/redo restores the document.
4. Add a revolute, spherical, prismatic or fixed joint between bodies or to
   ground. A world point and axis initialise two independently editable frames.
5. Optionally configure a revolute/prismatic motion law or a body load.
   Dialogs support constant, linear and tabulated laws.
6. Resolve diagnostics, set duration and time step, then run the design.
7. Inspect the result animation, time slider and curves; export CSV with
   provenance. Return to the design to edit its model.

Units are SI: m, kg, s, N, N·m and kg·m². Orientation fields use degrees with
`Rz × Ry × Rx`; documents store rotation matrices. Angular laws use radians.
A cylinder's axis is local Z.

**Moving a body does not solve its constraints.** Inconsistent anchors remain
visible and prevent a run; native initialisation rejects silently relocating
bodies. It can initialise velocities prescribed by motion laws. The grid does
not represent ground contact.

Force and moment components are in the world frame. The application point is
body-local, so its lever arm rotates with the body. Orange force and purple
moment arrows indicate direction; length is symbolic. Tabulated laws interpolate
linearly and extend endpoint values. No Python expression is evaluated.

## Documents and results

Objects have persistent UUIDs. JSON uses `format: vinkulum-studio-project` and
an immutable project. Projects without CAD use schema 1; those with BREP use
schema 2. G0 files remain importable. English labels preserve existing joint
and law identifiers in saved files. Deleting a body retains broken joint/load
references so they can be diagnosed and repaired.

Saves and exports replace a destination after the temporary file has been fully
written and `fsync` completed; this is not a universal power-loss guarantee.
Each run captures its own project. The design remains editable while it runs;
result properties and geometry belong to the captured, read-only project.
Failure, crash or cancellation preserves the previous result. One simulation
worker is allowed; closing the application terminates it.

Temporary results are bounded NPZ archives validated without pickle: digest,
project identity, headers before allocation, shapes, finite values, chronology
and rotation matrices. Provenance contains versions and the native binary hash.
CSV exports native samples, including the initial state, units and captured
project. Revolute angles use the principal range `[-π, π]`, not a turn counter.
Playback selects samples without dynamic interpolation. History remains in
session memory within its budget; JSON saves the design and CSV retains data.

## Domain and qualification

Limits: 32 bodies, 64 joints, 128 loads, 20,000 requested steps, 100,000
body/sample pairs and 64 MiB per result archive. These budgets do not guarantee
convergence or impose an OS memory limit. Workers are not security sandboxes.
Contact, flexible bodies, URDF import, collaboration and automatic mechanism
synthesis are not exposed by this Studio workflow. STEP import accepts one
solid as specified in the CAD contract.

A trajectory's scientific status remains **`NotAssessed`**. Integrity checks
and tested physical references are not error bounds for arbitrary models.
See the [3D qualification record](../../docs/STUDIO_3D.md).

```sh
PY="$VIRTUAL_ENV/bin/python" bash ci/studio.sh
```

Linux tests use a real Qt/X11/OpenGL context under Xvfb. On macOS, run from a
desktop session with screen access. Human qualification is still needed for
comfort, shortcuts and device-specific gestures on each platform.
