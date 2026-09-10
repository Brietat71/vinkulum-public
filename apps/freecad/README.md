# FreeCAD → Vinkulum → FreeCAD: integration experiment

This prototype tests using FreeCAD's existing parametric modelling interface
with Vinkulum's mechanics engine in a separate process. It edits a PartDesign
pad, captures its solid and physical properties, runs an explicit revolute
mechanism, and displays the returned native poses on a separate copy in FreeCAD.
The parametric source is unchanged by playback.

The [retained qualification](../../docs/bancs/freecad-bridge-2026/README.md)
contains four real desktop runs, editable `.FCStd` sources, STEP captures,
reopenable Vinkulum projects, native trajectory archives, screenshots and an
independent physical-pendulum reference. It is an executable feasibility
experiment; an installable workbench and a dedicated application remain possible
next interfaces for the same transport.

## Process and geometry contract

```mermaid
flowchart LR
    F[FreeCAD parametric source] --> C[Captured STEP and SI properties]
    C --> V[Separate Vinkulum Python process]
    V --> R[Validated native trajectory]
    R --> P[Captured-shape playback in FreeCAD]
```

The tested official **FreeCAD 1.1.3** AppImage uses Python 3.11.14 and OCCT 7.8.1.
The Vinkulum worker uses Python 3.14.7, Studio 0.6.0a2.dev6, kernel 0.20.0 and
OCCT 8.0.1.0. Each process keeps its own Python and native libraries. The bridge
removes inherited `PYTHONHOME` and library/plugin overrides before launching the
explicitly selected Vinkulum interpreter; the AppImage sets those for FreeCAD's
own runtime. FreeCAD already supports external Python extensions through its
[workbench interface](https://github.com/FreeCAD/FreeCAD-documentation/blob/main/wiki/Workbench_creation.md).

[bridge.py](bridge.py) runs inside FreeCAD and captures one valid solid in a
fresh directory. STEP bytes, a request UUID, body identity, local BREP and global
placement identify that capture. Volume, centre and full centroidal inertia are
also recorded from FreeCAD. Its millimetre quantities become SI using `10⁻⁹`
for volume and `ρ × 10⁻¹⁵` for the inertia integral in mm⁵.

[worker.py](worker.py) runs in the existing Studio CAD environment. The production
OCCT 8 importer reads the STEP and recomputes all physical properties. Volume,
mass, centre and every inertia component must agree with the FreeCAD capture
within scale-dependent absolute budgets (`εV`, `εm`, `εL`, `εmL²`, `ε = 10⁻⁸`;
L is the imported bounding-box diagonal). The normal mechanical adapter creates
the explicit world-origin, world-Y revolute joint and validates the native
trajectory archive. The selected allocation is two threads; CAD import and
mechanics run sequentially in this worker.

FreeCAD receives JSON poses in metres and row-major body-to-world matrices.
The bridge checks capture identity, increasing times, finite values, proper
rotations and the initial pose. It also rechecks the source geometry and its
global placement. Moving an enclosing Part can leave the local BREP unchanged;
that inherited placement is therefore part of the signature. Playback applies
the rigid transform about the captured centre to a separate `Part::Feature`.
It does not write the motion into the original pad or its design placement.
Replacing the input file also invalidates admission, even if a result carries
the hash of that replacement: identity is checked against the original capture.

The QProcess has a 60-second timeout and bounded diagnostic capture. This
experiment exercises one job at a time in a dedicated FreeCAD process. A
production host still needs application-wide CPU admission, ordinary workbench
controls and complete close/cancel lifecycle qualification.

## Reproduce

Install the [Studio CAD environment](../../docs/STUDIO_CAD.md), then use a
dedicated FreeCAD process. The qualification creates and closes its own document
and exits that process; it refuses to start if documents are already open.
Use a new output directory:

```sh
export VINKULUM_FREECAD_PROTOTYPE="$PWD/apps/freecad"
export VINKULUM_FREECAD_PYTHON=/absolute/path/to/studio-environment/bin/python
export VINKULUM_FREECAD_CAPTURE=/tmp/freecad-vinkulum-capture
FreeCAD "$VINKULUM_FREECAD_PROTOTYPE/qualify.FCMacro"
"$VINKULUM_FREECAD_PYTHON" apps/freecad/verify.py \
  "$VINKULUM_FREECAD_CAPTURE" --plot
```

`FreeCAD` above denotes the executable for a normal FreeCAD installation. On
headless Linux, prefix it with `xvfb-run -a -s '-screen 0 1600x1100x24'`.
The qualification used the extracted official AppImage's `AppRun`, with
`XDG_CONFIG_HOME`, `XDG_DATA_HOME` and `XDG_CACHE_HOME` pointing to separate
temporary directories. The release URL and checked checksum are in the record.

[qualify.FCMacro](qualify.FCMacro) changes the pad length from 800 to 1200 mm,
with native time steps of 5 and 2.5 ms for each length. It checks all displayed
centres against the returned poses, confirms that the source stays unchanged,
and rejects changed geometry, a moved parent, a corrupted rotation and a
replaced request. A Qt
timer records continued event-loop delivery during each external calculation.

[verify.py](verify.py) needs numpy, and matplotlib only for `--plot`. It imports
neither FreeCAD, Studio nor a geometry/physics engine. It checks input/result
hashes, exact agreement between the native archive and the replayed data, the
analytic box mass/inertia, pivot kinematics and convergence against an independent
scalar pendulum equation. It can also verify the extracted retained archive.

The normal Studio suite includes [two external-worker regressions](../studio/tests/test_freecad_bridge.py)
using the retained FreeCAD STEP. They require the Studio CAD environment but no
running FreeCAD: one checks an actual native trajectory against the independent
pendulum equation; the other changes the captured mass and checks rejection
before simulation. Run them directly with:

```sh
"$VINKULUM_FREECAD_PYTHON" -m unittest discover \
  -s apps/studio/tests -p test_freecad_bridge.py -v
```

## Scope of the result

The positive experiment is one uniform rectangular rigid body, one explicitly
defined revolute joint, two lengths and two seconds of motion from rest at 20°.
It establishes this parametric-edit → physical-capture → mechanics → native-pose
display path across the two runtimes. It does not establish general assembly
constraint conversion, persistent face attachment, arbitrary mechanisms,
multi-window scheduling, FEM or macOS/Windows packaging. STEP transfers evaluated
geometry; the original FreeCAD feature graph remains in `.FCStd`.

All bridge code, models, reference calculations and capture data are original
Vinkulum contributions under [Apache-2.0](../../LICENSE). The unchanged official
FreeCAD runtime remains available from its upstream project with its own licence.
