# Native Assembly Motion analyses — development qualification

This source branch brings native FreeCAD assemblies into the existing **Motion
analysis** task. The extension's menu opens an included editable double-pendulum
Assembly, selects it and offers the same density, time-step, thread and engine
controls as the single-solid workflow. Revolute joint frames come from FreeCAD's
Assembly objects; the task omits the single-solid pivot/axis controls.

This is a development contribution. The published 0.1.0a2 ZIP remains the
single-solid release; build/install the extension from this branch to try this
workflow. FreeCAD is the primary desktop interface. No standalone Studio GUI or
competing workbench is introduced.

![Actual native Assembly Motion task and captured motion](assembly-analysis.png)

## Use the native workflow

1. Open **Vinkulum → Open double-pendulum assembly**.
2. With the Assembly selected, open **Vinkulum → Motion analysis**.
3. Choose the separate engine Python and calculation directory. To reproduce the
   reference, set density to 2700 kg/m³, duration to 0.5 s and step to 0.005 s.
4. Run the calculation, then use the slider or **Play captured motion**.
5. Save the `.FCStd` file. Reopen its Motion analysis and choose **Reopen last
   calculation** to display the captured result without starting an engine.

The document stores the analysis UUID, Assembly link, inputs and last calculation
folder. Each linked component receives a body UUID derived from the analysis UUID
and its native object name, preserving identities across repeated calculations
and document reopening; each calculation still receives a fresh run UUID. Labels
remain editable. The current same-document source link is checked by the task.
Changing the source between a solid and an Assembly closes the old task and
configures the appropriate native property/task controls on reopening.

Playback shows copies of both moving bodies and the fixed component. It does not
apply the computed motion to native design placements or connector frames.
Saving, closing, source invalidation or deleting a preview object removes the
whole temporary group and restores source visibility. Source-joint edits
invalidate playback as well as shape and inherited-placement edits.

Only independent preview features are recomputed when creating the display. A
whole-document recompute was observed to execute pending native joint changes
after editing/restoring a connector, changing source frames during a later
reopen. The observed coordinate change is only 1.39 × 10⁻¹⁷ m, alongside a
changed BREP hash, but it invalidates the exact captured-source signature. The retained counterexample fails on that intermediate implementation;
the final implementation preserves the source state across repeated reopen/delete
operations. The [upstream feature recompute API](https://github.com/FreeCAD/FreeCAD/blob/main/src/App/DocumentObjectPyImp.cpp)
provides this separate feature operation; the actual behaviour is qualified on
the runtime below.

## Validation and limits

The installed package is checked against its manifest, including both bundled
FCStd examples. The recipe uses the real menu action and task buttons, not a
replacement UI or mock engine. It exercises native analysis creation Undo/Redo,
save/reopen, two actual calculations with stable project/body identities,
capture reopening without a working engine, both displayed body centres,
fixed-part placement, source preservation, moved Assembly and edited connector
rejection, removal of either a moving or fixed preview object, source-type
changes and document close while a real barrier process is running.

The independent [verifier](../../../apps/freecad/verify_assembly_analysis.py) uses
NumPy and the original finite-section double-pendulum Lagrange reference from
[the conversion record](../freecad-assembly-2026/README.md). It checks request/STEP
and native-archive hashes, exact native NPZ/JSON sample agreement, stable
identities and both recalculated trajectories. At 5 ms, each trajectory has
maximum angle error 3.410 × 10⁻⁴ rad against that reference, below the declared
2 × 10⁻³ rad budget. The 16/32-subdivision reference gap is below 4.32 × 10⁻¹³ rad
(budget 10⁻¹⁰). Displayed-centre error is below 10⁻⁹ m. Its recorded environment
contains no FreeCAD, Vinkulum, Studio, OCCT, build123d or Qt imports.

The Assembly recipe passes 26 checks and exits normally. The canonical
single-solid extension qualification passes all 18 checks, and five rapid
example-close cycles pass on the same installed candidate. The mandatory `ci/local.sh` passes, with three
optional Exudyn integrations skipped. No solver algorithm or geometry/physics
admission tolerance is changed. The numerical kernel and retained mechanical
campaigns are unchanged from the conversion contribution's `--bancs` run.

This qualifies a flat native Assembly with one grounded solid, uniform moving
body density and Revolute joints, starting from rest under gravity −Z. The
bounded converter's limits still apply: at most sixteen moving single solids,
ten seconds and 2,001 native samples. The included fixture has two moving
rectangular solids and detached native connector frames. Nested assemblies,
general face-attachment persistence, joint limits, other joint types, closed-loop
qualification, contact, FEM and other operating systems remain outside the
result. Native mechanical status remains `NotAssessed`.

## Reproduce

Build and install the extension into an isolated FreeCAD user profile, following
[the extension guide](../../../apps/freecad/README.md). The recipe checks the
installed files and launches their worker in the explicitly selected external
engine environment. Use a fresh output directory and an empty dedicated FreeCAD
process. The recipe closes that process when finished; a session that already
contains documents is refused without closing those documents.

```sh
export VINKULUM_ASSEMBLY_ANALYSIS_CHECK=/tmp/freecad-assembly-analysis
export VINKULUM_ASSEMBLY_RECORD="$PWD/docs/bancs/freecad-assembly-2026/record.zip"
export VINKULUM_FREECAD_PYTHON=/absolute/path/to/engine/bin/python
FreeCAD "$PWD/apps/freecad/qualify_assembly_analysis.FCMacro"
/absolute/path/to/numpy-only/bin/python apps/freecad/verify_assembly_analysis.py \
  "$VINKULUM_ASSEMBLY_ANALYSIS_CHECK"
```

On headless Linux prefix FreeCAD with
`xvfb-run -a -s '-screen 0 1600x1100x24'`. Check both the FreeCAD exit code and
`report.json`; a successful interpreter exit alone is insufficient. The report,
actual screenshot, captured inputs/native results, producer sources and installed
candidate ZIP are retained in `record.zip`. Its manifest explicitly marks the
uncommitted qualification candidate based on `4862923`; it is not a new release.
The archive also carries a package with expanded installation instructions;
byte comparison confirms that only the installation text and its manifest entry
change, with all runtime and example bytes unchanged.
A separate actual-FreeCAD probe verifies preservation of a pre-existing document
when the qualification refuses its session.

```sh
cd docs/bancs/freecad-assembly-analysis-2026
sha256sum -c SHA256SUMS
unzip record.zip -d /tmp/retained-assembly-analysis
/absolute/path/to/numpy-only/bin/python \
  ../../../apps/freecad/verify_assembly_analysis.py \
  /tmp/retained-assembly-analysis/record
```

The host is the unchanged official FreeCAD 1.1.3 Linux x86-64 AppImage, Python
3.11.14, Qt 6 and OCCT 7.8.1. The engine is Vinkulum 0.20.0, Studio 0.6.0a2.dev6,
Python 3.14.7 and OCCT 8.0.1.0. The independent reference uses Python 3.14.7 and
NumPy 2.5.3. Runtime provenance and the official AppImage checksum are in the
[conversion record](../freecad-assembly-2026/README.md#provenance). The source,
models, reference and captured data are original Apache-2.0 contributions; the
separately installed FreeCAD runtime retains its upstream licence.
