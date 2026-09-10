# FreeCAD to CalculiX: captured-face experiment

FreeCAD is the primary interface; custom Studio GUI development is paused.
This source-only experiment transfers a native solid and explicit boundary faces
through the existing OCCT 8 / Gmsh HXT / CalculiX adapters, then imports native
FreeCAD FEM mesh and displacement objects. It is **not an interactive FEM command**
and is not included in the published FreeCAD 0.1.0a2 extension ZIP.

## Run on Linux

Use a dedicated FreeCAD 1.1.3 Qt 6 process and the separate
[qualified engine environment](FREECAD_ENGINE.md). Supply an OCCT 8-enabled HXT
Gmsh executable and CalculiX. The host's OCCT 7.8.1 libraries stay separate from
the engine's OCCT 8.0.1 libraries. The runner creates and closes documents and
exits FreeCAD; do not load the qualification macro into a working session.

```sh
python3 ci/freecad_static.py \
  --freecad /path/to/FreeCAD/AppRun \
  --engine-python /path/to/engine/venv/bin/python \
  --gmsh /path/to/occt8/gmsh \
  --ccx /path/to/ccx \
  --output /path/to/new-qualification-directory

/path/to/engine/venv/bin/python ci/check_freecad_static_mapping.py \
  /path/to/new-qualification-directory
```

Without DISPLAY, the runner uses Xvfb. Each run captures source hashes, executable
paths, STEP files, raw mesh/solver evidence, source and result FCStd files, and a
screenshot. It imposes a 270-second process timeout and cleans up its dedicated
process group. This is test-runner cleanup, not an interactive job manager.

## What is checked

The two references are a 120 × 20 × 15 mm bar with E = 210 GPa, Poisson ratio
zero, one fully fixed end and 2 MPa tensile pressure on the other. The second
case is rotated 37 degrees about (1, 2, 3) and translated (250, −100, 500) mm.
For this particular problem, the exact displacement is affine:
`u = R e_x (p/E) x_local`, and energy is `p² V / (2 E)`.
Every nodal displacement, the net reaction and the total strain energy are checked
against that independent reference. The macro also checks metre/millimetre
conversion, saving/reopening native results, and rejecting changed geometry before
adding result objects. This affine reference does not demonstrate convergence
for a general elastic problem.

Selected FreeCAD face numbers are used only to export individual trimmed STEP
sheets. The worker matches their geometry against entire mesh surface groups,
checks membership of all boundary nodes, then checks quadratic surface area and
centroid. It never assumes that Gmsh surface numbers match FreeCAD face numbers.
A separate replay deliberately renumbers mesh surfaces and requires exactly
unchanged support and load arrays. It rejects duplicate geometric matches,
a displaced face and incorrect area coverage, and revalidates raw solver archives.

## Limits

One top-level solid; small-strain linear isotropic elasticity; quadratic tetrahedra;
at most 32 selected faces and 6,000 imported result nodes. Only the two planar-face
reference cases above are qualified here. Membership samples and area/centroid
checks are numerical guards, not a proof of correspondence for arbitrary curved
or topologically complex faces. A coarse curved mesh may be rejected. No general
finite-element error bound, contact, nonlinear material or assembly qualification
is claimed.

The result is a snapshot linked to its source with the captured geometry hash.
Geometry is rechecked on import, but editing it later does not automatically
invalidate the already imported result. Preview hashes bind the request identity;
they are not a cryptographic attestation of solver output against hostile edits.
Interactive face selection, boundary-condition editing, progress/cancellation,
post-import stale-result handling and FEM task-panel integration remain to be built.

The [retained qualification](bancs/freecad-static-2026/README.md) records the measured
errors, runtime versions and actual artifacts.
