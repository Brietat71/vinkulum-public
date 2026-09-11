# Capture solids inside FreeCAD compounds

Development extension **0.1.0a3.dev7** reads mass properties from the single solid
inside a `Compound` or `CompSolid`, including nested singleton containers.
Previously, a Boolean-cut plate could pass the one-solid admission check and
then fail with `AttributeError: 'Part.Compound' object has no attribute
'CenterOfMass'` before launching its calculation.

Motion and static capture now unwrap containers for volume, centre of mass and
the complete centroidal inertia tensor. The source document, its STEP export,
geometry fingerprint and static face numbering remain unchanged. Containers
with extra faces, edges, vertices or another solid are refused before creating
a capture directory. Additional geometry is not silently discarded. This
change concerns the single-source motion and static capture paths; native
Assembly capture has its own component extraction contract.

The [retained qualification record](bancs/freecad-compound-capture-2026/README.md)
contains the original failure, analytic checks, actual worker results and native
static result persistence checks.

## Reproduce on Linux

Use the [qualified patched host](FREECAD_HOST.md), the
[separate engine](FREECAD_ENGINE.md), OCCT 8-enabled Gmsh HXT and CalculiX.
From the repository root:

```sh
python3 ci/freecad_static.py --recipe compound-capture \
  --freecad /absolute/path/to/FreeCAD-Vinkulum \
  --engine-python /absolute/path/to/engine/bin/python \
  --gmsh /absolute/path/to/gmsh --ccx /absolute/path/to/ccx \
  --output /tmp/vinkulum-compound-capture

python3 ci/freecad_static.py --recipe static-curved-result --rotated --compound \
  --freecad /absolute/path/to/FreeCAD-Vinkulum \
  --engine-python /absolute/path/to/engine/bin/python \
  --gmsh /absolute/path/to/gmsh --ccx /absolute/path/to/ccx \
  --output /tmp/vinkulum-compound-static
```

Each command requires a new output directory. The runner builds and installs the
extension ZIP in an isolated FreeCAD profile and retains its source hashes,
console output, native report and captured files. The first recipe captures
inputs without running a solver; the second runs the native task, mesher and
CalculiX and checks the reopened FEM result.

The capture reference is a 120 × 80 × 10 mm plate minus four through cylinders
of radius 3 mm. Signed analytic volumes and the parallel-axis theorem provide
an independent mass, centre and full inertia reference, including a translated
and rotated nested container. Element comparisons use relative tolerance
`1e-9` and absolute tolerance `1e-13` in their recorded SI quantity. Both motion
and static capture are checked, together with rejection of extra geometry.

The static reference uses a parametrically linked cylinder inside a native
`Part::Compound`, rotated in space, with zero Poisson ratio under axial tension.
Its displacement field is checked against the analytic solution, its result
must remain current after save/reopen, and a `1e-10` mm height edit must invalidate
it. This is a bounded regression, not a general FEM error estimate or a guarantee
for all imported topology.
