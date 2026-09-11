# Explicit geometry admission for native Assembly components

Development extension **0.1.0a3.dev8** applies the same single-solid admission
rule to native Assembly components as to single-source motion and statics.
A valid Solid, or nested singleton Compound/CompSolid wrappers, can be captured.
A container with an extra face, edge, vertex or another solid is rejected with
the component's label. Previously, Assembly extraction selected `Solids[0]` and
silently discarded additional non-solid geometry.

The conversion still composes the native parent/component placement on a copy.
It preserves the source design and reads volume, centre of mass and the complete
centroidal inertia tensor from the admitted solid. This does not add support for
nested mechanisms, mixed-dimensional bodies or additional joint types.

The [retained qualification record](bancs/freecad-assembly-geometry-2026/README.md)
contains the original counterexample, passing recipes, worker output and
recomputation diagnostic.

## Reproduce the native regression

Use the [qualified Linux FreeCAD host](FREECAD_HOST.md) and
[separate engine](FREECAD_ENGINE.md). From the repository root:

```sh
python3 ci/freecad_extension.py --recipe assembly-geometry \
  --freecad /absolute/path/to/FreeCAD-Vinkulum \
  --engine-python /absolute/path/to/engine/bin/python \
  --output /tmp/vinkulum-assembly-geometry

python3 ci/freecad_extension.py --recipe assembly-analysis \
  --freecad /absolute/path/to/FreeCAD-Vinkulum \
  --engine-python /absolute/path/to/engine/bin/python \
  --output /tmp/vinkulum-assembly-analysis
```

Each output directory must be new. The runner installs the extension in an
isolated native profile; the first recipe checks geometry admission and the
second exercises the real motion worker, playback, persistence and cancellation.

The geometry recipe compares a 10 × 20 × 30 mm box against independent analytic
volume and full centroidal volume-moment formulas. It composes a translated
parent rotated 37° around Z with a translated component rotated 23° around X.
Plain solids, singleton and nested wrappers are tested directly and via native
App::Link. Centre and inertia comparisons use relative tolerance `1e-10` and
absolute tolerance `1e-8` in mm and mm⁵ respectively; volume uses relative
`1e-12` in mm³. Source BREP and placements must remain unchanged during extraction.
Extra faces, edges, vertices and solids are rejected through both representations.

The complete double-pendulum Assembly additionally rejects extra geometry on
its grounded component before creating any capture files. Restoring that ground
geometry permits a new capture. Repeated snapshots without editing are identical.

## Native recomputation and exact identity

Restoring one component's geometry does not promise bit-identical placements of
all other components after FreeCAD recomputes and solves the Assembly. During
this qualification, a BREP location term on `Lower` changed from zero to
approximately `-2.19047482849133e-29 mm`. The native Assembly source calls its
solver during recomputation when `SolveOnRecompute` is enabled.

The diagnostic before/after snapshots and BREP files are retained. No tolerance,
rounding or normalization of nonzero placement values was introduced into
Vinkulum's production identity checks. A changed snapshot can therefore require
a new calculation, even when the native perturbation is physically negligible.
These tests establish the specified capture cases, not general topology or
bit-exact reversibility of FreeCAD's Assembly solver.
