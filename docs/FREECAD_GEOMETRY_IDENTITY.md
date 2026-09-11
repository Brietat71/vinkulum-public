# Static geometry identity and host persistence

Development candidate, not yet a released extension. New static analyses use
schema 2 and the `brep-v1-zero-checked0-local-world-v1` fingerprint. Existing schema-1
analyses keep the original raw BREP/world-placement rule. There is no Boolean
shape-equivalence fallback and no numerical tolerance in fingerprint admission.

The fingerprint covers the cleaned, root-local BREP representation and the
object's global placement. It ignores triangulation used for rendering and the
`Checked` bookkeeping bit in each recognized OCCT V1 topology record. It retains
geometry numbers, tolerances, topology order, orientation and the other six
record flags. Literal numeric `-0` tokens in geometry records and signed zero in the global
placement are normalized; topology-link orientation signs remain untouched; nonfinite
placement values and unknown fingerprint/schema versions are rejected. This is
identity of a defined serialized representation, not universal mathematical
identity of physical shapes or bit-exact persistence of every original double.
The V1 recognizer processes BREP emitted by Part; it is not an arbitrary BREP
file validator. Broader CAD-domain qualification remains in progress.

Each capture records its fingerprint kind. Import uses that same kind and the
captured request hash. Reading an old capture without this field uses the legacy
rule. An old analysis is migrated when all its boundary rows have been removed
and the user selects new boundary faces. That transaction updates the schema and
clears the captured-input fingerprint. It retains the old result as historical,
out-of-date data; a new calculation is required. Undo/Redo preserves these states.

## Host requirement for stable curved-shape reopening

The qualified Linux host is FreeCAD 1.1.3 with the
[stream-format patch](../ci/patches/freecad-1.1.3-brep-floatfield.patch) and
[placement round-trip patch](../ci/patches/freecad-1.1.3-placement-roundtrip.patch). This host
uses OCCT 7.8.1; Vinkulum's separate engine continues to use OCCT 8. The patched
host is a source build, not a published Vinkulum binary or an upstream release.
The extension does not silently change the user's document preferences.

Stock FreeCAD 1.1.3 can write ASCII BREP through a stream that inherited fixed
floating-point formatting. On the tested curved shapes this changes serialized
coordinates relative to a fresh export. The patch clears the inherited
floatfield for this export and restores the caller's flags on exit, including
exceptional exit. With the stock writer, a reopened result can therefore become
out of date; the fingerprint must not mask this geometry change to preserve it.

## Reproduce the qualified host build on Linux

The [automated host recipe](FREECAD_HOST.md) performs the pinned clone, patch
checks, locked SDK installation and build, and writes a launcher and provenance
record. Its full build, resume and native qualification are retained in the
[host build evidence](bancs/freecad-host-build-2026/README.md). The manual steps
below document the same upstream source and patch requirements.

Prerequisite: Pixi 0.80.0, Git, and enough space for the full isolated SDK/build.
The qualified Linux x86-64 Pixi executable SHA-256 is
`387a2d3052e656f61ccf735e6750255451366f45635a2da09116b1f8394b2837`.
Use the upstream FreeCAD repository and its committed lock file; do not update
the lock file to resolve an old-format warning.

From the Vinkulum checkout, select a new directory outside it:

```sh
set -eu
vinkulum_checkout="$PWD"
freecad_checkout="/tmp/freecad-vinkulum-host"
git clone --depth 1 --branch 1.1.3 https://github.com/FreeCAD/FreeCAD.git "$freecad_checkout"
cd "$freecad_checkout"
test "$(git rev-parse HEAD)" = 145529fe741292ff0b3977a01195bf0247425794
printf '%s  %s\n' f9a9c1f90cd7aa37655c593793c67dacec3e7b532b2403a29ff2575f013947f9 pixi.lock | sha256sum -c -
git submodule update --init --recursive
pixi install --locked
git apply --check "$vinkulum_checkout/ci/patches/freecad-1.1.3-brep-floatfield.patch"
git apply "$vinkulum_checkout/ci/patches/freecad-1.1.3-brep-floatfield.patch"
git apply --check "$vinkulum_checkout/ci/patches/freecad-1.1.3-placement-roundtrip.patch"
git apply "$vinkulum_checkout/ci/patches/freecad-1.1.3-placement-roundtrip.patch"
CFLAGS= CXXFLAGS= pixi run --locked cmake --preset conda-linux-release '-DCMAKE_JOB_POOLS=compile_jobs=4;link_jobs=1'
pixi run --locked cmake --build build/release --parallel 4
pixi run --locked build/release/bin/FreeCAD
```

FreeCAD and its dependencies retain their own licenses; this patch does not
change the host's licensing. The build uses the same source, lock, preset and
patch as the local qualification. It is not an installation into the system
FreeCAD directory.

## Verification scope

The rebuilt unmodified host reproduced the four curved-shape counterexamples.
The rebuilt patched host passed the five-case direct stream comparison and 26
native Part tests. The versioned extension has also passed a real CalculiX
static workflow, independent reopening of its result, invalidation after a
1e-10 mm length edit, and explicit legacy migration. These are scoped regression
checks, not certification of every CAD shape or solver error bound.

`ci/freecad_static.py --recipe static-identity` exercises the installed extension's
metadata mask, unsupported formats, curved-shape reopen/orientation and legacy
migration Undo/Redo. `--recipe static-fingerprint` retains the schema-1 expression
comparison. Both accept the same runtime paths as `--recipe static-task`.

### Curved-boundary transfer regression

The `static-curved-result` recipe exposed a separate CAD–mesh transfer defect on
a radius-10 mm, height-30 mm cylinder with 5 mm mesh size. The old face-mapping
check compared the discretized cap area with CAD area at a fixed 1e-6 relative
threshold. Its 16-segment quadratic circular boundary has relative area error
approximately 4.9318e-5. The independent formula
`A_n = n R² sin(π/n) (4 − cos(π/n)) / 3` agrees with the measured quadratic mesh
area to about 6e-19 m². This is a discretization error, not a persistence error.

The candidate now exports one CAD BREP witness per surface in the same Gmsh
process as the mesh. The request records this export mode; the mesh result
records each witness hash. Replay requires the complete surface set and matching
bytes. Old mesh request/result payloads retain their original fields and hashes.
The export temporarily selects automatic format detection: Gmsh's global
`msh2` option otherwise also affects `Save` commands targeting BREP files.

After the existing trimmed-face node-membership test, numerical CAD intersection
coverage and CAD area/centroid checks use the existing 1e-6 coverage budget.
Successive CAD subtraction also bounds the remaining uncovered area: adding
intersection areas alone can accept duplicated partial witnesses with holes in
their union. Repeated surface identifiers are rejected.
Discretized mesh area and its relative CAD-area difference are reported
separately. This numerical transfer check is not the source-identity contract;
it does not establish exact equality of arbitrarily close CAD shapes.

On the qualified patched host, the cylinder now solves with 1,680 nodes and
912 quadratic tetrahedra. Maximum displacement error against the uniaxial
zero-Poisson reference is 1.639e-11 m. Its result remains current after render,
save and reopen; a height change of 1e-10 mm makes it stale. Missing, modified,
swapped, extra witnesses and a missing witness payload are rejected on replay.
The executable `ci/verify_cad_coverage.py` additionally checks a rectangular
face and a two-face partition under identity and composed X/Z rotations with
translation. It rejects missing halves, duplicated halves, overlap, displaced
planes and disjoint faces (16 cases). This is a bounded OCCT transfer test, not
an end-to-end rotated FreeCAD solve or qualification of arbitrary partitions.

The parametric cylinder also exposed a geometry-axis coefficient changing from
`0` to `-0` on reopen. The current fingerprint normalizes these literal numeric
zero tokens in geometry data, while preserving topology-link orientation signs.

### Rotated parametric cylinder: placement persistence correction

The native `static-curved-result --rotated` qualification composes X=23° and
Z=37° rotations and translates the cylinder by (123, -41, 78) mm. The actual
mesh/CalculiX solve agrees with the rotated uniaxial displacement field within
3.463e-11 m. The initial save/reopen failed despite identical normalized local
BREPs: four matrix entries changed by a few floating-point units.

FreeCAD 1.1.3 `PropertyPlacement::Save` inherited fixed, 16-decimal-place
formatting. `PropertyPlacement::Restore` also reconstructed the quaternion from
axis-angle fields instead of retaining the quaternion saved in the file.
Correcting restoration alone reproduced the failure because writing had already
lost component precision.

The additional `freecad-1.1.3-placement-roundtrip.patch` scopes significant-digit
formatting with `max_digits10` to placement serialization and restores the
original stream flags and precision afterward. Reading retains all four stored
quaternion components without renormalization, alongside the original raw axis
and angle, including full turns. Files without all quaternion fields retain the
existing reconstruction path. No tolerance was added to Vinkulum's identity
check. Lost digits in previously written documents cannot be recovered.

With both host patches, the rotated cylinder passes calculation, rendering,
save/reopen and invalidation after a 1e-10 mm height edit. The separate
`static-pose` recipe checks 12 rotations inside a placed parent over three
save/reopen cycles: local quaternion, raw axis/angle, translation and composed
world matrix all compare exactly (36 checks). Cases include tiny angles,
180°, full turns, negative angles and composed rotations. The host's 21 C++
Rotation/Placement tests also pass. These tests do not establish persistence
for every FreeCAD property, legacy malformed document or CAD operation.

These are local host patches, not fixes available in the standard FreeCAD
1.1.3 distribution. Further packaging and compatibility qualification are required
before publishing this Vinkulum candidate.

Retained reports, captures and hashes: [development qualification](bancs/freecad-geometry-identity-2026/README.md).
