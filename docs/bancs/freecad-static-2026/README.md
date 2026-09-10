# Native FreeCAD / CalculiX qualification — 2026-09-10

Source-only experiment; see [scope and reproduction](../../FREECAD_STATIC_EXPERIMENT.md).
The retained `record.zip` contains both original and result FCStd documents,
trimmed-face and whole-solid STEP captures, raw Gmsh and CalculiX archives,
preview data, console output, screenshot, source hashes, exact executed sources,
engine runtime inventory and the mapping counterexample report. Verify with
`sha256sum -c SHA256SUMS`.

The run was executed before commit: `provenance.json` intentionally records the
then-current base commit and dirty source status. Each recorded source digest was
checked against the committed source before archiving. Local absolute paths are
execution provenance, not portable installation locations.

Host: FreeCAD 1.1.3 / Python 3.11 / Qt 6 / OCCT 7.8.1, Linux x86_64.
Separate worker: Python 3.14.7 / Vinkulum 0.20 / OCCT 8.0.1.
Mesher: Gmsh 5.0.0-git-91b4154 with OCCT 8.0.1; solver: CalculiX 2.21.
Two engine threads; quadratic tetrahedra.

| Measured result | Original bar | Rotated and translated bar |
| --- | ---: | ---: |
| Nodes / elements | 611 / 260 | 613 / 262 |
| Maximum displacement error [m] | 4.762e-13 | 7.067e-14 |
| Relative total strain-energy error | 5.000e-8 | 5.000e-8 |
| Net reaction vector error [N] | 1.000e-5 | 5.634e-10 |
| Native result saved and reopened | pass | pass |
| Changed source rejected before document mutation | pass | pass |

Independent archive replay passes for both cases: raw solver evidence is
revalidated, deliberately renumbered mesh surfaces preserve support/force arrays
exactly, and ambiguous geometry, a displaced face, and doubled reference area
are all rejected.

Development exposed two host-API/lifecycle defects: `FemMesh.addVolume` requires
a list, and a deleted Qt process must no longer be retained for final cleanup.
Both were corrected before this retained run; the dedicated runner exited zero.
Earlier intermediate reports written before cleanup were not accepted as a
successful end-to-end qualification. No numerical acceptance budget was relaxed.

This affine elastic reference exercises transfers and conservation checks.
It does not certify arbitrary face matching, mesh convergence, nonlinear physics,
or an interactive FEM workflow. Already imported results are snapshots and do not
automatically become invalid when the design is subsequently edited.
