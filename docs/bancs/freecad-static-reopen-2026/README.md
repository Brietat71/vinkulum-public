# Static result persistence after reopening — 2026-09-11

Development **0.1.0a3.dev6** fixes a FreeCAD persistence defect in the native
linear-static task. Before this change, FreeCAD could add an identity location
or round a curved BREP coordinate while reopening an unchanged FCStd. The raw
geometry hash then differed, so a saved result was marked stale even though its
source, selected faces and native FEM mesh were still present.

New analyses capture the source BREP, its exact global placement and each selected
face BREP when the first boundary condition is created. The raw hash remains the
normal fast check. Only if that hash differs does the task require OCCT boolean
differences of exactly zero for the captured/current solids and every selected
face. Face, edge, vertex and solid counts must also match, and the placement must
match exactly. Any failed import, malformed snapshot, topology change, face change
or nonzero difference remains **Out of date — recalculate**. The checks operate on
temporary shapes and do not mutate the source.

The installed dev6 package passes nineteen grouped checks in FreeCAD 1.1.3 / Qt 6
on Linux. It runs the actual Gmsh/CalculiX tension case, saves and reopens its
result, recomputes the document, and verifies that the result remains current
before testing deletion, material changes, source geometry changes, cancellation
and late-result refusal. The maximum affine-bar displacement error is below
4.762e-13 m. The native process exits zero without a Python traceback.

An FCStd from the earlier development package, which has no BREP snapshot
properties, is opened using the final dev6 package. `inputs()` raises the explicit
geometry-reselection error and recomputation marks its result stale. Dev6 never
silently upgrades an old result: a legacy document continues to use its raw
geometry hash, and stays stale when reopening changes that representation.

`boolean/` records zero two-way volume differences for a saved/reopened box and
cylinder, while 0.001 mm translations and size changes produce nonzero volumes.
`faces/` similarly records zero two-way areas for restored selected faces and
positive areas for a different face. `placement/` records exact stability of a
top-level rotated and translated source over two reopen cycles. `cost/` records
that the solid fallback costs 11.73–13.18 s for a 262-face, 256-hole plate, so it
is deliberately used only after a raw-hash mismatch. `canonical-copy/` and
`normalization/` retain rejected alternatives: root-location cleanup does not
stabilize a cylinder, and an identity geometric transform changed that cylinder's
volume. They are not runtime paths.

`record.zip` contains the accepted task, legacy-refusal probe, exact dev6 package,
producer hashes, executed sources, documents, raw engine evidence and the stated
diagnostic records. Verify `sha256sum -c SHA256SUMS` before extracting. The
archive describes one FreeCAD/OCCT version and these cases; it does not establish
a canonical BREP representation, arbitrary shape equivalence, general FEM
accuracy or a performance bound for all models.
