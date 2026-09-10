# FreeCAD static fingerprint cost — 2026-09-10

Development **0.1.0a3.dev5** computes the static fingerprint from one fresh geometry
serialization and the inputs validated against it. Dev4 serialized the same BREP
once directly and once while reading the inputs. The public inputs and fingerprint
formats are unchanged. There is no cache across calls, and the geometry/refusal
checks are preserved. Solver execution and physical tolerances are unchanged.

An installed FreeCAD experiment compares the frozen dev4 fingerprint expression
with the new implementation on four valid, one-solid plates. Paired timing order
alternates across nine samples per shape. All hashes agree, material and pressure
changes alter both fingerprints, and both implementations refuse stale face
references after a placement change. Restoring the placement restores the original
fingerprint. Instrumentation observes two geometry serializations previously and
one now per fingerprint call.

| Through-holes | BREP faces | Previous median | Current median |
| ---: | ---: | ---: | ---: |
| 0 | 6 | 0.460 ms | 0.288 ms |
| 16 | 22 | 2.145 ms | 1.259 ms |
| 64 | 70 | 7.109 ms | 4.084 ms |
| 256 | 262 | 27.339 ms | 15.592 ms |

These are local fingerprint timings, excluding geometry creation, rendering,
meshing and solver execution. They are not a general GUI latency bound or a
comparison with another CAD package. The plate is 180 × 180 × 10 mm, with
radius-2 mm holes on a 10 mm grid. Raw samples and BREP files are retained.

The installed static task separately passes seventeen grouped checks, including
fingerprint equality/single serialization, native input/result persistence,
pressure editing, context-menu routing, changed-geometry refusal, cancellation,
real accepted and deliberately stale CalculiX runs. Maximum nodal displacement
error against the affine bar reference is 4.762e-13 m. Both FreeCAD processes exit
zero and their native consoles contain no Python traceback.

`record.zip` retains both exact installed ZIPs, producer hashes and executed
sources, raw timing samples, BREP plates, FCStd inputs/results, STEP captures and
mesh/solver evidence. Every packaged file has identical bytes across the two runs;
ZIP container metadata need not be identical. Producer source hashes match the
checkout. Verify `sha256sum -c SHA256SUMS` before extracting. Provenance deliberately
records an uncommitted candidate based on 24ec62ba9eb55077a1f3a90a99d65615bc432aaf,
not an official extension release.

Runtime: FreeCAD 1.1.3 / Qt 6 / Python 3.11 / host OCCT 7.8.1 on Linux x86-64;
separate Python 3.14.7 / Vinkulum 0.20 / OCCT 8.0.1; Gmsh HXT
5.0.0-git-91b4154 / OCCT 8.0.1 and CalculiX 2.21, two engine threads. The cost
experiment itself launches no external engine. Repeat using
`ci/freecad_static.py --recipe static-fingerprint` with the executable/output
arguments from the [static task guide](../../FREECAD_STATIC_TASK.md), then use
`--recipe static-task` for the solver/lifecycle checks.

The new helper remains synchronous in FreeCAD's document context. The tested
geometry/input states establish equivalence for these cases, not arbitrary custom
Python objects or concurrent unsupported mutation of a document. This change does
not solve other potential costs such as STEP export or native result rendering.

## Separate pre-existing reopening limitation

A supplementary cross-version test found that a saved dev4 result initially reads
`Current capture`, but recomputing its fingerprint after reopening marks it stale.
The same document is refused by **both dev4 and dev5**, with identical current
geometry digests and identical stale-boundary errors. This is not introduced or
fixed by the single-serialization change. The saved mesh and result data exist;
that is weaker than proving that their stored capture remains current after reload.

`reopen-counterexample/` retains the original FCStd, old package, both-version
probes/reports and a BREP diagnostic. Comparing a freshly constructed box with the
restored shape reveals an explicit identity location in the latter serialization.
This demonstrates a representation difference; it is not a complete canonicalization
or migration solution. Copying the shape changes location representation further.
No comparison tolerance or geometry guard was weakened. Future stable hashing must
be tested separately across old documents and all users of the shared CAD bridge.
