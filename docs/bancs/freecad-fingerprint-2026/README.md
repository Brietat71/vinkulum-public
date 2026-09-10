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
| 0 | 6 | 0.466 ms | 0.292 ms |
| 16 | 22 | 2.167 ms | 1.270 ms |
| 64 | 70 | 7.130 ms | 4.117 ms |
| 256 | 262 | 27.214 ms | 15.513 ms |

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
ZIP container metadata need not be identical. Producer source hashes matched the
checkout when recorded. Verify `sha256sum -c SHA256SUMS` before extracting. Provenance deliberately
records an uncommitted candidate based on bf414f8b6b464605dccb058f84ceab278bd302227,
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

## Process-observation regression

The initial mandatory pre-push run stopped because a process disappeared after
its `/proc/<pid>/stat` file was opened: Linux returned ESRCH (`ProcessLookupError`)
rather than ENOENT. The test observer now treats exactly those two disappearance
errors as dead processes; other errors still propagate. The native cancellation
fixture uses the same rule. The production guardian is unchanged.

Three supervision tests pass, covering real worker/grandchild termination, worker
exit-code propagation and injected disappearance errors with a permission-error
countercheck. The failed log and final test source/log are retained under `guardian/`.
Both native qualifications were rerun after this test correction.

## Composition with the standard FEM adapter

After merging main 8485065, the shared runner preserves `static`, `static-task`
and `static-fingerprint`. `--native-inputs` is accepted only for `static`; both
incompatible combinations are refused before creating the output directory.

All three actual FreeCAD recipes pass on the composition: four cost/equivalence
cases, seventeen native task checks and two standard FEM input cases. All three
processes exit zero without Python tracebacks. `integration-record.zip` retains
all three runs and their exact producer sources/hashes. The two installed ZIPs
have identical file contents; every extension file except its provenance manifest
also matches the earlier qualified dev5 package. The runtime fingerprint change
and the reopening limitation above are unchanged by this integration.

The composition provenance records the uncommitted merge based on
68eb5805ded2e72426a83723240336701cf961ca. No new official release is published by
this record. The CLI merge resolution adds no new physics or acceptance tolerance.
