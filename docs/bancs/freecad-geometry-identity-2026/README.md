# FreeCAD geometry identity: development qualification

Local Linux qualification of extension 0.1.0a3.dev6 and adapter 0.6.1.dev2.
This is not a released binary. See [the contract and host build recipe](../../FREECAD_GEOMETRY_IDENTITY.md).

- `fresh-engine.zip`: clean OCCT 8 installation without Qt/VTK; 39 identity and legacy migration checks; independent result reopening; rotated cylinder solve and persistence with the installed engine.
- `placement.zip`: two scoped FreeCAD host patches; 12 placements over three round trips; 21 native C++ Rotation/Placement checks; rotated cylinder persistence.
- `curved-witnesses.zip`: earlier development capture establishing CAD witness export, separate discrete area error, tampered witness refusal, and the signed-zero counterexample. Its recorded producer versions predate the final dev6 build.
- `static-task.zip`: final dev6 native UI, actual box solve, lifecycle and stale-input checks (18 checks).
- Coverage reports accompany the executable `ci/verify_cad_coverage.py`: 16 bounded rigid-pose/partition cases and replay of the retained cylinder mapping after adding the uncovered-area guard.

Archive hashes are listed in `sha256.json`; archives retain producer hashes and runtime records. Some captures contain absolute temporary paths as provenance. They are evidence, not relocatable runtime installations.

The qualified host is a rebuilt FreeCAD 1.1.3 with both supplied patches. Stock FreeCAD does not provide these fixes. New source identities use a versioned serialized representation; numerical CAD coverage retains an explicit 1e-6 budget. Neither establishes universal CAD equality or solver error bounds. Compatibility of arbitrary older/malformed placement records remains outside this qualification.
