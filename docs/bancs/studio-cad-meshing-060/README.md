# Studio 0.6.0.dev4 — captured CAD-to-statics qualification

This is a Linux **installed-source-wheel** qualification, not a new standalone
release or a general engineering accuracy certificate. The downloadable Linux
application remains Studio 0.5.0.

## Reopen the actual example

Download and extract [mesh-and-static-example.zip](mesh-and-static-example.zip).
The archive is about 2.3 MB and includes a hash manifest and two complete captures:

- `mesh/result.json`: open with **Run → Static study from CAD… → Open mesh…**.
  Keep the six companion mesh/input/log files alongside it.
- `calculation/study.json`: open with **Open CAD study…** to recover the material,
  clamped face and top pressure. The study includes the original solid and its
  boundary-condition identities.
- `calculation/result.json`: open with **Run → Linear statics · CalculiX… →
  Open result…**. Keep `study.json`, `study.inp`, `study.dat` and `solver.log`
  alongside the result.

Reopening these files requires Studio 0.6.0.dev4 but does not execute Gmsh,
CalculiX or a CAD engine. To enter the CAD study window, select the included
machined plate body first; opening the saved capture then replaces that input
snapshot explicitly. Recalculation requires the separate qualified engines.

SHA-256 of the ZIP:
`ccd0f11e309d8df1a23cce97f1351e96e472de2dac8dc4068607b3aa0eb3df2b`.

## What was checked

The complete desktop suite ran from outside the repository against an installed
wheel: **125 tests, 115 passed, 10 skipped, exit 0**. The ten skipped tests need
Pinocchio in the main test environment. Desktop composition tests used the
already qualified separate Studio 0.6.0.dev2 / Pinocchio 4.1.0 worker.
Its existing operator protocol was not changed here. See [studio-tests.log](studio-tests.log).

An additional [fresh-process archive test](saved-example-tests.log) checked the
ZIP manifest, mesh/raw-data consistency, physical-condition identities and
CalculiX result reopening. It removed the executable search path, prohibited
subprocess execution, verified no files changed and checked that Qt, VTK, OCCT,
build123d and Gmsh were not imported. **One test passed, exit 0.**

The new tests cover actual mouse face picking, additive selection, camera drags,
filtered selection, true quadratic edge display, material/load persistence,
condition identity, undo/redo, transactional failure, recapture and cancellation
during threaded checks. Independent backend references cover pressure force and
moment balance on closed curved tetrahedra, analytic volume, total-force
assembly, Gmsh/CalculiX node ordering and decimal input transport against a
compiled Fortran reader.

The GUI [recipe](../../../ci/studio_mesh_recipe.py) opened the command from the
normal editor, generated a fresh mesh, added conditions through their dialogs,
saved the study, ran CalculiX and reopened the result. The screenshots are
actual application captures:
[mesh workspace](../../assets/studio-cad-mesh.png),
[calculated displacement](../../assets/studio-cad-mesh-static.png).

## Observed plate calculation

The captured 120 × 75 × 20 mm perforated, filleted plate used a 40 mm target
size and quadratic elements. It produced **4,082 nodes, 2,113 C3D10 elements
and 29 boundary groups**. Face 1 fixes all three displacements at the minimum-X
end. Face 11 carries **100,000 Pa inward pressure** on the top surface.
The material is E = 210 GPa, ν = 0.3.

The captured elastic energy is **0.0022972887863252937 J**. This is an observed
workflow result, not an independently certified reference solution. Mesh volume
differs from CAD by **0.07855%**; that difference does not bound stress error.
The saved source and actual input-deck geometry both passed the local element
admission checks. Decimal transport changed 397 coordinate values by at most
`8.673617379884035e-19 m`, and 13 nodal-force values by at most
`4.336808689942018e-19 N`. Material values did not change.

The [recipe report](recipe.json) records this run. Earlier repetitions produced
slightly different connectivity (4,084 nodes / 2,115 elements). The parallel
mesher is not claimed bitwise deterministic. A finer 12 mm request exceeded
the current mesh budget and was refused; it was not silently coarsened.

The [qualification manifest](qualification.json) records the final wheel,
all 47 package modules, test logs and engine fingerprints. Source revisions for
Gmsh and OCCT are pinned in the [build recipe](../../../ci/build_mesher.py).
A fresh end-to-end rebuild of that script is still in progress; the GUI and
numerical checks above used the completed local Gmsh/OCCT 8 build.
