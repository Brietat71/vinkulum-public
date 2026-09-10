# Studio 0.6.0.dev3 — installed-package tetrahedral qualification

The [manifest](qualification.json) identifies the exact installed wheel and
every packaged Python source file, compared byte for byte with the source tree.
Checks ran outside the repository on Linux x86-64. No standalone binary is
qualified by this record.

| Check | Result | Evidence |
|---|---|---|
| Full Studio desktop suite, OCCT 8 and external CalculiX/Pinocchio | 105 executed, 95 passed, 10 Pinocchio numerical-engine cases skipped in the GUI interpreter | [Log](studio-tests.log) |
| Six C3D10 elements, 27 nodes, analytic pure bending | Relative strain-energy error 1.025 × 10⁻⁷ with CalculiX 2.21 | [Recipe](recipe.json), [raw calculation](calculation/result.json) |
| Displayed mesh, stress tables and source/results separation | Real Qt/VTK tests passed | [Application capture](../../assets/studio-tetra-bending.png) |
| Public Studio 0.5.0 C3D8 archive | Reopened by the new installed wheel, original schema preserved | Archive hash and check recorded in the manifest |
| New public bending archive, fresh process with empty executable search path | 1 passed; no GUI/CAD/solver import | [Saved-example log](saved-example-tests.log) |

The eight tetrahedral test methods include independent boundary-integrated
traction patches on linear, quadratic and curved elements; pure-bending
quadratic exactness; linear-element refinement; unequal physical quadrature
weights; schema/edge-identity rejection; exact Bernstein bounds and bounded
subdivision; and an inversion missed by every integration point. The full
desktop run also includes a real curved C3D10 display/save/table transaction.
The saved public archive test was added after the full pass and checked
separately. No packaged runtime source changed afterward; CI includes the new
test in its subsequent full suite.

An initial bending test assumed less than 20% error after a coarse linear mesh
refinement. Measurements contradicted that assumption. The final test uses
geometrically balanced refinements and verifies decreasing energy error and
the variational ordering against the continuum reference. The measured 20.46%
error with 1,920 C3D4 elements remains explicit in the
[guide](../../STUDIO_TETRAHEDRA.md); no production solver code was modified to
force an accuracy threshold.

Pinocchio's engine is absent from the GUI interpreter. Its existing separate
installed Studio 0.6.0.dev2 / Pinocchio 4.1.0 environment served the controller
and composition checks in this run; this record does not claim a new local run
of its separate numerical suite. Public CI installs the current wheel in both
environments and runs that numerical suite separately.

To reproduce the screenshot and inspectable calculation from an installed
Studio environment with `ccx`, run from the repository root:

```sh
xvfb-run -a -s '-screen 0 1600x1100x24' \
  python ci/studio_tetra_recipe.py /tmp/vinkulum-tetra-demo
```

Use a new output directory. The [captured study](calculation/study.json),
[input deck](calculation/study.inp), [raw output](calculation/study.dat) and
[solver log](calculation/solver.log) accompany the result. The archive can be
reopened without launching CalculiX. Scientific status remains `NotAssessed`.

The exact curved-element check establishes a positive local Jacobian
determinant for admitted binary64 geometry. It does not certify global mesh
injectivity, CAD approximation or solution error. Automatic CAD meshing and
interactive support/traction assignment remain future integration work.
