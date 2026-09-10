# Studio 0.6.0.dev2 — installed-package CAD qualification

These checks ran locally on Linux x86-64 with an installed Studio wheel outside
the source tree. The CAD/GUI environment uses OCCT 8.0.1 and the reviewed build123d
adaptation; Pinocchio 4.1.0 runs in a separate installed environment. The
[manifest](qualification.json) records the wheel, packaged Python source and
evidence hashes. **No standalone archive is qualified by this record.**

| Check | Result | Evidence |
|---|---|---|
| Full installed Studio desktop suite with CAD and external CalculiX/Pinocchio | 86 passed; 10 Pinocchio numerical-engine tests skipped in the CAD/GUI environment | [Desktop log](studio-tests.log) |
| Pinocchio numerical suite in its separate installed environment | 11 passed, none skipped | [Operator log](pinocchio-tests.log) |
| Actual feature-editor preview of the included plate, 120 → 150 mm | Passed; four features, recomputed cut/fillet and mass properties | [Recipe report](recipe.json), [application capture](../../assets/studio-cad-history.png) |

The desktop suite includes independent solid-volume/inertia references, stable
feature identities and dependencies, frozen cutter geometry, design-frame
preservation under centre-of-mass shifts, schema compatibility, distinct worker
failure causes, preview/apply separation and an actual editor undo/redo transaction.
It also includes the CAD-to-Pinocchio composition case: a schema-3 CAD rod is
passed to the separate worker, and its mass matrix, gravity effort and centre
of mass are compared with analytic pendulum expressions.

The first full pass exposed a mistake in that composition test's reference:
the model retained default gravity 9.80665 m/s², while the expected effort used
9.81 m/s². The fixture now sets gravity explicitly to 9.81 m/s². No production
physics code was changed to make this check pass. The logs here are the final
successful runs. The shared Pinocchio admission test runs in both environments;
the table records separate executions, not 97 distinct test cases.

The [public example](../../../examples/studio/platine-parametrique.vinkulum.json)
was reopened from its copied repository path with the installed package and
regenerated to 150 mm. Its resulting mass was 1.676901907865507 kg, matching the
preview report. Run [the capture recipe](../../../ci/studio_cad_history_recipe.py)
in a fresh output folder to reconstruct the before/after projects and screenshot.

The manifest verifies that every packaged Python file matches both the source
tree and the two installed environments. This identifies the tested source
without retroactively claiming a clean Git build. Documentation, CI wiring and
the explicit-gravity test fixture were finalized after the wheel was built;
the packaged runtime files did not change.

These are whole-solid parametric features. Constrained sketches, persistent
face/edge naming, topological split/merge resolution, CAD-to-FEM meshing and
automatic attachment to evolving material points remain unimplemented. The
worker response checks detect inconsistent payloads; they are not independent
certification of arbitrary BREP properties. See the [workflow and limits](../../STUDIO_CAD_HISTORY.md).

The public [Studio 0.5.0 Linux download](https://github.com/Brietat71/vinkulum-public/releases/tag/studio-v0.5.0-linux)
predates this editor and the Pinocchio workspace.
