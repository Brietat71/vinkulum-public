# Installed sketch, CAD and Pinocchio qualification — 10 September 2026

**Studio 0.6.0a2.dev2**, built and installed locally on Linux with Python 3.14.7.
This record qualifies a source wheel. The public standalone Linux download
remains **0.6.0a1**; no replacement archive is implied by these captures.

The [user guide](../../STUDIO_SKETCH.md) explains the admitted line profiles and
affine dimensions. This record includes the actual inputs, real worker outputs
and unmodified Qt window captures.

## Inspect or reproduce the demonstration

- [80 mm sketch](sketch-80.png) and [100 mm sketch](sketch-100.png).
- [Regenerated solid before Apply](extruded-bracket.png).
- [Before project](bracket-80.vinkulum.json) and [after project](bracket-100.vinkulum.json), both schema 4.
- [Raw project/operator archive](sketch-and-operators.zip), including the captured
  Pinocchio request, project, state, result and worker log.
- [Recipe report](recipe.json), [execution log](recipe.log) and
  [qualification manifest](qualification.json).

Open the projects with this source version. To reopen the captured operators
without installing Pinocchio, extract the archive, open **Articulated operators /
Pinocchio… → Open result…** and choose `pinocchio/operators`.

The [executable recipe](../../../ci/studio_sketch_recipe.py) uses real Qt actions:
create a sketch extrusion, edit Width from 80 to 100 mm in the upstream sketch,
preview, apply, undo, redo, save and reopen. It then optionally runs the same
captured body through a separate Pinocchio worker.

From an installed CAD wheel, outside the checkout:

```sh
VINKULUM_STUDIO_PY=/absolute/path/to/studio-env/bin/python
VINKULUM_PIN_PY=/absolute/path/to/pinocchio-env/bin/python
VINKULUM_SOURCE=/absolute/path/to/vinkulum-public
cd /tmp
PY="$VINKULUM_STUDIO_PY" VINKULUM_PINOCCHIO_PYTHON="$VINKULUM_PIN_PY" \
  bash "$VINKULUM_SOURCE/ci/studio.sh"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 xvfb-run -a \
  -s '-screen 0 1600x1100x24' "$VINKULUM_STUDIO_PY" \
  "$VINKULUM_SOURCE/ci/studio_sketch_recipe.py" /fresh/output/directory \
  --pinocchio-python "$VINKULUM_PIN_PY"
```

Both environments must install a Studio version that reads schema 4. Pinocchio
4.1.0 is separate from the OCCT 8.0.1.0 / adapted build123d 0.11.1 environment.
Omit `--pinocchio-python` to run only the sketch/CAD transaction. The output
directory must not already exist. A `PYTHONPATH` overlay is insufficient to
qualify CAD and Pinocchio workers because their environment strips that override.

## Results and independent references

| Check | Result |
|---|---|
| [Installed Studio suite](studio-tests.log) | 149 tests, 139 passed and 10 skipped because Pinocchio is deliberately in a separate environment |
| [Separate Pinocchio references](pinocchio-tests.log) | 11 passed, including the 10 skipped engine cases and one repeated admission case |
| [200% sketch interactions](hidpi-tests.log) | Six passed on Qt/X11 with a 3200×2200 physical display |
| Wheel/source/install identity | All 54 Python modules match byte-for-byte in the source, wheel and both installed environments |
| Dependency compatibility | Both installed environments pass `uv pip check` |
| Real transaction | OCCT regeneration, preview isolation, one document undo entry and schema-4 reopening passed |
| Identity and placement | Point and feature UUIDs, body identity and design origin preserved |
| Cross-engine capture | Pinocchio retained the schema-4 CAD graph and matched the two-cuboid mechanical reference |

The bracket is the union of a foot and web with disjoint interiors. At a
10 mm height and density 7800 kg/m³, its mass is
`7800 × (Width × 20 + 30 × 40) × 10 × 10⁻⁹ kg`, with lengths in mm.
The CAD tests independently sum the two cuboids' centres of mass and full
inertia tensors, including the parallel-axis terms. They also check upstream
edits with a rotated body frame, invalid profiles and schema downgrade refusal.

For the final 100 mm bracket, the Pinocchio recipe places a Y revolute joint at
the design origin and evaluates `q = 0.2 rad` with explicitly prescribed gravity
`(0, 0, −9.81) m/s²`. Independent cuboid integrals give:

| Quantity | Reference |
|---|---|
| Mass | 0.2496 kg |
| Scalar joint mass matrix | 0.0005564 kg·m² |
| Gravity effort | −0.09092371137244647 N·m |
| World COM | (0.03713330171187109, 0.02125, −0.002425598683861675) m |

The recipe compares the mass matrix at relative tolerance `10⁻¹¹` and absolute
tolerance `10⁻¹⁵ kg·m²`, gravity effort at `10⁻¹¹` / `10⁻¹⁴ N·m`, and COM position
at absolute tolerance `10⁻¹³ m`. These are test criteria for this case, not
universal accuracy bounds for arbitrary CAD geometry or dynamics.

The pure sketch tests compare graph rank and nearest-seed projection with an
independent dense NumPy reference. They exercise exact decimal cycles and a
`10⁻²⁰ mm` contradiction, entity ordering, bounded literals, crossing/touching
loops, adjacent overlap and near-collinear profiles. GUI tests cover dimension
double-click, drawing/closure, constrained dragging, undo, deletion dependencies,
retained invalid text, solved coordinates on import and repair of out-of-range
imported drafts. Qt callback exceptions are checked explicitly.

The current scope is planar simple line loops and affine X/Y equations. It does
not establish nonlinear constraint solving, general sketching, persistent CAD
face naming, a certified CAD kernel or parity with a mature commercial suite.
