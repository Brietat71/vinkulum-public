# CalculiX integration — first linear statics adapter

**Status: experimental Python/CLI adapter on the `feat/calculix-static` branch.**
It is not yet connected to Studio's graphical workspaces and is not included in
the published Studio 0.4.1 binary. The native multibody kernel remains unchanged.

This first adapter generates a CalculiX input deck from a validated mesh study,
runs an installed `ccx` executable and reads the final displacement, external
force, integration-point stress and energy tables. Every run retains its input,
raw output, solver log and versioned result metadata.

## Reproduce a study

Install Studio from this branch in an editable Python 3.14 environment. Install
CalculiX separately; on Ubuntu 24.04 its package is `calculix-ccx`. The local
qualification used **CalculiX 2.21**, as reported by the actual executable.
Other versions and macOS have not been qualified for this adapter.

```sh
python -m vinkulum_studio.calculix examples/studio/fem/tension.ccx.json \
  --output /tmp/vinkulum-tension-study
```

Choose a new output directory for every run. An existing directory is refused,
so previous results cannot be overwritten. `--ccx /path/to/ccx` selects the
installed executable. `--timeout 60` limits the solve; version identification
has a separate ten-second limit. Failure preserves input and diagnostic files
without writing a successful `result.json`.

The example is a 1 × 0.1 × 0.1 m specimen with eight affine C3D8 elements,
Young modulus 210 GPa, Poisson ratio 0.3 and a total axial load of 1,000 N.
Its sides are traction-free and supports allow Poisson contraction. Independent
uniaxial elasticity gives:

- axial stress: `F / A = 100,000 Pa`;
- end displacement: `F L / (E A) ≈ 4.76190476 × 10⁻⁷ m`;
- elastic energy: `F² L / (2 E A) ≈ 2.38095238 × 10⁻⁴ J`.

This uniform-strain patch case is exactly representable by the element family
in exact arithmetic. It does not establish accuracy for bending, general meshes
or arbitrary CAD parts. The native `.dat` tables round values to seven
significant digits; the reference comparisons allow relative error `5 × 10⁻⁷`
plus explicit scale-dependent absolute tolerances for nominally zero channels.

Python callers can use `load_study(path)` and `run_static(study, new_directory)`
from `vinkulum_studio.calculix` without creating a Qt application.

## Supported contract

The JSON format is `vinkulum-static-study`, schema 1, containing a `study` object:
`nodes`, `elements`, `fixed_dofs`, `forces`, `young_pa` and `poisson`.

| Field | Meaning |
|---|---|
| `nodes` | World X/Y/Z coordinates in metres; array position defines node ID 1..N |
| `elements` | Eight node IDs in CalculiX C3D8 ordering; element IDs are 1..M |
| `fixed_dofs` | `(node_id, axis)` rows; axes 1/2/3 mean world X/Y/Z, prescribed displacement zero |
| `forces` | `(node_id, Fx, Fy, Fz)` rows in newtons; combine loads on the same node first |
| `young_pa`, `poisson` | One homogeneous, isotropic, linear elastic material; E > 0 and −1 < ν < 0.5 |

Only affine hexahedra are accepted. The mesh must be connected through shared
faces, without repeated elements, unused nodes or faces with more than two
owners. Inverted or degenerate elements and unresolved global rigid modes are
rejected. The scaled determinant threshold is `10⁻¹²`; affine geometry is
checked within `10⁻¹⁰` in dimensionless element coordinates, including the thin
direction of an anisotropic element. Support rank is checked
on dimensionless rigid-mode rows with a singular-value tolerance of `10⁻¹⁰`.
These are numerical admission rules, not a proof that an arbitrary mesh is
geometrically valid or well conditioned. Intersections between remote elements
are not detected; nearly incompressible materials can exhibit locking.

Budgets are 6,000 nodes, 5,000 elements, a 4 MiB input JSON document and 32 MiB
per output/log file. Solver and output limits are not OS memory limits or a
security sandbox. The user chooses a trusted installed executable.

This version has no tetrahedra, warped elements, nonlinear materials, contact,
nonzero prescribed displacements, MPCs, distributed-load cards, dynamics or
automatic CAD meshing. Nodal loads can represent consistently integrated face
tractions, as in the example.

## Result meaning and verification

`result.json` retains displacement in metres, forces and reactions in newtons,
stress in pascals, energy density in J/m³ and total strain energy in joules.
Stresses remain at the eight integration points of each element; they are not
silently extrapolated or averaged at nodes. Component order is
`xx, yy, zz, xy, xz, yz` in the world frame.

CalculiX `RF` contains applied loads as well as support reactions. For this
adapter's nodal-load-only contract, reactions are obtained by subtracting the
captured applied forces. The parser checks complete unique identities, finite
values, zero supports, absence of unsupported reactions, resultant force
and moment balance, nonnegative energy density and the linear elastic work identity.
Moment balance uses undeformed coordinates, with centred and scaled lever arms
so its tolerance is independent of the world-frame origin. Force and moment
residuals allow relative rounding of `2 × 10⁻⁶` against the summed force scale.
The work comparison permits `5 × 10⁻⁶` relative mismatch to account for printed
output rounding; it is not a displacement error estimator.

`load_factor: 1.0` denotes the completed static loading step, not physical time.
Metadata records engine/version, executable/input/raw-output SHA-256 hashes,
requested thread limit and adapter schema. Executable identity is checked
before and after execution. Original raw files remain available for inspection.
The scientific status stays **`NotAssessed`**: these checks and patch references
do not certify a general finite-element discretisation error.

## Qualification and next work

Run the dedicated tests with an installed `ccx`:

```sh
python -m unittest discover -s apps/studio/tests -p test_calculix.py -v
```

The tests include tension, compression, zero load, multiple mesh resolutions,
zero and negative Poisson ratios, analytic displacement/stress/energy,
reaction balance, a rotation by axis permutation plus translation, invalid
geometry/material/supports, an inverted thin-element counterexample,
incomplete/corrupted output, an unbalanced couple with zero resultant force,
failed execution, timeout and protection of previous results. The public Linux
CI installs `calculix-ccx` so these reference tests run against a real solver;
local runs without `ccx` explicitly skip the three solver-dependent test methods.

The next integration steps are an asynchronous GUI controller with cancellation,
mesh and support selection, a read-only static-result view, and explicit CAD
mesh provenance. These need their own qualification before a graphical solver
choice is enabled. Refinement studies for bending and validated meshing should
precede a broader element/geometry contract.

CalculiX's [official project](https://www.calculix.de/) and
[upstream documentation](https://www.dhondt.de/) describe its much wider scope.
Original adapter code follows Vinkulum's Apache-2.0 licence. CalculiX retains
its own licence; this change does not bundle its executable.
