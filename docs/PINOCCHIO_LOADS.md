# Derivatives of applied world loads

The experimental articulated adapter 0.2.0 adds the configuration derivative
missing in [issue #6](https://github.com/Brietat71/vinkulum-public/issues/6).
It retains the intrinsic RNEA channels and reports the applied-load and loaded
inverse derivatives separately. The domain is the existing fixed-base tree of
one-coordinate revolute/prismatic joints, with world forces and free moments
at points fixed in body coordinates. No follower loads, floating bases, contact,
closed loops or trajectory integration are introduced. Results stay `NotAssessed`.

## Derivation and axes

Let `Jv` and `W` be the body COM Jacobian's linear and angular rows, expressed in
world axes. Let `a = R p` be the rotated body-local application-point offset.
The point Jacobian column is `V_i = Jv_i + W_i × a`. For a world force `F` and
free moment `M`, virtual work gives

```text
tau_external_i = V_i · F + W_i · M
d a / d q_j = W_j × a
d V_i / d q_j = d Jv_i / d q_j
                 + (d W_i / d q_j) × a + W_i × (W_j × a)
D_external[i,j] = (d V_i / d q_j) · F + (d W_i / d q_j) · M
```

Thus both the COM Jacobian derivative and point transport are needed. Holding
`F`, `M`, `p` and load evaluation time fixed does **not** make `V` constant.

For each loaded body, restrict the columns to its ordered root-to-body path.
Every earlier coordinate is an ancestor of every later coordinate on that path.
If `j <= i`, perturbing revolute coordinate `j` rotates the point-Jacobian column
`i` as a world vector, giving `d V_i / d q_j = W_j × V_i`. For a prismatic `j`,
`W_j=0` and the same expression vanishes: the downstream point and joint origin
translate together. Equality of mixed derivatives of the Euclidean point
position supplies the other triangle. For angular columns, only a **strict**
ancestor rotates the later axis:

```text
H_point[:,i,j] = W_min(i,j) × V_max(i,j)
H_angular[:,i,j] = W_j × W_i if j < i, else 0
```

Here `H[k,i,j]` means `d J[k,i] / d q_j`: world component, Jacobian column,
differentiated coordinate. Columns outside that body's ancestor path have zero
entries. The implementation contracts these cross products directly with the
loads. It does not allocate a Hessian for every body or call a library tensor
binding. This is an original assembly using the admitted tree and public
Jacobians. Background on differential kinematics:
[Haviland and Corke (2020)](https://arxiv.org/abs/2010.08696).

The force contribution is symmetric because `-F · r(q)` is a scalar potential.
An arbitrary constant world free moment does not have that property. The moment
contribution is generally nonsymmetric; the implementation never symmetrizes it.
Reversed declared joint edges retain their coordinate signs through the
Jacobians. Reordering unrelated branches cannot introduce a cross-branch term.

## Independent references

The [reference module](../apps/studio/tests/fixtures/applied_loads/reference.py)
constructs original fixtures and differentiates explicit rotation matrices.
It does not obtain its expected derivatives from the implementation's
Jacobian cross-product formula, Pinocchio or another solver.

For `R = Rx(q1) Ry(q2)`, point `r = R [0,0,L]`, and `F = [0,0,f]`:

```text
tau_F = -f L [sin(q1) cos(q2), cos(q1) sin(q2)]
D_F = f L [[-cos(q1) cos(q2), sin(q1) sin(q2)],
           [ sin(q1) sin(q2),-cos(q1) cos(q2)]]
```

For the free moment `M = [0,0,mu]`, the angular Jacobian is
`[ex, Rx(q1) ey]`. Therefore

```text
tau_M = [0, mu sin(q1)]
D_M = [[0, 0], [mu cos(q1), 0]]
```

The second case detects a swapped derivative axis or accidental symmetrization.
Further references use a general body-local point and COM offset, arbitrary
world force/moment components, and `r = Rz(q1) (r0 + q2 ex)` for a
revolute/prismatic pair. The numerical tests cover 32 seeded spatial/mixed states,
single-coordinate shapes, zero loads, multiple loads, time laws, independent
branches, rotated/translated worlds, body-coordinate changes about the COM,
and reversal of both root and internal declared edges.

The qualification script also sweeps central-difference increments from `1e-2`
to `1e-8` in each coordinate's unit while holding velocity, requested acceleration
and load time fixed. This is a secondary check. The report retains each absolute
error and the analytic matrix, norm, rank and condition number. Norms and
conditioning use characteristic coordinate units of 1 rad / 1 m and effort
units of 1 N m / 1 N, making the scaled matrix dimensionless. A null condition
number denotes an infinite condition number, not a well-conditioned result.

## Captured result contract

New files use `schema=2`, `adapter_version=0.2.0`:

| Channel | Meaning |
|---|---|
| `intrinsic_inverse_derivatives` | Existing RNEA derivatives including gravity, excluding applied world loads |
| `external_effort_derivatives` | Derivatives of the summed point-force and free-moment generalized effort |
| `loaded_inverse_derivatives` | Intrinsic minus external, the derivative of `inverse_effort` |

Each channel contains `q`, `velocity` and `acceleration` matrices. Row `i` is
effort `i`; column `j` is the differentiated variable `j`. Entry units are row
effort divided by column variable: N m or N over rad or m, with the appropriate
seconds for velocity/acceleration. For the admitted load laws the external
velocity and acceleration matrices are exactly zero, so those loaded matrices
equal their intrinsic counterparts. The acceleration derivative remains `M(q)`.

The engine-free reader validates all dimensions and finite values, the loaded
subtraction identity, zero velocity/acceleration load derivatives, captured
conventions and the external configuration derivative reconstructed from the
captured world Jacobians/loads and tree. That reconstruction shares the assembly
with the producer and is a consistency check, not an independent scientific
reference or producer authentication. Analytic references supply the independent
evidence. A coherently fabricated archive can still pass consistency checks.

The reader also accepts original schema-1 / adapter-0.1.0 archives with their
unchanged convention strings. Their missing load derivatives stay absent.
Studio adds **External dτ/dq** and **Loaded inverse dτ/dq** tables with full
precision CSV export. Opening an older result shows that these channels are
unavailable and disables their export; it does not substitute zeros or silently
recalculate the archive. Existing intrinsic tables keep their meaning.

## Binding observation and reproducibility

The installed Pinocchio 4.1.0 / EigenPy 3.13.0 binding exposes
`computeJointKinematicHessians` and `getFrameKinematicHessian`. On the qualified
Linux/NumPy stack, the latter returns a C-contiguous array of shape `(6,nv,nv)`
whose direct indexing disagrees with the independently differentiated COM
Jacobian. Flattening it in C order and reshaping in Fortran order recovers the
expected `[component,i,j]` values for the recorded spatial two-revolute case.
The [qualification](../ci/qualify_applied_loads.py) retains the raw array, strides,
expected array and both errors. This observation is limited to the recorded
bindings and fixture; it is not a portable tensor-conversion rule. The production
load derivative uses neither this tensor nor the diagnostic unpacking.

Use the separate [qualified worker environment](PINOCCHIO_OPERATORS.md):

```sh
PYTHONPATH="$PWD/apps/studio" .venv-pinocchio/bin/python -m unittest discover \
  -s apps/studio/tests -p test_applied_loads.py -v
PYTHONPATH="$PWD/apps/studio" .venv-pinocchio/bin/python ci/qualify_applied_loads.py \
  --output /tmp/vinkulum-applied-loads --samples 7
```

Every output directory must be new. The script retains inputs and result archives
for the analytic fixtures and 2-, 16- and 32-coordinate mixed trees, at the current
32-body document limit. It separates warmed operator evaluation, load-derivative
assembly with precomputed Jacobians, model construction, Python startup, engine
imports, JSON encoding, file writing and result reading. A separate full CLI run
includes process startup, provenance hashing, computation and file exchange.
File-write timing does not include `fsync`; intervals overlap in what they
measure and must not be added as disjoint phases. No timing ratio is a CI gate.
BLAS/OpenMP/MKL thread settings, CPU affinity, versions, source hashes and all
timing samples are recorded. These measurements make no general speedup claim.

The [retained qualification](bancs/studio-applied-loads-2026/README.md) includes
the numerical tables, raw result archives and an actual Studio capture.
