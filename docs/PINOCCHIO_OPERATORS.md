# Pinocchio operators on a captured Studio mechanism

**Experimental source workflow.** The first adapter reads a Studio project and
computes operators at one articulated state in a separate Python environment.
It is not yet a Studio analysis window, a trajectory integrator or part of the
published Studio 0.5.0 standalone archive. Results remain `NotAssessed`.

## Install and run

The qualified engine is [Pinocchio 4.1.0](https://github.com/stack-of-tasks/pinocchio/releases/tag/v4.1.0)
on Linux x86-64 / CPython 3.14. Its numerical dependencies are pinned in
[`ci/pinocchio-linux-requirements.txt`](../ci/pinocchio-linux-requirements.txt).
Keep this environment separate from the Studio CAD environment.

From a source checkout, build the native wheel once if it is not already available:

```sh
maturin build --release --locked --interpreter python --out dist
python -m venv .venv-pinocchio
.venv-pinocchio/bin/python -m pip install --only-binary=:all: \
  -r ci/pinocchio-linux-requirements.txt dist/*.whl
.venv-pinocchio/bin/python -m pip install -e ./apps/studio
.venv-pinocchio/bin/python -m pip check
```

The current Studio distribution also installs its Qt/VTK dependencies in this
worker environment; the worker imports none of Qt, VTK, OCCT or the native
Vinkulum kernel. No CAD extra is required. A smaller independent worker package
is distribution work still to do. Pinocchio and its dependencies retain their
upstream licences; no external engine is vendored into the repository.

Save the **Double pendulum** example from Studio, then run:

```sh
.venv-pinocchio/bin/python -m vinkulum_studio.pinocchio_backend mechanism.vinkulum.json \
  --output /tmp/vinkulum-articulated-01
```

Every run requires a new output directory. It keeps the captured project,
requested state and result. An existing result is never overwritten. A failure
can leave the input files for diagnosis but does not produce a successful result.
`--state state.json` optionally supplies any of these fields:

```json
{
  "q": [0.2, -0.3],
  "velocity": [0.1, 0.4],
  "acceleration": [0.0, 0.0],
  "effort": [0.0, 0.0],
  "time_s": 0.0
}
```

The vector order is the `coordinate_order` recorded in the result: breadth-first
from ground, siblings ordered by joint UUID. Reordering the project's arrays
does not change that ordering. Omitted `q` uses the declared initial joint
coordinates; omitted velocities, accelerations and actuator efforts are zero.
`time_s` selects the value of load laws, within the captured project duration.
It is not a simulation time step.

Revolute coordinates are in radians and efforts in N m; prismatic coordinates
are in metres and efforts in N. Velocity and acceleration divide the coordinate
unit by seconds and seconds squared. Matrix entries can have mixed units: entry
`M[i,j]` maps acceleration `j` to effort `i`. Do not label the entire matrix as
kg or kg m² when the mechanism mixes revolute and prismatic joints.

## Mechanical conversion

All bodies must form one tree attached to ground (`None`). Multiple branches
attached to ground are allowed. This first adapter admits revolute and prismatic
joints, with homogeneous or explicit physical inertia, and world forces/moments
at body-local application points. It rejects disconnected bodies, closed loops,
spherical/fixed joints, imposed motion and incompatible anchors/frames. These
are adapter limits, not limits of Pinocchio itself.

Studio body origins are centres of mass. The full local inertia tensor is
passed with the body's placement relative to its incoming joint; Pinocchio
performs the corresponding rotation and parallel-axis transformation.

Let `T_PA` locate an attachment frame in the parent body and `T_CB` locate the
other attachment frame in the child body. The body's placement in its incoming
joint is `inverse(T_CB)`. Its incoming joint placement is the parent's body
placement multiplied by `T_PA`. For a declared A→B revolute joint, the relative
attachment rotation is `Rz(q)`; for a prismatic joint the translation is
`q * ez`. When the tree traversal goes from B to A, the joint axis is reversed
and the declared coordinate keeps its sign. No constraint is removed to make
a mechanism fit the tree representation.

Initial revolute coordinates are principal angles, with no inferred turn count.
`joint_at` creates coincident attachment frames at the authored pose, so its
initial coordinate can be zero even when a body is inclined in world space.
The conversion reconstructs the initial body poses and checks them before use.
The numerical world origin is shifted to the first ground attachment for
kinematic calculations; output positions are restored to the document world.
This avoids subtracting large world lever arms when calculating Jacobians.

## Operator meaning

The intrinsic dynamics are

```text
M(q) a + h(q,v) = tau_intrinsic
```

`h` includes gravity and velocity-dependent terms. With applied world loads,

```text
tau_external = sum(J_point.T F_world + J_angular.T M_world)
inverse_effort = M(q) a + h(q,v) - tau_external
forward_acceleration = ABA(q, v, effort + tau_external)
```

Body Jacobians use world-aligned axes **at the body centre of mass**, in the
row order `[vx, vy, vz, wx, wy, wz]`. Thus `J @ velocity` returns m/s and rad/s.
This is Pinocchio's `LOCAL_WORLD_ALIGNED` convention; it does not move the
velocity's reference point to the world origin. See the
[upstream explanation of Jacobian frames](https://github.com/stack-of-tasks/pinocchio/discussions/1837).
For a body-local point `p`, the linear Jacobian is
`J_point = J_linear - skew(R p) J_angular`.

`intrinsic_inverse_derivatives` contains the partial derivatives of intrinsic
RNEA with respect to configuration, velocity and acceleration. It **excludes the
configuration derivative of applied world loads**. Subtracting `tau_external`
from an inverse effort does not make those intrinsic derivatives derivatives
of the complete loaded problem. `d tau_intrinsic / d a = M`.

Kinetic energy is `0.5 v.T M v`. Potential energy is
`-sum(mass * gravity dot world centre-of-mass position)`; translating the world
origin can change its additive constant.

## Independent evidence

Run the eleven dedicated tests in the separate environment:

```sh
.venv-pinocchio/bin/python -m unittest discover \
  -s apps/studio/tests -p test_pinocchio.py -v
```

The suite checks an inclined single link with a full inertia tensor, a two-link
mechanism at 64 seeded states, and a revolute/prismatic mechanism against
independently derived Lagrange equations. It checks analytic inverse-dynamics
derivatives, body positions/Jacobians, edge reversal, array reordering, rotated
body frames, a rotated/translated world, and point-load virtual work. It also
runs the CLI outside the source directory, verifies that the worker does not
import the GUI/native kernel, and rejects unsupported models and invalid states.

For the two-link planar reference, let link lengths be `l1,l2`, centre-of-mass
offsets `c1,c2`, masses `m1,m2`, and transverse inertias `I1,I2`. Define
`k = m2*l1*c2`. With downward vertical as `q1=0` and `q2` relative to link 1:

```text
M11 = I1 + I2 + m1*c1² + m2*(l1²+c2²) + 2*k*cos(q2)
M12 = M21 = I2 + m2*c2² + k*cos(q2)
M22 = I2 + m2*c2²

h1 = -k*sin(q2)*(2*v1*v2+v2²)
     + g*((m1*c1+m2*l1)*sin(q1)+m2*c2*sin(q1+q2))
h2 = k*sin(q2)*v1² + g*m2*c2*sin(q1+q2)
```

Differentiating these expressions supplies the independent derivative reference
in the tests. Agreement between RNEA, CRBA and ABA alone is not an independent
physical reference.

Internal operator residuals use an absolute tolerance of `1e-12` in each row's
effort unit plus a relative tolerance of `1e-10` against its term scale. The
analytic tests impose their own tighter tolerances; these are explicitly coded
in [`test_pinocchio.py`](../apps/studio/tests/test_pinocchio.py).
A successful consistency check does not certify an arbitrary model's accuracy.

Each result records the project/state hashes, installed package versions,
Python/platform identity, adapter-file hashes and installed Pinocchio native
file hashes. Pinocchio file identities are compared before and after the run.
These fingerprints identify files; they do not authenticate their publisher or
constitute an operating-system execution sandbox.

## Remaining integration work

The GUI needs a supervised optional-worker controller, captured-state editing,
operator/derivative inspection and comparisons with native-kernel observables.
Trajectory simulation additionally needs an explicit integrator and its own
convergence/energy campaign. Closed loops, fixed-joint aggregation, floating
bases, contact and derivatives of applied world loads require further work.
macOS and packaged-worker distribution have not been qualified.
