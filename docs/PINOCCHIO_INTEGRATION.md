# Pinocchio integration

Status: experimental source adapter with a separate Linux Pinocchio 4.1.0
environment. The [operator CLI and mathematical contract](PINOCCHIO_OPERATORS.md)
are implemented and checked against independent references. The Studio GUI
workflow, trajectory integrator and standalone distribution remain to be built.
Pinocchio is not bundled or selectable in the published Studio 0.5.0 archive.

## Role in Vinkulum

Pinocchio belongs alongside the native Vinkulum kernel as a specialised engine
for articulated rigid bodies. Its upstream implementation offers kinematics,
Jacobians, recursive rigid-body dynamics and analytical derivatives through
C++ and Python interfaces. These are useful building blocks for mechanism
analysis, robotics, optimisation and control. See the
[upstream project](https://github.com/stack-of-tasks/pinocchio) and
[algorithm overview](https://stack-of-tasks.github.io/pinocchio/).

The first integration should expose kinematics, inverse dynamics (RNEA),
generalised mass matrices (CRBA), forward dynamics (ABA) and selected
derivatives. Pinocchio computes dynamics operators; a Studio trajectory also
needs an explicitly chosen time integrator, initialisation and output contract.
Adding the library alone does not deliver a complete simulation workflow.

Vinkulum keeps its general engineering scope. CAD remains with OCCT/build123d;
the native kernel and future MBDyn, CalculiX, DUST and NeuralFoil connectors
retain their own capabilities and assumptions. The expected benefit from
Pinocchio is access to established articulated-body algorithms and independent
comparisons. Any speed claim requires measurements of the complete workflow.

## First supported model contract

Start with fixed-base trees containing rigid bodies and revolute or prismatic
joints, without imposed motion or contact. Define supported loads explicitly;
reject an input when a load or joint law cannot be represented. A disconnected
body, loop, unsupported joint or ambiguous assembly must produce a diagnostic
before execution. Never remove a constraint to make a model fit.

This restriction belongs to the first Vinkulum adapter. It is not a claim that
Pinocchio itself has no constrained-dynamics capabilities. Closed loops, floating
bases and contact each need a later qualification campaign and result semantics.

The adapter must carry:

| Quantity | Required conversion and evidence |
|---|---|
| Body identity | Stable Studio UUID ↔ engine body/joint/frame identifiers |
| Units | SI throughout; the CAD millimetre conversion stays at the geometry boundary |
| Poses | Document body-to-world rotations and translations; explicit parent/child joint transforms and axis signs |
| Inertia | Mass, local centre of mass, and full tensor about the centre of mass; parallel-axis transforms only where required |
| Configuration | Explicit mapping from document poses to joint coordinates; reject incompatible initial poses |
| Velocity | Distinguish configuration dimension `nq` from tangent dimension `nv`; never differentiate quaternion coefficients as angular velocity |
| Spatial vectors | Declare linear/angular ordering and the reference frame for every Jacobian, velocity and wrench |
| Results | Engine and adapter versions, captured inputs, units, frames, integrator, tolerances and native timestamps |

Reuse Studio's existing document and result provenance. Keep engine-specific
objects out of project files and the GUI process. Start with a separate optional
worker environment so a dependency mismatch cannot prevent native Studio/CAD
from starting. A worker is a failure boundary, not a security sandbox.

## Implementation gates

1. **Dependency qualification.** Select a released Pinocchio version, verify
   Python 3.14 and Linux x86-64 availability, pin the complete environment, and
   run upstream examples locally. Audit transitive libraries and licensing.
   Then qualify Apple Silicon separately. No fake availability button in Studio.
2. **Document adapter.** Implement the supported tree contract and rejection
   diagnostics. A one-link pendulum and a two-link arm must retain poses,
   masses, inertia, gravity, joint axes and body identities through conversion.
3. **Numerical qualification.** Compare pendulum gravity torque and acceleration
   with analytic equations; compare a two-link mechanism with independently
   derived Lagrange equations. On varied configurations, check
   `RNEA(q,v,ABA(q,v,tau)) ≈ tau` and `M(q)a + h(q,v) ≈ tau`, within stated
   scale-aware tolerances and with identical supported loads. These consistency
   checks supplement, rather than replace, independent references.
4. **Derivative qualification.** Check Jacobians and selected analytical
   derivatives using tangent-space perturbations, a step-size sweep and rotated
   frames. Specify the reference frame before comparing arrays. Publish failures
   and conditioning alongside errors.
5. **Studio workflow.** Expose only qualified operations, show why a model is
   ineligible, and archive outputs with their engine identity. Compare physical
   observables on captured inputs; retain each trajectory's native time grid.
6. **Performance and distribution.** Measure cold start, model conversion,
   repeated calls, result transfer and full trajectories separately. Use warmups,
   repetitions, fixed thread counts and matching accuracy. Bundle only after
   extracted-executable tests pass on the target platform.

No global kernel certification or performance advantage follows automatically
from agreement between two implementations. Analytic cases, convergence studies
and adversarial frame/unit tests remain necessary.

## Upstream relationship

Pinocchio uses the [BSD-2-Clause licence](https://github.com/stack-of-tasks/pinocchio/blob/devel/LICENSE).
Preserve its notices and attribution, together with each dependency's terms,
when distributing it. The Vinkulum adapter can remain Apache-2.0.

Prefer a thin adapter over an external fork. If a demonstrated compatibility
or performance issue requires a patch, record the upstream revision, retain a
minimal reproducible failure, version the patch and pursue an upstreamable fix.
Rust work should target measured conversion, data-transfer or other bottlenecks;
Pinocchio's existing numerical algorithms are already implemented in C++.

This is a concrete contributor project: the first reviewable deliverable is
the pinned Linux environment, the one-link adapter and its analytic comparison.
See [contributor projects](CONTRIBUTOR_PROJECTS.md).
