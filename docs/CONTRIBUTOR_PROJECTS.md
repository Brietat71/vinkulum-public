# Build something that advances open engineering

Vinkulum welcomes students, researchers, engineers and people who build for the
pleasure of understanding how things work. A small, reproducible contribution
is a useful starting point. Larger projects should begin with a scoped proposal.

These are proposed projects, not funded positions or scheduled commitments.

| Project | Useful background | A reviewable first result |
|---|---|---|
| CAD regression corpus | CAD, mechanical design | Three redistributable STEP parts with known units, volume, inertia and failure cases |
| Responsive CAD service | Python, processes, profiling | Measure a persistent worker against fresh processes; retain timeout, cancellation and recovery after an OCCT failure |
| Rust mesh/data path | Rust, Python FFI, VTK | Profile a reproducible large model; improve a demonstrated bottleneck with numerical equivalence and memory measurements |
| Mechanism validation | Dynamics, experimental mechanics | An analytic or published benchmark, convergence study and independent reference |
| Verifiable numerics | Numerical analysis, exact arithmetic, Lean | A sharply stated guarantee, executable verifier and adversarial counterexamples |
| [Pinocchio connector](PINOCCHIO_INTEGRATION.md) | Rigid-body dynamics, robotics, Python/C++ | A pinned Linux environment and a one-link adapter checked against analytic gravity torque and acceleration |
| MBDyn connector | Multibody simulation | Map a small mechanism, units, frames, solver version and result channels through an explicit adapter |
| [CalculiX workflow](CALCULIX_INTEGRATION.md) | FEM, Python/Qt | Extend the first static-study workspace with mesh/support selection or archived-result reopening; qualify mesh provenance before CAD meshing |
| DUST connector | Aerodynamics | A minimal reproducible coupling contract and an uncoupled reference case first |
| NeuralFoil workflow | Aerodynamics, Python | A traceable polar workflow, including model version and limits, connected to Studio |
| Accessible desktop | Qt, UX | One keyboard-complete workflow, tested at normal and high DPI |
| Teach with Vinkulum | Teaching, technical writing | A runnable lab exercise with expected observations and a worked explanation |

Start an [issue](https://github.com/Brietat71/vinkulum-public/issues/new) with:
the problem, a small example, what success means, the reference you will compare
against, and the part you want to take on. Check existing issues before starting.
The [contribution guide](../CONTRIBUTING.md) describes validation and licensing.

For external components, keep patches/forks explicit, preserve authorship and
licences, and document a route for useful changes to return upstream. Every
solver retains its own identity and scientific assumptions. Integration should
make those assumptions easier to inspect.

There is no blanket promise of supervision, publication credit or funding.
Agree substantial research collaborations and paid deliverables separately.
