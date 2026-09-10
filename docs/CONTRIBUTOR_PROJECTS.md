# Build something that advances open engineering

Vinkulum welcomes students, researchers, engineers and people who build for the
pleasure of understanding how things work. A small, reproducible contribution
is a useful starting point. Larger projects should begin with a scoped proposal.

These are proposed projects, not funded positions or scheduled commitments.

Scoped contribution briefs cover [an analytic STEP fixture (#1)](https://github.com/Brietat71/vinkulum-public/issues/1),
[keyboard and high-DPI behaviour (#2)](https://github.com/Brietat71/vinkulum-public/issues/2),
[a CalculiX bending refinement study (#4)](https://github.com/Brietat71/vinkulum-public/issues/4),
and [applied-load derivatives for Pinocchio (#6)](https://github.com/Brietat71/vinkulum-public/issues/6).
Each records the deliverable, acceptance evidence and current status.

| Project | Useful background | A reviewable first result |
|---|---|---|
| CAD regression corpus | CAD, mechanical design | Three redistributable STEP parts with known units, volume, inertia and failure cases |
| Responsive CAD service | Python, processes, profiling | Measure a persistent worker against fresh processes; retain timeout, cancellation and recovery after an OCCT failure |
| Rust mesh/data path | Rust, Python FFI, VTK | Profile a reproducible large model; improve a demonstrated bottleneck with numerical equivalence and memory measurements |
| Mechanism validation | Dynamics, experimental mechanics | An analytic or published benchmark, convergence study and independent reference |
| Verifiable numerics | Numerical analysis, exact arithmetic, Lean | A sharply stated guarantee, executable verifier and adversarial counterexamples |
| [Pinocchio operators](PINOCCHIO_OPERATORS.md) | Rigid-body dynamics, robotics, Python/C++ | Extend independent references for captured operators and [applied-load derivatives](PINOCCHIO_LOADS.md), or compare supported observables with the native kernel |
| MBDyn connector | Multibody simulation | Map a small mechanism, units, frames, solver version and result channels through an explicit adapter |
| [CalculiX workflow](STUDIO_CAD_MESHING.md) | FEM, Python/Qt | Qualify a CAD-derived study under mesh refinement, with an independent reference and explicit locations for reported displacement/stress |
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
