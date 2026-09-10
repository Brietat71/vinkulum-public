# Build something that advances open engineering

Vinkulum welcomes students, researchers, engineers and people who build for the
pleasure of understanding how things work. A small, reproducible contribution
is a useful starting point. Larger projects should begin with a scoped proposal.

These are proposed projects, not funded positions or scheduled commitments.

New interface work targets **FreeCAD on Linux**. The standalone Studio GUI is
paused. [Browse current issues](https://github.com/Brietat71/vinkulum-public/issues?q=is%3Aissue+is%3Aopen)
before starting; the table below proposes directions rather than claiming an
issue is unassigned or a feature is still missing.

| Project | Useful background | A reviewable first result |
|---|---|---|
| [CAD regression corpus](../apps/studio/tests/fixtures/analytic_step_corpus/README.md) | CAD, mechanical design | Extend the three analytic parts with a redistributable vendor export or trimmed curved surface; include independent volume/inertia references and a failure case |
| [FreeCAD first run](FREECAD_FIRST_RUN.md) | Python, Qt, usability | Reproduce installation on a clean Linux machine; improve one measured setup or keyboard obstacle with a native-session regression |
| [Capture performance](CAD_PROCESS_REUSE.md) | CAD, processes, profiling | Measure a representative FreeCAD model from capture through worker completion; report GUI latency, memory and interruption behaviour |
| Rust mesh/data path | Rust, Python FFI, VTK | Profile a reproducible large model; improve a demonstrated bottleneck with numerical equivalence and memory measurements |
| Mechanism validation | Dynamics, experimental mechanics | An analytic or published benchmark, convergence study and independent reference |
| Verifiable numerics | Numerical analysis, exact arithmetic, Lean | A sharply stated guarantee, executable verifier and adversarial counterexamples |
| [Pinocchio operators](PINOCCHIO_OPERATORS.md) | Rigid-body dynamics, robotics, Python/C++ | Extend independent references for captured operators and [applied-load derivatives](PINOCCHIO_LOADS.md), or compare supported observables with the native kernel |
| MBDyn connector | Multibody simulation | Map a small mechanism, units, frames, solver version and result channels through an explicit adapter |
| [CalculiX workflow](FREECAD_STATIC_TASK.md) | FEM, Python/Qt | Qualify a CAD-derived study under mesh refinement, with an independent reference and explicit locations for reported displacement/stress |
| DUST connector | Aerodynamics | A minimal reproducible coupling contract and an uncoupled reference case first |
| NeuralFoil workflow | Aerodynamics, Python | A traceable polar workflow, including model version and limits, connected through an explicit external-engine contract |
| Accessible FreeCAD tasks | Qt, UX | One keyboard-complete native workflow, tested on Linux at normal and high DPI |
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
