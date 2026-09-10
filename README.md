<p align="center">
  <img src="docs/assets/vinkulum-banner.svg" alt="Vinkulum — open engineering. Design. Simulate. Verify." width="100%">
</p>

**Open mechanics, inside FreeCAD.**

Vinkulum combines a **Rust mechanics kernel**, a Python API and independently
checkable numerical work. **FreeCAD is now the primary desktop interface.**
Keep its modelling tools and feature tree; run Vinkulum in a separate process
and inspect captured motion in the same FreeCAD document.

[Try the FreeCAD extension](apps/freecad/README.md#install) ·
[Contribute](docs/CONTRIBUTOR_PROJECTS.md) ·
[Discuss](https://github.com/Brietat71/vinkulum-public/discussions) ·
[Scientific guarantees](docs/CERTIFICATION_NOYAU.md) ·
[Support the project](docs/FUNDING.md) · [Technical archive in French](README.fr.md)

![The Vinkulum extension in the actual FreeCAD Linux interface](docs/bancs/freecad-extension-010/freecad-extension.png)

*An editable PartDesign pendulum, a native FreeCAD task panel, and a captured
Vinkulum trajectory. The motion uses a temporary copy; the design keeps its
original geometry and placement.*

## Try it in FreeCAD

**FreeCAD extension 0.1.0a2 · Kernel 0.20.0 · Linux first.**

1. [Install the extension and configure its separate engine](apps/freecad/README.md#install).
2. Open **Vinkulum → Open pendulum example** in FreeCAD.
3. Select the body, open **Vinkulum → Motion analysis**, then run and inspect
   the captured native samples. No Vinkulum workbench switch is required.
4. Change the PartDesign pad, capture it again and compare the physical result.

The first extension covers **one rigid solid and one explicit revolute joint**,
with configurable world pivot, axis, density, native step and CPU allocation.
Calculations run outside FreeCAD's GUI process. Cancellation and document close
retire the job; stale geometry is refused for playback. Saving the FreeCAD file
removes the temporary motion shape before serialization. Saved calculations
reopen without running the engine.

Motion analyses live in FreeCAD's document tree. Save their settings with the
design, reopen them by double-click or right-click, and undo input edits using
FreeCAD's normal controls. Each analysis remembers its last calculation folder
for explicit replay; results remain external files.

This is a research alpha. General FreeCAD Assembly conversion, multiple-body
host models and FreeCAD controls for FEM and other engines remain upcoming work.
The [document and lifecycle qualification](docs/bancs/freecad-analysis-010a2/README.md) states
exactly which runtime, physical cases and lifecycle behaviours were exercised.

The qualified Linux host is FreeCAD 1.1.3 / Qt 6. Its own OCCT 7.8.1 stays in
its process. Vinkulum reimports the captured STEP with **OCCT 8.0.1** and checks
volume, centre and the full inertia tensor against FreeCAD before calculating.
The two Python environments keep separate native libraries.

## One project, several scientific engines

The long-term goal is a common engineering environment for open mechanics,
finite elements and other scientific engines. Each engine keeps its identity,
licence, units, physical assumptions and reference cases.

| Component | Available today | FreeCAD interface |
|---|---|---|
| **Vinkulum** | Native Rust mechanics, Python API, captured rigid-body trajectories | First single-solid revolute workflow |
| **OCCT 8 + build123d** | STEP import, solid geometry and SI mass properties | Checked capture boundary with the FreeCAD host |
| **CalculiX** | External linear-elasticity adapter: C3D4, admitted curved C3D10 and affine C3D8; captured results and references | Planned |
| **Gmsh** | External OCCT 8 mesher and captured boundary conditions | Planned |
| **Pinocchio** | Fixed-base rigid-tree operators, derivatives and independent references through a separate worker | Planned |
| **MBDyn** | Existing comparison work | Connector planned |
| **DUST** | Future aerodynamic workflows | Planned |
| **NeuralFoil** | Optional validation tooling | Connector planned |

The broader kernel API includes work beyond the initial FreeCAD host. See the
[kernel API](docs/API.md), [CalculiX contract](docs/CALCULIX_INTEGRATION.md),
[CAD meshing guide](docs/STUDIO_CAD_MESHING.md) and
[Pinocchio operators](docs/PINOCCHIO_OPERATORS.md). These backend capabilities
are not all exposed by the first FreeCAD extension.

## Evidence you can inspect

Vinkulum publishes units, frames, convergence observations, failure cases and
bounded guarantees. Passing a test or matching another solver does not certify
arbitrary trajectories or models.

- [FreeCAD extension: actual GUI, saved files, process lifecycle and frame transport](docs/bancs/freecad-extension-010/README.md)
- [Persistent FreeCAD analyses: save/reopen, Undo/Redo and numeric precision](docs/bancs/freecad-analysis-010a2/README.md)
- [Independent finite-section pendulum reference and four retained FreeCAD captures](docs/bancs/freecad-bridge-2026/README.md)
- [Kernel guarantees and remaining obligations](docs/CERTIFICATION_NOYAU.md)
- [Numerical benchmarks and historical comparisons](README.fr.md#résultats-mesurés-et-comparaison-externe)
- [OCCT 8 adaptation and independent CAD properties](docs/STUDIO_CAD.md)
- [Quadratic tetrahedra, analytic bending and exact local Jacobian bounds](docs/STUDIO_TETRAHEDRA.md)
- [Lean proof workspace](preuves/README.md)

Performance work starts with a reproducible workload and a profile. Rust and
native multithreading belong where measurements justify them. Comparisons must
retain accuracy and include transfer costs, memory use and regressions.

## Development

The [FreeCAD extension guide](apps/freecad/README.md) explains installation,
engine configuration, packaging and the actual FreeCAD qualification recipe.
The extension uses the existing CAD and mechanical adapter modules in the
`vinkulum_studio` Python package, in a separate process. It never opens Studio's
GUI. The engine environment needs Python 3.14, Vinkulum 0.20 and the
[OCCT 8 CAD dependencies](docs/STUDIO_CAD.md).

Development of the custom Studio GUI is paused. Its code and
[previous desktop releases](https://github.com/Brietat71/vinkulum-public/releases)
remain available for reproducibility; new interface work targets FreeCAD.

## Help build it

Students, researchers, engineers and interested contributors can help with a
failing physical example, a CAD regression, a numerical proof or a reproduced
result. Useful priorities for the FreeCAD direction include:

- Convert a small assembly with explicit units, frames and independent mechanics.
- Exercise cancellation, save/reopen and source edits in real FreeCAD sessions.
- Add a STEP solid with independently known mass properties.
- Connect a supported CalculiX study to native FreeCAD controls.
- Profile the capture/process boundary before proposing a performance change.

Read [CONTRIBUTING.md](CONTRIBUTING.md) and the
[contributor projects](docs/CONTRIBUTOR_PROJECTS.md), then
[propose a scoped contribution](https://github.com/Brietat71/vinkulum-public/issues/new).
Several older desktop tasks refer to the paused Studio GUI; prefer FreeCAD
for new interaction work.

If this direction matters to you, **star the repository**, share a real use case
or help reproduce a benchmark. For labs and organisations interested in funding
maintenance or a public milestone, see [Support Vinkulum](docs/FUNDING.md).

## Licence

Original Vinkulum code is **[Apache-2.0](LICENSE)**. FreeCAD, external libraries,
solvers and models retain their own licences and attribution. See
[third-party notices](THIRD_PARTY_NOTICES.md). Compatibility patches and reference
cases remain explicit so their provenance can be reviewed and useful changes
can return upstream.
