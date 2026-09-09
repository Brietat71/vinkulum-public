<p align="center">
  <img src="docs/assets/vinkulum-banner.svg" alt="Vinkulum — The open engineering workbench. Design. Simulate. Verify." width="100%">
</p>

**Build the mechanism. Run the physics. Inspect the evidence.**

Vinkulum is an open engineering project bringing **CAD, multibody simulation
and scientific verification** into a shared workflow. It combines a Rust
mechanics kernel, a Python API and a native 3D desktop application — with the
ambition of becoming a home for the open solvers engineers and researchers rely on.

[Get started](#get-started) · [Contribute](docs/CONTRIBUTOR_PROJECTS.md) ·
[Scientific guarantees](docs/CERTIFICATION_NOYAU.md) ·
[Support the project](docs/FUNDING.md) · [Documentation technique en français](README.fr.md)

<p align="center">
  <img src="docs/assets/studio-cad.png" alt="Vinkulum Studio: a perforated, filleted CAD plate, model tree, mass and inertia inspector" width="100%">
</p>

*An actual Studio session: create a plate, subtract a cylinder, fillet the edges,
then obtain mass and inertia from the exact solid. The
[example project](examples/studio/platine-percee.vinkulum.json) and
[reproducible CAD recipe](ci/studio_cad_recipe.py) are included.*

## What you can do today

| Layer | Available in the source tree |
|---|---|
| **Design** | OCCT **8.0.1** and adapted **build123d**: primitives, extrusions, solid booleans, all-edge fillets and single-solid STEP import/export. Exact BREP geometry, SI mass properties and a separate display mesh. |
| **Model** | A Qt/VTK workbench for rigid mechanisms: bodies, joints, loads, numerical properties, 3D manipulation, undo/redo and project files. |
| **Simulate** | The native Rust kernel computes in a separate process. Studio captures the model and settings associated with each run. |
| **Examine** | Animate results, inspect curves and samples, compare captured runs on their native time grids, and export with units and provenance. |
| **Research** | Use the broader Python kernel API for rigid/flexible mechanics, contact and analysis. Explore explicit numerical contracts, independent references and selected Lean proofs. |

**Kernel 0.19.0 · Studio 0.4.0 · Research software under active development.**
Studio currently exposes a subset of the kernel. CAD is an initial solid-modelling
workflow; interactive constrained sketches and a regenerating feature tree are
future work. The complete kernel is not certified. Each guarantee has a stated
domain and its own evidence. See the [CAD contract](docs/STUDIO_CAD.md),
[Studio guide](apps/studio/README.md) and [kernel API](docs/API.md).

## One workbench, several scientific engines

The long-term goal is a common place to prepare models, choose an appropriate
engine and inspect traceable results. Each engine keeps its identity, licence
and physical assumptions.

| Component | Place in the project | Integration status |
|---|---|---|
| **Vinkulum** | General-purpose mechanics and verifiable numerical research | Native kernel; rigid-mechanism Studio adapter available |
| **OCCT 8 + build123d** | Exact CAD and mass properties | Integrated; local compatibility patches and qualification corpus included |
| **Pinocchio** | Articulated-body algorithms, Jacobians and derivatives | [Connector planned](docs/PINOCCHIO_INTEGRATION.md); first target: qualified rigid trees |
| **MBDyn** | Multibody workflows and independent reference calculations | Existing comparison work; Studio connector planned |
| **CalculiX** | Finite-element workflows | Planned |
| **DUST** | Aerodynamic workflows and future coupling | Planned |
| **NeuralFoil** | Airfoil polar workflows | Used by optional validation tooling; Studio workflow planned |

A connector must establish units, frames, supported physics and reproducible
reference cases before it becomes a selectable engine. Future integrations in
this table are not part of the current desktop binary.

## Get started

For development, use **Python 3.14**, Rust/Cargo, a C++17 compiler, a system
linker and [uv](https://docs.astral.sh/uv/). Linux desktop prerequisites and the
standalone packaging recipe are in the [Studio guide](apps/studio/README.md).

```bash
git clone https://github.com/Brietat71/vinkulum-public.git
cd vinkulum-public
uv venv --python 3.14 .venv
source .venv/bin/activate
uv pip install 'maturin>=1.15,<2'
maturin develop --uv --release --extras verification
uv pip install -e './apps/studio[test]'
python -m vinkulum_studio
```

To add CAD, install the reviewed OCCT 8 adaptations (`patch` is required):

```bash
python ci/prepare_cad.py build/cad-sources
uv pip install build/cad-sources/build123d-0.11.1 \
  build/cad-sources/ocpsvg-0.6.0 -e './apps/studio[cad,test]'
python -m vinkulum_studio
```

After a Python/UI edit, restart Studio. Rebuild the native extension only after
Rust changes; build a standalone archive when preparing a delivery.
No Vinkulum package is currently published on PyPI. Check version and platform
in the [public releases](https://github.com/Brietat71/vinkulum-public/releases):
older binaries may predate the features shown in this source tree.

## Evidence you can inspect

Vinkulum publishes the assumptions behind its results: units and frames,
convergence studies, failure cases and bounded guarantees. Passing a test or
agreeing with another solver does not certify every trajectory.

- [Kernel guarantees and remaining obligations](docs/CERTIFICATION_NOYAU.md)
- [Numerical benchmarks and historical comparisons](README.fr.md#résultats-mesurés-et-comparaison-externe)
- [CAD adaptation, analytic checks and upstream test subset](docs/STUDIO_CAD.md)
- [Studio interaction and rendering qualification](docs/STUDIO_GUI_2026.md)
- [Lean proof workspace](preuves/README.md)

Performance contributions start with a reproducible workload and a profile.
Rust is welcome where it improves a measured bottleneck; comparisons must retain
accuracy, account for transfer costs and report regressions as well as gains.

## Help build it

There is room here for **students, researchers, engineers and curious builders**.
A useful contribution can be a failing physical example, a CAD regression part,
a better interaction, a numerical proof or an independently reproduced result.

The [contributor projects](docs/CONTRIBUTOR_PROJECTS.md) describe concrete first
deliverables across dynamics, CAD, Pinocchio, other solver connectors, Rust
performance, Qt and teaching. Read [CONTRIBUTING.md](CONTRIBUTING.md), then
[propose a scoped project](https://github.com/Brietat71/vinkulum-public/issues/new)
or submit a reproducible fix.

If this direction matters to you, **star the repository**, share a real use case,
or help reproduce a benchmark. For labs and organisations interested in funding
maintenance or a public milestone, see [Support Vinkulum](docs/FUNDING.md).

## Licence

Original Vinkulum code is **[Apache-2.0](LICENSE)**. External libraries, solvers
and models retain their own licences and attribution. See
[third-party notices](THIRD_PARTY_NOTICES.md). The CAD compatibility patches are
kept explicit so their provenance can be reviewed and useful changes can return
upstream.
