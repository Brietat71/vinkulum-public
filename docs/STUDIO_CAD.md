# Studio CAD — OCCT 8 minimum

**FreeCAD users:** start with the [dedicated Linux engine installation](FREECAD_ENGINE.md).
It creates and checks the separate environment used by the FreeCAD extension.
The Studio GUI workflows and release milestones below are historical; new
interface development targets FreeCAD.

Studio 0.4.0 introduced solid modelling in the multibody editor; 0.4.1 translates
the interface into English. The mechanics kernel remains **Vinkulum 0.19.0**.
[Parametric solid features](STUDIO_CAD_HISTORY.md), regeneration previews and an
undoable apply transaction are included in the **0.6.0a1 Linux alpha**.
The current **0.6.0a2.dev2 source version** adds [constrained line sketches](STUDIO_SKETCH.md)
and editable profile extrusions; that increment requires a source installation.

## Available workflow

The **CAD** button and the Create menu open boxes, cylinders, spheres, XY
rectangle/disc and constrained line-profile extrusions, subtraction, union, intersection, all-edge fillets
and STEP import/export. Changes participate in document undo/redo. Boolean
operations retain tool body B; it remains a mechanical body until removed.

Operations accept numerical parameters. New solids retain an immutable feature
graph: edit an upstream dimension or placement, preview dependent operations,
then explicitly apply the result. Imported and older solids start as captured
BREP inputs; their old prose journals are not replayed as construction programs.
Sketches admit horizontal/vertical segments, fixed points and signed X/Y dimensions.
Curves, nonlinear constraints, persistent face/edge selection, TNaming
integration and multi-part STEP assemblies remain future work. The dialogs
are modal; CAD runs in a separate process with cancellation and a 60-second limit.

## Reproducible installation

In the Python 3.14 Studio environment with its native kernel installed:

```sh
python ci/prepare_cad.py build/cad-sources
uv pip install build/cad-sources/build123d-0.11.1 \
  build/cad-sources/ocpsvg-0.6.0 -e './apps/studio[cad,test]'
python -m vinkulum_studio
```

Source distributions are checked by SHA-256 before extraction. The explicit
[compatibility patches](../ci/patches) have distinct local package versions.
Dependency constraints remain enforced and OCCT 7 is rejected.

| Component | Version | Adaptation |
|---|---|---|
| OCCT through cadquery-ocp-novtk | 8.0.1 / 8.0.1.0.0 | Unmodified upstream binary |
| build123d | 0.11.1+vinkulum.occt8 | OCCT 8 collection imports, Bnd_Box bounds, dependencies |
| ocpsvg | 0.6.0+vinkulum.occt8 | Point collection and OCCT 8 dependency |
| ocp_gordon | 0.3.1 | Upstream version compatible with OCCT 8 |

Old TopTools/TColgp/TColStd aliases are replaced with corresponding concrete
`OCP.collections` types. `Bnd_Box.Get()` now returns an unregistered binding type;
minimum and maximum corners provide the same six coordinates. Shared libraries
remain separate and replaceable in the bundle. Sources and licences are
identified in [third-party notices](../THIRD_PARTY_NOTICES.md).

## Geometry and mechanics

- BREP is stored in millimetres. OCCT reads the units declared in STEP.
- The mechanical body frame is recentered at the solid's centre of mass.
  Positions and display mesh coordinates are converted to metres.
- For homogeneous density ρ, `m = ρ V_mm³ × 10⁻⁹` and
  `I_kg·m² = (I_mm⁵ / V_mm³) × 10⁻⁶ × m`, about the centre of mass.
  Mass properties come from the BREP, independently of tessellation.
- A CAD edit rebases anchors, joint frames and load points to preserve their
  world positions and orientations.
- The mechanical worker consumes the captured document, mass and SI tensor.
  It loads neither OCCT nor build123d and does not derive inertia from a mesh.
- Sketch extrusions use schema 4; earlier parametric CAD features use schema 3.
  Schema 2 BREP projects and schema 1 mechanisms remain readable. Projects keep
  their earlier schema when possible; a sketch cannot be written below schema 4
  and a feature graph cannot be written below schema 3.

Import accepts one valid solid, with at most 8 MB of STEP input. Captured data
is limited to 2 MB of BREP and 30,000 vertices/triangles per part, and 4 MB of
BREP, including captured recipe inputs, and 100,000 mesh elements per document.
Feature graphs admit at most 100 nodes and 2 MB of captured input BREP.
Rejection preserves the previous
document. These bounds define the initial integration's scope.

## Qualification

The 0.4.0 qualification baseline passed **48 Studio tests**, including CAD
workflows and Inspector precision. Its versions, logs and recipe measurements
are in [the Studio CAD 0.4.0 record](bancs/studio-cad-040/README.md).
The English interface adds a compatibility check for displayed joint/law labels
and their saved identifiers. The [0.4.1 local record](bancs/studio-041/README.md)
contains the 49-test log and English demo recipe. Run the current suite with
`bash ci/studio.sh` in the installed Studio CAD environment.

Studio tests compare volumes and inertias with independent analytic formulas,
including an [independently written STEP L bracket](../apps/studio/tests/fixtures/independent_step/README.md)
and an [eccentric curved bore and enclosed void](../apps/studio/tests/fixtures/analytic_step_corpus/README.md).
The three-part corpus checks full rotated inertia, metre/millimetre declarations,
lost-cavity detection and rejection of disjoint solids through the real worker.
Tests also check rotations and centres of mass, booleans, fillets, STEP metre/millimetre
conversion, BREP persistence, attachment rebasing, worker failure and the
CAD → mechanics → Qt/VTK rendering workflow.

Studio dev6 implements [supervised process reuse](CAD_PROCESS_REUSE.md) following
the [CAD service experiment](CAD_SERVICE_EXPERIMENT.md). One compatible idle worker
is retained for 15 seconds by default. Each operation receives its own CPU lease,
captured request and validated result. The retained GUI qualification compares
fresh and reused workers through document application and rendering, with
independent mass/inertia references and memory measurements.

**96 targeted upstream tests** pass for build123d 0.11.1's `test_bound_box`,
`test_mass_properties`, `test_location` and `test_build_part`. This qualifies
the tested paths rather than the whole build123d or OCCT API. `uv pip check`
checks the installed dependency graph. Bundle checks run their own CAD worker,
verify a perforated part and STEP round-trip, then run a double pendulum and
check rendering.

Integral properties are numerical and subject to OCCT's geometric tolerances.
Topological checks and analytic references do not certify arbitrary mechanical
trajectories.

Upstream: [OCCT 8.0.1](https://github.com/Open-Cascade-SAS/OCCT/tree/V8_0_1),
[build123d](https://github.com/gumyr/build123d),
[OCP](https://github.com/CadQuery/OCP),
[ocpsvg](https://pypi.org/project/ocpsvg/0.6.0/),
[ocp_gordon](https://pypi.org/project/ocp-gordon/0.3.1/).
