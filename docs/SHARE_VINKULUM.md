# Help people discover Vinkulum

Vinkulum is building an open engineering workbench for CAD, multibody dynamics
and inspectable scientific results. A useful introduction is a real example
and an invitation to contribute something specific.

Repository: [Brietat71/vinkulum-public](https://github.com/Brietat71/vinkulum-public)

## A demonstration you can reproduce

Download the [Linux x86_64 preview](https://github.com/Brietat71/vinkulum-public/releases/tag/studio-v0.5.0-linux)
or install Studio with CAD using the [guide](../apps/studio/README.md), then open
[the machined plate project](../examples/studio/platine-percee.vinkulum.json).
Inspect its dimensions, mass and inertia, or export the solid as STEP from the
CAD dialog. Load the double-pendulum example separately to try simulation,
animation and result curves.

To reproduce the CAD screenshot from real operations:

```sh
# From the repository root, in the Studio CAD environment, on Linux with Xvfb:
xvfb-run -a -s '-screen 0 1600x1100x24' \
  python ci/studio_cad_recipe.py /tmp/vinkulum-demo
```

The recipe creates a 120 × 75 × 20 mm plate, subtracts an offset cylindrical
hole and applies 2 mm fillets. It saves the project, a screenshot and a JSON
report of geometry versions, mass, volume and timings. The displayed mass comes
from the BREP with a homogeneous density of 7,800 kg/m³.

The [README screenshot](assets/studio-cad.png) is an actual application capture.
Use it with a source link and state the Studio version. A screenshot demonstrates
a workflow; it does not establish general solver accuracy.

## Try a captured finite-element result

The Studio 0.5.0 Linux release also provides
[`Vinkulum-CalculiX-tension-example.zip`](https://github.com/Brietat71/vinkulum-public/releases/download/studio-v0.5.0-linux/Vinkulum-CalculiX-tension-example.zip).
Extract it, open **Run → Linear statics · CalculiX…**, then **Open result…**
and select its `calculation/result.json`. This inspection needs no installed
CalculiX: the archive contains the mesh study, input deck, solver log and raw
output required to recheck the stored values. Install `ccx` to run a new study.

The [actual workspace capture](assets/studio-static.png) shows an eight-element
tension specimen. Its expected axial stress is 100,000 Pa and elastic energy
is approximately 0.000238095238 J. The
[reference and assumptions](CALCULIX_INTEGRATION.md#reproduce-a-study) explain
this limited patch test; it does not establish general FEM accuracy.

To regenerate the demonstration in the source environment with `ccx` installed:

```sh
xvfb-run -a -s '-screen 0 1600x1100x24' \
  python ci/studio_static_recipe.py /tmp/vinkulum-static-demo
```

Invite a concrete contribution alongside the demonstration:
[a CAD fixture (#1)](https://github.com/Brietat71/vinkulum-public/issues/1),
[a desktop interaction check (#2)](https://github.com/Brietat71/vinkulum-public/issues/2),
or [a FEM convergence study (#4)](https://github.com/Brietat71/vinkulum-public/issues/4).

## A short announcement to adapt

Suggested title: **Vinkulum: an open engineering workbench for CAD and inspectable simulation**.

> Vinkulum brings OCCT 8 / build123d CAD, a Rust multibody kernel and a first
> CalculiX workspace into a native 3D desktop. The Linux preview is ready to try;
> a small saved FEM example can be inspected without installing the solver.
>
> We are looking for contributors who enjoy making engineering tools reliable:
> independent physical references, CAD regression parts, Qt interactions and
> numerical verification. The README has real screenshots, runnable examples
> and three scoped starting issues. It is early research software under
> Apache-2.0; each external engine keeps its own licence and assumptions.
>
> Source and download: https://github.com/Brietat71/vinkulum-public

## Short introduction

> Vinkulum is an open engineering workbench combining OCCT 8 CAD, a Rust
> multibody kernel and a Python/Qt 3D desktop. Build a mechanism, run the physics
> and inspect the model and settings behind the result. It is early research
> software, with concrete projects for CAD developers, dynamics researchers,
> numerical analysts and Qt contributors. Explore the source and try a small
> example: https://github.com/Brietat71/vinkulum-public

## Introduction for a technical community

> We are building Vinkulum, an Apache-2.0 engineering project with a native
> Rust mechanics kernel, a Python API and a Qt/VTK desktop workbench.
>
> The current CAD workflow uses OCCT 8 and explicit build123d compatibility
> patches: primitives, booleans, all-edge fillets, single-solid STEP exchange,
> and SI mass and inertia from the BREP. The desktop also edits rigid mechanisms,
> runs the native solver in a separate process and compares captured results.
>
> We welcome redistributable CAD test parts, independent dynamics references,
> help with keyboard/DPI behaviour, and carefully scoped solver adapters.
> Studio 0.5.0 also adds an experimental CalculiX statics workspace with
> captured mesh studies, displacement display and checked result reopening.
> Its first domain is affine C3D8 linear elasticity; it requires an installed
> solver and does not automatically mesh CAD parts. Pinocchio and other engine
> workflows remain integration work.
>
> This is research software. Numerical guarantees have explicit domains; the
> complete kernel is not certified. The source, limits and contributor projects
> are public: https://github.com/Brietat71/vinkulum-public

These are draft introductions for maintainers and community members to adapt.
This document does not mean they have been posted or sent anywhere.

## First technical audiences

Reviewed on 10 September 2026:

- **build123d users and contributors:** its
  [Show and tell category](https://github.com/gumyr/build123d/discussions/categories/show-and-tell)
  is a relevant place for the CAD demonstration and a request for regression
  parts. Upstream [OCCT 8 support is already being worked on](https://github.com/gumyr/build123d/discussions/1439).
  Present Vinkulum's patches as a local adaptation and share the qualification
  evidence; do not imply upstream endorsement or a replacement for that work.
- **Show HN:** introduce the overall project with a downloadable application
  and runnable example when the maintainer can answer questions. Its
  [guidelines](https://news.ycombinator.com/showhn.html) call for something people
  can try, discourage fundraising pages and prohibit soliciting votes. Use a
  project introduction rather than a routine version announcement.

Suggested build123d title: **Vinkulum Studio: OCCT 8 solids and SI mass properties
in an open multibody workbench**.

> We have added a small build123d/OCCT 8 workflow to Vinkulum Studio, a Qt/VTK
> engineering workbench. The included example creates a plate, cuts an offset
> hole and fillets its edges. The resulting BREP supplies the body's SI mass
> and inertia, while a separate mesh is used for display.
>
> Our build123d 0.11.1 and ocpsvg adaptations are explicit, versioned patches.
> We checked 96 targeted upstream tests and added independent volume/inertia,
> frame and STEP-unit checks. This is a limited qualification corpus, not a
> claim of compatibility with the entire API. We know upstream OCCT 8 work is
> underway and would welcome guidance on making useful fixes upstreamable.
>
> We would especially value small redistributable STEP parts with known mass
> properties or a reproducible failure. Source, limits and the demonstration:
> https://github.com/Brietat71/vinkulum-public

## Make an introduction useful

- For CAD users: share a reproducible part and ask for a difficult STEP example
  whose redistribution is permitted.
- For researchers: share one reference case, its assumptions and a question
  about the numerical method or verification.
- For desktop developers: show one complete interaction and link to a scoped
  Qt, accessibility or high-DPI task.
- For students: point to a small [contributor project](CONTRIBUTOR_PROJECTS.md)
  with a result they can run and explain.

Before posting in a community, check its current rules and adapt the message
to its interests. Prefer a substantive demonstration and a specific question.
Do not imply endorsements, production readiness, blanket certification or
unmeasured superiority over other solvers.

Stars help discovery. Reproduced examples, useful issues and returning
contributors tell us whether the project is becoming useful. Financial support
is a separate route described in [Funding](FUNDING.md).
