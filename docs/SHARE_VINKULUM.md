# Help people discover Vinkulum

Lead with something people can try: a real part, its calculation and the files
behind the displayed values. Vinkulum is an open engineering workbench for CAD,
multibody dynamics and inspectable scientific results.

[Repository](https://github.com/Brietat71/vinkulum-public) ·
[Linux 0.6.0a1 preview](https://github.com/Brietat71/vinkulum-public/releases/tag/studio-v0.6.0a1-linux) ·
[Try the included examples in five minutes](../apps/studio/packaging/EXAMPLES.md) ·
[Contributor projects](CONTRIBUTOR_PROJECTS.md)

## Choose one demonstration

| Audience | Show | Invite this contribution |
|---|---|---|
| CAD developers | [Edit the stock length; regenerate the cut and fillets](STUDIO_CAD_HISTORY.md) | [An independent STEP fixture (#1)](https://github.com/Brietat71/vinkulum-public/issues/1) |
| CAD/CAE and desktop developers | [Pick faces on a meshed plate; inspect its captured calculation](#from-a-cad-part-to-an-inspectable-calculation) | [A keyboard/high-DPI workflow (#2)](https://github.com/Brietat71/vinkulum-public/issues/2) |
| FEM and numerical-analysis researchers | [Six tetrahedra and an independent elasticity solution](#six-tetrahedra-and-an-independent-elasticity-solution) | [A bending convergence study (#4)](https://github.com/Brietat71/vinkulum-public/issues/4) |
| Dynamics and robotics researchers | [Inspect captured Pinocchio operators](#inspect-the-operators-behind-a-mechanism) | [Applied-load derivatives (#6)](https://github.com/Brietat71/vinkulum-public/issues/6) |

All four workspaces and their examples are in the Linux **0.6.0a1 alpha preview**.
CAD regeneration and native Vinkulum simulation use bundled components. New
meshes, CalculiX calculations and Pinocchio operators need their separate engines.
Inspecting the shipped captures does not execute those engines.

Keep the actual screenshot, version and contribution link together when sharing.
The screenshots show workflows; numerical claims require the associated reference
and stated domain. Linux packaging does not qualify other operating systems.

## A short announcement to adapt

Suggested title: **Vinkulum: CAD, open solvers and the evidence behind a result**.

> Vinkulum brings OCCT 8 / build123d CAD, a Rust mechanics kernel and CalculiX
> and Pinocchio workspaces into a native 3D desktop. The Linux alpha includes
> an editable plate, a saved CAD-to-FEM study and a captured two-link mechanism.
> Download it, change a CAD dimension, or inspect the inputs behind a result.
>
> We are looking for contributors who enjoy making engineering tools reliable:
> CAD regression parts, independent physical references, Qt interactions and
> numerical verification. Four scoped starting issues link to runnable examples.
> It is research software under Apache-2.0; each external engine keeps its own
> licence, assumptions and installation requirements.
>
> Source, examples and Linux download: https://github.com/Brietat71/vinkulum-public

Short version:

> A CAD part, a simulation and the files behind the result. Vinkulum is an open
> engineering workbench with OCCT 8 CAD, a Rust mechanics kernel and experimental
> CalculiX/Pinocchio workspaces. Try the Linux alpha and help with CAD parts,
> independent physical references or desktop interactions:
> https://github.com/Brietat71/vinkulum-public

## A demonstration you can reproduce

Open `Examples/parametric-plate.vinkulum.json` from the Linux archive. Select
the plate and use **Create → Edit CAD features…**. Change the stock length
from **120 to 150 mm**, preview the dependent cut and fillets, then apply one
undoable change. The BREP supplies the resulting SI mass and inertia.

For a source installation, the [CAD recipe](../ci/studio_cad_recipe.py) creates
the original 120 × 75 × 20 mm plate with its offset hole and 2 mm fillets.
The [feature-history guide](STUDIO_CAD_HISTORY.md) adds the regeneration recipe.
Use the [actual application capture](assets/studio-cad-history.png) and invite
a redistributable STEP part with independently known mass properties.

## Try a captured finite-element result

Follow steps 2–3 in [Examples/README.md](../apps/studio/packaging/EXAMPLES.md)
to inspect the plate study and its displacement in the 0.6.0a1 Linux preview.
The application keeps the raw solver files beside the result and rechecks them
when it opens. New static calculations require an installed `ccx` executable.

For a smaller analytic reference, the earlier
[eight-element tension capture](https://github.com/Brietat71/vinkulum-public/releases/download/studio-v0.5.0-linux/Vinkulum-CalculiX-tension-example.zip)
also opens in the new preview. It has a 1,000 N axial load, expected stress
100,000 Pa and elastic energy approximately 0.000238095238 J. See its
[reference and assumptions](CALCULIX_INTEGRATION.md#reproduce-a-study).
This patch test does not establish general finite-element accuracy.

## From a CAD part to an inspectable calculation

Show the [actual mesh workspace](assets/studio-cad-mesh.png) alongside the
[calculated displacement](assets/studio-cad-mesh-static.png). The
[captured example](bancs/studio-cad-meshing-060/mesh-and-static-example.zip)
keeps the solid, quadratic tetrahedra, face conditions, solver input and raw
output together. It is also shipped in `Examples/cad-statics` in the Linux alpha.

Suggested introduction for CAD/CAE developers:

> A CAD part, a clamped face, a pressure load — and the files behind the result.
>
> Vinkulum Studio takes an OCCT 8 solid through Gmsh meshing, 3D face selection
> and a CalculiX static study. Supports and loads remain editable physical
> conditions. A saved calculation reopens without executing either solver,
> including the captured inputs behind the displayed values.
>
> The Linux alpha includes a real perforated-plate example. The repository
> provides screenshots, qualification records and a pinned OCCT 8/Gmsh build
> recipe. We welcome CAD/CAE contributors to challenge face selection and add
> independent FEM references. This is experimental linear statics; the plate
> demonstrates a workflow, with no general stress-accuracy or speed claim.
>
> Try it and contribute: https://github.com/Brietat71/vinkulum-public

## Six tetrahedra and an independent elasticity solution

Open `Examples/tetra-bending/result.json` through **Run → Linear statics ·
CalculiX… → Open result…**. Share the [actual bending capture](assets/studio-tetra-bending.png)
with its [independent solution and convergence table](STUDIO_TETRAHEDRA.md).

Suggested research introduction:

> Six quadratic tetrahedra can represent this particular pure-bending solution.
> What should a simulation workbench let you inspect to verify that result?
>
> Vinkulum publishes the analytic field, a real CalculiX calculation and its
> raw files. The same reference exposes large energy errors in coarse linear
> tetrahedral meshes. The source also checks curved tetrahedra with exact
> arithmetic bounds on their local Jacobian determinant. Geometry validity
> and solution accuracy have separate evidence.
>
> We welcome independent FEM references, adversarial meshes and numerical
> analysis contributions. Try the captured example in the Linux alpha:
> https://github.com/Brietat71/vinkulum-public

This is a polynomial patch reference, not a general accuracy or speed comparison
against other solvers. Preserve that distinction when shortening the introduction.

## Inspect the operators behind a mechanism

Open `Examples/pinocchio/result.json` through **Run → Articulated operators ·
Pinocchio… → Open result…**. Show the [actual workspace](assets/studio-pinocchio.png)
and its [independent two-link Lagrange references](PINOCCHIO_OPERATORS.md).
The example opens without the engine; recomputation uses its separate environment.

A useful research contribution is to derive and qualify the configuration
derivative of applied point loads. The existing intrinsic RNEA derivative must
remain distinguished from the derivative of a loaded problem. This workspace
is limited to fixed-base rigid trees and does not integrate Pinocchio trajectories.

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

For a development update, show the [signed boundary-direction example](STUDIO_CAD_MESHING.md#boundary-direction-symbols):
changing the multiplier from +1 to −1 reverses the pressure arrows while the
support axes stay fixed. The screenshots come from the installed **0.6.0a2.dev1
source wheel**; the downloadable Linux alpha remains **0.6.0a1**. A useful
contribution is to reproduce the nine boundary-workspace tests on a different
GPU or display scale and report the platform and failing image. Keep the
directions' meaning explicit: these sampled arrows do not represent individual
nodal forces or a converged stress solution.

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
