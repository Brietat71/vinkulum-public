# Change a CAD dimension and inspect the regenerated part

**Included in Studio 0.6.0a1; introduced in source 0.6.0.dev2.** Solid features
can be edited after creation. The
[Linux alpha](https://github.com/Brietat71/vinkulum-public/releases/tag/studio-v0.6.0a1-linux)
includes OCCT 8 / build123d and the parametric plate. For a source installation,
use the [qualified CAD environment](STUDIO_CAD.md#reproducible-installation).

## A complete editing transaction

1. Create a CAD part, or open the [parametric plate example](../examples/studio/platine-parametrique.vinkulum.json).
2. Select its body and choose **Edit CAD features…** in the Inspector or
   **Create → Edit CAD features…** (**Ctrl+Shift+H**).
3. Select a feature. Edit its name, dimensions in millimetres, or local placement.
   **Preview** (**Ctrl+Return**) regenerates the dependent solid operations in
   a separate worker. Numerical fields retain full precision when displayed compactly.
4. Inspect **Original part** or **Last preview**. **Show assembly** restores
   the other bodies and attachments for context. The model is unchanged at this point.
5. **Apply to model** commits one undoable change. **Cancel** closes the editor
   without applying its preview. New edits disable Apply until a matching preview succeeds.

![Real CAD editor: stock length changed from 120 to 150 mm, downstream cut and fillet regenerated](assets/studio-cad-history.png)

The demonstration changes the stock length from 120 to 150 mm while retaining
an offset cylindrical cut and a 2 mm all-edge fillet. Its mass changes from
approximately 1.326705377 kg to 1.676901908 kg at 7,800 kg/m³. These are actual
OCCT integral properties, independently of the display tessellation.

```sh
# In the Studio CAD environment, from the repository root:
xvfb-run -a -s '-screen 0 1600x1100x24' \
  python ci/studio_cad_history_recipe.py /tmp/vinkulum-cad-history-demo
```

The recipe constructs each solid through the real CAD worker and saves before/
after projects, the application capture and a JSON report. Use a new output folder.

## What a feature refers to

Every feature has a stable UUID, explicit input feature identities and one
expected solid output. Dependencies determine evaluation order; the array order
and displayed names are not reference keys. Duplicate identities, missing inputs,
cycles, disconnected features and invalid parameters are rejected before geometry
execution. Changing a feature's parameters retains its identity.

Supported features are boxes, cylinders, spheres, XY rectangle/disc extrusions,
captured solids, subtraction, union, intersection and all-edge fillets. The
feature editor changes parameters and names; operations are added through the
existing CAD design dialog. It does not yet offer drag-to-reorder operations,
feature suppression or arbitrary dependency rewiring.

A boolean tool is a **captured solid input**, with its placement in the target
part's design frame, source body UUID and source body fingerprint. Its geometry
remains available if the original tool body is later moved or deleted. Editing
that body does not silently change the cut. The captured tool's placement is
editable in the feature graph; changing its shape currently requires a new
operation. The original tool remains a separate mechanical body until removed
from the mechanism, as in the previous CAD workflow.

Older schema-2 parts and imported STEP solids begin with a captured BREP feature.
Their construction history cannot be recovered from a filename or a prose journal.
They can acquire new dependent operations, but a captured solid does not expose
its unavailable original sketch dimensions.

These are **whole-solid feature references**. Face/edge selection, constrained
sketches, topology split/merge resolution and TNaming persistence are still
unimplemented. This milestone does not satisfy the complete CAD-02–CAD-06
topological-reference requirements of the engineering specification.

## Geometry, mass and attachment frames

Feature geometry uses millimetres in a persistent design frame. The mechanical
body frame is at the current centre of mass, with SI position, mass and inertia.
The recipe stores the design origin in that centred body frame. Regeneration
preserves the design frame's world placement even when a dimension moves the
centre of mass. For an extrusion, changing height keeps the starting profile
plane fixed rather than translating the entire part to retain its old mass centre.

Applying a preview recomputes **homogeneous** mass and inertia using the entered
density, including when the previous body used an explicitly specified inertia.
Joint attachment positions/orientations and load points are rebased so their
world locations remain unchanged. They are not automatically reassigned to a new
face, edge or material point. Such an association needs a separate explicit
geometric/physical intent contract.

Schema 3 persists both the BREP/display mesh and the feature graph. Earlier
schemas remain readable, and a graph cannot masquerade as schema 2. Limits remain
explicit: 100 features, 2 MB of captured input BREP per recipe, 4 MB of total BREP
per document including captured inputs, and the existing display-mesh budgets.

## Failure and verification

The [installed-package qualification](bancs/studio-cad-history-060/README.md)
includes test logs, the demonstration report and hashes of the exact packaged
Python sources. It qualifies this Linux source workflow, not a standalone archive.

The CAD worker reports the input-file hash with its result. The parent validates
the body payload, target identity, homogeneous mass/density relation and, for
regeneration, the requested feature graph and design-frame placement. These
checks detect inconsistent output; they do not independently certify an arbitrary
BREP's volume or authenticate a malicious producer.

The controller distinguishes start failure, native crash, geometry-operation
failure, timeout, cancellation and invalid output. It stops the actual worker
on cancellation or close. Failed previews preserve the original model and the
last successful preview; Apply remains disabled after a failed attempt.

Tests cover upstream dimension changes through cut/fillet operations, independent
volume and inertia formulas, rotated design frames, centre-of-mass shifts,
frozen cutters, stable identities under renaming/array reordering, legacy files,
schema downgrade rejection and attachment preservation. Real Qt tests cover
preview versus apply, exact numeric input, invalid fields, worker failure causes
and a single undo/redo transaction through the keyboard shortcut.

A composition test passes a schema-3 CAD rod through the separate Pinocchio
worker and compares its mass matrix, gravity effort and COM position with
analytic pendulum formulas. Gravity is explicitly fixed to 9.81 m/s² in that
fixture; the document's general default remains 9.80665 m/s².
