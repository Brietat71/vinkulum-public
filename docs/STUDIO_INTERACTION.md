# Studio: one project window

Studio 0.6.0a2.dev7 revises interaction after user testing of the dev5 Mac DMG
identified missing context actions and slow transitions. This is an initial
correction of those paths; it does not establish that the whole UX is finished.

Click an object to select it. Drag to orbit without changing selection. Right-click
an object for properties, CAD features, framing, visibility or document actions;
right-click the background for creation and scene actions. The browser uses the
same object commands. Result views retain their read-only document semantics.
A click on empty space clears selection. A manipulator gesture stays separate
from a selection click.

The project browser stays on the left. Select a body to edit its geometry and
properties; select an entry under **Analyses** to show that analysis and its
parameters in the same main window. **+ Analysis** adds or selects motion,
CAD mesh and conditions, linear statics or articulated analysis. Motion results
appear in the same browser after a calculation. There are no separate workshop
windows or Back buttons in this path.

An analysis view is created on first use and retained while selecting other
items. Its pending inputs, captured results and camera survive navigation. The
CAD-to-statics handoff also stays in the project window. Background completion
never replaces an unrelated analysis being edited. A CAD-to-statics preparation
also refuses to replace a target whose inputs, captured result or running work
changed while validation was in progress. Right-click an analysis to
show or close it; closing retains the existing cancellation and unsaved-input guards.

Save targets the displayed mechanism or analysis. Analysis inputs still use their
own file formats: saving the mechanism JSON does **not** save every open analysis
or its results. Results keep their captured inputs and provenance. Editing the
mechanism does not silently update an analysis capture. Existing standalone CLI
analysis entry points remain available.

F frames the selection and Shift+F frames all visible objects in the active view
while retaining the viewing direction. Motion/model navigation retains the camera.

Mouse wheels support fractional steps. On macOS, continuous two-finger trackpad
scrolling pans and a native pinch zooms. Secondary click or Control-click opens
a context menu. This follows the [Qt wheel-event contract](https://doc.qt.io/qt-6/qwheelevent.html)
and [native gesture events](https://doc.qt.io/qt-6/qnativegestureevent.html).
The optional navigation hint explains these gestures.

## Work removed from the interaction path

Selecting an object previously reconstructed the complete browser and diagnostics.
It now updates selection and the inspector while retaining existing scene actors.
Project diagnostics are cached by immutable document identity: a changed document
is checked again. This cache does not replace solver admission or result checks.

Geometry, editability, selection and pose updates previously each rendered a
frame synchronously. Their paint requests are now coalesced into the final frame
by Qt. Screenshot capture remains synchronous. Native geometry and numerical
algorithms are unchanged. The dev6 worker reuse also remains available for CAD.

## Evidence and remaining work

Real Qt/OpenGL tests exercise selection, context actions, camera preservation,
fractional wheel input and synthetic native pinch events. Mac Control-click and
continuous pixel pan are also exercised through Qt events. Project navigation
tests click browser entries, retain a modified material input and camera, save
the active analysis, preserve a view during background refresh, and check the
CAD-to-statics handoff and refusal to discard modified inputs. Embedded views
must share the editor's top-level window and keep its browser visible.

The frozen application recipe checks the same embedding for saved Pinocchio,
CAD mesh and static captures, and retains full-window screenshots. Native macOS
qualification remains required for each DMG; Linux tests alone do not establish
Mac latency or trackpad ergonomics.

The retained [initial Linux diagnostic](bancs/studio-interaction-dev7/README.md)
is from the earlier interaction increment, before unified navigation. Its
286.2 → 99.2 ms plate-command observation measures synchronous return only.
It is neither a final-source benchmark nor a Mac speedup claim. The current
`ci/studio_interaction_recipe.py` records installed source hashes, platform,
OpenGL driver, cold openings and repeated transitions with explicit boundaries.

Cold view construction and large ordinary project reads still contain GUI-thread
work. Each retained view has its own VTK context; this is a single user window,
not a claim of one underlying renderer. This increment does not establish that
the entire UX or arbitrary engineering models are qualified.
