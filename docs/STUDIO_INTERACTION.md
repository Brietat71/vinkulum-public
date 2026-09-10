# Studio interaction and workspace transitions

Studio 0.6.0a2.dev7 revises interaction after user testing of the dev5 Mac DMG
identified missing context actions and slow transitions. This is an initial
correction of those paths; it does not establish that the whole UX is finished.

Click an object to select it. Drag to orbit without changing selection. Right-click
an object for properties, CAD features, framing, visibility or document actions;
right-click the background for creation and scene actions. The browser uses the
same object commands. Result views retain their read-only document semantics.
A click on empty space clears selection. A manipulator gesture stays separate
from a selection click.

F frames the selection and Shift+F frames all visible objects while retaining
the viewing direction. Model/Simulate transitions retain the camera and browser
state. The Workspaces button exposes CAD-to-statics, linear statics and articulated
analysis. Back returns to the preceding window without closing the study. Opening
that analysis again retains its inputs, camera and render context. Closing a
window still uses its existing cancellation and unsaved-input rules.

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

Real Qt/OpenGL tests exercise the pointer transaction, context actions, camera
preservation, fractional wheel input and synthetic native pinch events. Mac Control-click and continuous pixel pan are also exercised through Qt events. A separate
workspace test edits an input, returns through the visible Back button and opens
the same study again. Existing mechanics, CAD and archive checks remain required.

The first local Linux/X11 diagnostic measured the parametric plate open command's
synchronous block at 286.9 ms before the changes and 100.1 ms after render
coalescing. This is one diagnostic run, not a performance guarantee. First openings
of the three analysis windows still took about 128–150 ms locally after the
change; they need further work. Those figures do not establish Mac latency or
first-frame latency. The raw recipe records its exact measurement boundary.

Native macOS qualification and an updated frozen DMG remain required before
claiming that this increment resolves the user's Mac experience. Cold module
construction and large ordinary project reads still contain GUI-thread work.

Local installed dev7 validation passed with 189 executed tests and 21 optional
skips (210 discovered), followed by the additional Mac input-mapping regression.
All 63 installed application modules match the source bytes. The original
viewport/manipulator, CAD transaction and stored-result checks remain in the suite.
