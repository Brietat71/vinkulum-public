# Keyboard-only CAD workflow — Linux/X11 qualification

Qualification for [issue #2](https://github.com/Brietat71/vinkulum-public/issues/2),
10 September 2026, based on public commit `bd0c4a0` plus the focus fix described
below. This covers one editing workflow, not general accessibility certification.

## Reproduce

Use the [Studio CAD development environment](STUDIO_CAD.md). On Ubuntu 24.04,
install `xvfb`, `xauth`, `openbox` and `x11-utils` in addition to Studio's Qt/X11
libraries. From the repository, with this checkout installed editable, run:

```sh
PY="$VIRTUAL_ENV/bin/python" bash ci/studio_keyboard.sh /tmp/studio-keyboard
```

For an existing environment without changing its installed Studio package:

```sh
PYTHONPATH="$PWD/apps/studio" PY=/path/to/cad-venv/bin/python \
  bash ci/studio_keyboard.sh /tmp/studio-keyboard
```

The script creates separate X11 desktops at 100% and 200%, starts Openbox,
waits for its window-manager registration and runs real Qt keyboard events.
CAD creation goes through Studio's actual out-of-process OCCT worker. No CAD
result, file dialog or document transaction is mocked. Initial window/focus
setup and dock resizing use Qt APIs; all editing, dialogs, undo/redo and saving
use keys. Screenshots are captured from X11, including the native VTK surface.

Openbox matters: bare Xvfb supplies no window manager and can leave **no active
window** after a modal closes. That effect is separate from the application
focus defect below. The desktop recipe is explicitly enabled with
`VINKULUM_KEYBOARD_TESTS=1` by this script; the ordinary Studio suite keeps its
focused regression and does not require Openbox.

## Exact key sequence

Start a blank project with the viewport focused. Run each sequence with an
inspector requested at 420 and 280 logical pixels. The recipe logs actual window,
dock and device-pixel-ratio values, as well as the Tab counts between targets.
`Ctrl+A`, followed by text, replaces the contents of a focused field.

1. **Alt+G**, then **Escape**: open and cancel CAD. Focus returns to the viewport
   and the document remains empty.
2. **Alt+G**: the Operation field starts on Box. **Tab**, **Shift+Tab**, **Tab**
   checks forward/reverse order and returns to Name. Enter `Keyboard box 137`.
3. **Tab**, enter `20`; **Tab**, enter `30`; **Tab**, enter `40`, for the box
   dimensions in millimetres. **Tab ×5** reaches Create part, through the three
   position fields and density. **Shift+Tab**, **Tab** compares the button with
   and without focus; **Space** creates the actual CAD body.
4. From the restored viewport, **Tab ×6** reaches Mass [kg]. Enter
   `1.3267053771417465`. **Tab** reaches Position X [m]; enter
   `1.234567890123456`. **Tab ×7** reaches Apply properties. Check its visible
   focus with **Shift+Tab**, **Tab**, then press **Space**.
5. Focus resumes at Name. The test changes the dock width to 320 and back and
   verifies that no edit or numeric change occurred. Enter `Renamed box 137`.
   **Tab ×10**, **Shift+Tab**, **Tab**, **Space** applies the rename.
6. **Ctrl+K**, type `Undo`, **Enter** restores the prior body exactly.
   **Ctrl+K**, type `Redo`, **Enter** restores the renamed body exactly.
7. **Ctrl+S**, enter a fresh absolute `.json` path, **Enter** saves. Reopening
   the file through the document reader must reproduce the whole project exactly.
8. **Ctrl+K**, type `Fit selection`, **Enter** frames the moved body for the
   capture. Tab back to Position X to record the visible focus and inspector.

This sequence also checks that numeric/text characters such as `1`, `3`, `7`
and `e` are accepted by the fields rather than consumed as camera/tool shortcuts.
Unchanged CAD data, orientation and other body properties must survive the rename
exactly; the check compares immutable body values, not rounded display strings.

## Observed defect and focused fix

**Reproduction before the fix:** reach Apply properties with Tab and press Space.
The editor commits successfully, but its form refresh destroys the focused button.
The replacement Apply button is disabled; keyboard focus is lost, so the next Tab
no longer has an inspector field as its starting point. This occurred with both
inspector widths on a window-managed X11 desktop.

**Fix:** when Apply properties held focus, restore focus to the first property
after rebuilding the form. For a body this is Name; for project settings it is
Gravity X. Applying from another panel does not take its focus. The existing
validation-error path remains unchanged.

The ordinary workspace regression covers exact mass preservation, continued Tab
navigation, the project/vector case and absence of focus stealing. The complete
CAD recipe independently detects the original failure after a real worker result.

## Evidence and limits

Studio `0.6.0a2.dev1`, kernel `0.19.0`, Python `3.14.7`, Qt/PySide6 `6.11.2`,
VTK `9.7.0`, OCCT `8.0.1.0.0`, adapted build123d `0.11.1+vinkulum.occt8`,
Openbox `3.6.1`, Linux x86-64 with glibc 2.39. Windows are 1280 × 860 logical
pixels, on desktops of 1600 × 1100 and 3200 × 2200 pixels. Both inspector widths
pass at both device pixel ratios. Logs and native captures are in
[the qualification folder](bancs/studio-keyboard-2026/).

The dark-theme inspector labels, SI units and focus outline were visually checked
in the native captures at both scales. Long numeric values scroll horizontally
while being edited and use a compact display when unfocused; the saved values
remain exact. The automated button check requires a rendered difference between
focused and unfocused states, not just a `hasFocus()` flag.

This is automated keyboard qualification on Linux/X11 with Openbox. It does not
qualify macOS, Wayland, a screen reader, every desktop/window manager, all themes
or a human user's comfort. No physics, meshing, CAD feature semantics or native
solver implementation is changed by this contribution.
