# Static boundary keyboard menu — 2026-09-10

The native static task opened its pressure context menu with the mouse, but a
QtTest Menu key sent to the focused boundary list did not open it. The retained
baseline exits with a failed report before starting the calculation. This is a
native Qt event regression, not a claim about every desktop's physical keyboard.

The list now handles Menu and Shift+F10 explicitly. It scrolls the current item
into view and opens that item's existing context menu. Other keys retain the
standard QListWidget handling; an empty selection opens no menu.

The corrected installed extension passes seventeen grouped checks in actual
FreeCAD 1.1.3 / Qt 6 on Linux. The mouse context-menu check remains, followed by
Menu and Shift+F10 on the selected pressure. Both focus the pressure value through
the existing menu action. The remaining checks cover native Undo, save/reopen,
an actual Gmsh/CalculiX calculation, stale result refusal and process cancellation.
FreeCAD exits zero without Python tracebacks.

The same seventeen checks pass when launched with `QT_SCALE_FACTOR=2`. This
records the requested scale; the recipe did not measure the device pixel ratio
or qualify arbitrary display layouts. The maximum displacement errors in both
accepted runs are below 4.762e-13 m for the existing affine tension reference.
No physical model or numerical tolerance changed.

`record.zip` retains the failing keyboard baseline and both passing runs,
including exact executed sources, source hashes, installed extension ZIPs,
reports, console logs, native documents and raw solver artifacts. Verify it with
`sha256sum -c SHA256SUMS`. Producer hashes were checked against the executed
sources before archival. The baseline macro is a diagnostic variant replacing
the mouse event with a Menu key; the final recipe preserves the mouse check.

A third passing run (`integrated/`) applies the same correction to main
8485065 after the native FEM input adapter merged. All seventeen checks pass,
including the updated bounded wait for the displayed native Close button.
Its executed source hashes and complete calculation record are also retained.

Run `ci/freecad_static.py --recipe static-task` with the executable arguments in
the [static task guide](../../FREECAD_STATIC_TASK.md). The external engine here
uses Python 3.14.7 / OCCT 8.0.1, Gmsh 5.0.0-git-91b4154 and CalculiX 2.21.
This qualifies the context-menu shortcut, not an entirely keyboard-operated CAD
workflow or other platforms. No published release asset is changed.
