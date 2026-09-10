Vinkulum Studio — current Apple Silicon research preview

This DMG includes the native Vinkulum kernel, Qt/VTK, OCCT 8.0.1 and the
adapted build123d. Studio can edit and regenerate CAD, run native multibody
simulations, and inspect saved CAD/FEM/Pinocchio captures without a Python
installation. The release title, DMG filename and build-info.json identify the
exact Studio version and source commit; build-info.json also lists dependencies.

Right-click opens actions for the pointed object or the scene. Dragging orbits
without changing selection; F frames without resetting the viewing direction.
The project browser keeps geometry, analyses and motion results in one window.
Select an analysis to show its parameters; returning preserves its inputs and
camera. CAD-to-statics preparation stays in the same window. Save targets the
current mechanism or analysis; analysis files remain separate from mechanism JSON.
Mac trackpad scrolling pans and pinching zooms. Scene updates are coalesced,
and compatible CAD operations reuse one supervised worker with a short idle expiry.

The current desktop runs CAD admission and saved-result validation in background
threads. Its shared CPU budget bounds admitted native engine jobs. The packaged
examples cover parametric CAD, boundary conditions, a CalculiX displacement and
stress capture, a tetrahedral bending reference, and articulated operators.

Drag the app into Applications. Copy Examples from the DMG to a writable folder
and start with Examples/README.md. Copies also travel inside the app's Resources.
macOS 14 or later on Apple Silicon is required. This build uses ad-hoc signing;
it has no Apple Developer ID signature or notarisation. See INSTALLATION.txt for
macOS's explicit approval procedure.

The application copied from this DMG is checked on a native ARM64 macOS runner:
normal startup, nonblack 3D rendering, native dynamics, OCCT/STEP operations,
parametric CAD regeneration, and unchanged-file reopening of delivered captures.
The attached reports and screenshots record those checks. Live external Gmsh,
CalculiX and Pinocchio execution is not qualified by this macOS build; those
engines remain external and are not bundled. Saved results open without them.

This is an alpha research application. The complete application and arbitrary
engineering models are not certified. Original code is Apache-2.0; dependency
licences and CAD adaptation patches are included.

Contributor projects: https://github.com/Brietat71/vinkulum-public/blob/main/docs/CONTRIBUTOR_PROJECTS.md
