Vinkulum Studio 0.4.2 — an English CAD and multibody workbench.

This patch fixes normal startup with PySide6 6.11.2. Version 0.4.1 referenced
an unavailable QLocale enum; its diagnostic mode bypassed that initialisation.
Startup and bundle checks now enter through the same application/window setup
as a normal launch, and the extracted Linux archive must pass both checks.

The interface, installation guides and contribution guide are now in English.
Saved joint and law identifiers remain compatible with existing projects.
Compact numerical fields preserve unedited values at full precision.

The CAD workflow introduced in 0.4.0 uses OCCT 8.0.1 and adapted build123d
0.11.1: primitives, extrusions, booleans, all-edge fillets and single-solid
STEP exchange. Mass, centre of mass and inertia come from the BREP. Projects
retain exact geometry and a separate display mesh. CAD runs in a separate
process; a failed operation preserves the document. Interactive constrained
sketches and a regenerating feature tree are future work.

Model, Simulate and Inspect workspaces provide a dense Inspector, X/Y/Z fields,
selection/move/rotate tools, a command palette and comparison of captured runs
on their native time grids. The kernel remains Vinkulum 0.19.0.

Standalone builds include Python 3.14, Qt/PySide6, VTK and CAD. Every distributed
package must pass checks on its own extracted executable: CAD operations,
STEP round-trip, native double-pendulum simulation, OpenGL rendering and version
identification. Read the accompanying report for the actual tested platform;
these checks do not establish general scientific certification.

Linux x86-64 on Ubuntu 24.04 / glibc 2.39 / X11 is the current packaging target.
The Apple Silicon recipe is prepared; Linux qualification does not qualify
macOS. Any attached DMG needs its own macOS report and instructions. Planned
signing is ad hoc, without Apple notarisation.

Pinocchio and future solver connectors remain roadmap work and are not bundled
in this release. Several technical research reports remain in French.

Original code is Apache-2.0; third-party licences and adaptation notices are
included. The tag, commit and build provenance identify the distributed source.
