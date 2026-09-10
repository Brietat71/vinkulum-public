Vinkulum Studio 0.5.0 — CAD, mechanisms and a first CalculiX workspace.

Open Run → Linear statics · CalculiX to load a mesh study, adjust its material
and load multiplier, and calculate with a separately installed ccx executable.
Inspect displacement colours, explicitly amplified deformation, reactions,
integration-point stresses and strain energy. Export the captured values to CSV.

Each calculation keeps its mesh, settings, input deck, raw output and result
metadata. Open result… rechecks an existing calculation without running the
solver or replacing the study being edited. The checks compare hashes, raw
values and equilibrium/energy balances. They establish internal consistency,
not archive authorship or general finite-element accuracy.

The initial CalculiX domain is linear isotropic elasticity on affine C3D8
elements, with zero supports and nodal loads. CalculiX 2.21 is the qualified
Linux executable; it is not bundled. CAD-to-FEM meshing, arbitrary elements,
nonlinear materials and contact are future work. Results remain NotAssessed.

Scene framing now includes visible joint and load symbols. Reference grids do
not influence clipping, and glyph sizes remain stable when a mechanism is
translated far from the world origin. No physical model is changed by framing.

The OCCT 8.0.1 / adapted build123d 0.11.1 CAD workflow provides primitives,
extrusions, booleans, all-edge fillets and single-solid STEP exchange. BREP
geometry supplies SI mass properties; a separate mesh supplies display geometry.
Interactive constrained sketches and a regenerating feature tree are future work.

Model, Simulate and Inspect workspaces provide a compact Inspector, 3D tools,
a command palette and comparison of captured native-kernel runs. The interface
is in English; several technical research reports remain in French. The native
kernel remains Vinkulum 0.19.0. Pinocchio and further engine connectors are planned.

Linux x86-64 / Ubuntu 24.04 / glibc 2.39 / X11 is the packaging target. The
standalone archive includes Python, Qt, VTK and CAD. Its extracted executable
must pass normal startup, rendering, CAD/STEP, native pendulum and external
CalculiX checks. See the attached reports for the tested source and platform.
These bounded checks do not certify the entire kernel or arbitrary models.

The Apple Silicon packaging recipe is prepared; this Linux release does not
qualify macOS and provides no new qualified DMG.

Original code is Apache-2.0. Dependency licences and adaptation notices are
included. Each external solver retains its own licence and physical assumptions.
