Vinkulum Studio 0.6.0a1 — a CAD part, its loads and the files behind the result.

This Linux alpha preview brings the developing 0.6 workspaces into one
standalone application. Python, the native Vinkulum 0.19.0 kernel, Qt/VTK,
OCCT 8.0.1 and the adapted build123d are included. No Python installation or
compilation is needed to launch Studio, edit CAD or inspect the shipped examples.

Start with Examples/README.md inside the extracted application folder:

- Edit the plate's stock length and regenerate its dependent cut and fillets.
- Open a captured CAD mesh, select boundary faces in 3D, and edit supports,
  inward pressure or a total force. Undo and redo preserve condition identity.
- Inspect the saved quadratic-tetrahedron CalculiX result, including displacement,
  integration-point stress, reactions, energy and the actual numerical inputs.
- Open the two-link Pinocchio capture and inspect its mass matrix and body
  Jacobians with explicit units, frames and derivative scope.

To compute a new CAD mesh, install Gmsh built with OCCT 8 or newer; the repository
provides a pinned source-build recipe. New static calculations require the
external ccx executable (qualified here with CalculiX 2.21). New Pinocchio
operators require its documented separate Python environment. These three
engines are not bundled. Saved captures open without executing them. Native
Vinkulum simulation and the OCCT/build123d CAD worker are bundled.

The current statics domain is linear isotropic elasticity on admitted C3D4,
curved C3D10 and affine C3D8 meshes. Mesh sizes are bounded. Face identities
belong to a captured mesh; successful remeshing clears previous conditions.
The plate is a workflow demonstration, not a certified stress solution. The
separate tetrahedral pure-bending example has an independent analytic reference.
Local Jacobian bounds, numerical input transport and balance checks have
explicit domains. The complete application and arbitrary models are not certified.

Packaging targets Linux x86-64, Ubuntu 24.04 / glibc 2.39, X11 and Mesa OpenGL.
The system supplies its graphics driver and Qt/X11 platform libraries. Older
glibc, ARM64 and native Wayland remain unqualified. See INSTALLATION.txt.
The extracted archive must pass startup, rendering, CAD/STEP, feature
regeneration, native dynamics, external Gmsh/OCCT 8 and CalculiX execution,
and unchanged-file reopening of the delivered CAD/FEM/Pinocchio captures.
The accompanying build-info.json and check reports identify the actual build.
This Linux release does not qualify macOS or provide a new Apple Silicon DMG.

Contributors: try an example and help with an independent CAD part, a desktop
interaction check, FEM convergence or articulated-load derivatives. Concrete
starting tasks: https://github.com/Brietat71/vinkulum-public/blob/main/docs/CONTRIBUTOR_PROJECTS.md

Original code is Apache-2.0. Dependency notices, licences and CAD adaptation
patches are included; external engines retain their own licences and assumptions.
