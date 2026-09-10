# Linear static analysis inside FreeCAD

Development **0.1.0a3.dev2** adds **Vinkulum → Static analysis…** to the current
FreeCAD workbench. It uses native document objects, face selection, properties
and task controls. FreeCAD 1.1.3 / Qt 6 on Linux is the qualified host. It is not
yet included in the public 0.1.0a2 extension release; build the development
extension with `python3 apps/freecad/package.py /tmp/Vinkulum-FreeCAD.zip` from
this checkout, close FreeCAD, and install its `Vinkulum` folder in your user Mod
directory. Keep the host and [OCCT 8 engine](FREECAD_ENGINE.md) interpreters separate.

## Workflow

1. Select one top-level solid or PartDesign Body, then open **Static analysis…**.
   A Static object is created in the document. Double-click it to reopen its task.
2. Set Young modulus in MPa, Poisson ratio, mesh size in mm and engine threads.
3. Select a face in the 3D view and click **Fix selected faces**. Ctrl-select can
   collect several faces. This fixes all three translational components.
4. Select the loaded face, enter **Pressure [MPa]** and click **Apply pressure to
   selected faces**. Positive pressure acts inward; negative pressure pulls outward.
   Select a pressure condition in the list and use **Update selected pressure** to
   edit it. Remove conditions with the button or the list's context menu.
5. Configure the separate engine Python, an OCCT 8-enabled Gmsh HXT executable,
   CalculiX and an absolute output directory. These paths are user preferences.
6. Click **Run static analysis**. Mesh generation and CalculiX run outside FreeCAD's
   GUI process. Cancel stops that job and its child processes. Closing the source
   document or task also cancels the job; Motion and statics share one active-job slot.
   A separate guardian watches the FreeCAD process and terminates the static process
   group if the host disappears abruptly.
7. The result is a native FEM mesh and displacement object saved with the FCStd
   document. Its colour uses displacement magnitude; native mesh coordinates and
   displacement vectors are in mm. The source is hidden while the result is shown.
   Save the FCStd to retain both source and result; raw captures remain in the
   configured run directory.

The material and face conditions are stored in the native document, with Undo/Redo.
Document properties retain full precision; merely reopening the task does not
round their values. `YoungPa` and condition `PressurePa` properties use Pa, while
the task displays MPa. A face condition records its source and the geometry hash
at selection time. If geometry changes, remove and explicitly reselect the faces;
the adapter never silently reuses a face number for new topology.

A result is admitted only if the captured geometry and analysis inputs still
match. A subsequent material/boundary/geometry change marks it **Out of date** and
hides it on document recomputation (and immediately for tracked task edits).
The original source visibility is restored. Reverting inputs does not automatically
show an old result. Deleting or manually altering native objects is not a general
tamper-proof audit mechanism; saved results and raw archives remain snapshots.

## Scope and evidence

One valid top-level solid, fixed supports, pressure, isotropic small-strain linear
elasticity, quadratic tetrahedra, at most 32 captured boundary faces and 6,000
imported result nodes. The UI does not add contact, nonlinear physics, stresses,
load cases, assemblies, adaptive convergence, or a general finite-element error
bound. Geometry membership and area/centroid checks are conservative numerical
guards; arbitrary curved-face correspondence is not certified.

The [underlying transfer experiment](FREECAD_STATIC_EXPERIMENT.md) checks the exact
affine solution on a bar and on a rotated/translated copy. The
[native task qualification](bancs/freecad-static-task-2026/README.md) exercises the
actual installed menu, face buttons, native persistence, calculation, stale inputs
and cancellation. Run it in an isolated empty session using:

```sh
python3 ci/freecad_static.py --recipe static-task \
  --freecad /path/to/FreeCAD/AppRun \
  --engine-python /path/to/engine/venv/bin/python \
  --gmsh /path/to/occt8/gmsh --ccx /path/to/ccx \
  --output /path/to/new-qualification-directory
```

The recipe creates and closes documents and exits FreeCAD. Use a dedicated process,
never load it into a working session. It checks the FreeCAD exit status and rejects
Python tracebacks in the native console in addition to checking its report.
