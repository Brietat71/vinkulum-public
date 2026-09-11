# Your first Vinkulum calculation in FreeCAD

The published **0.1.0a2** extension contains the single-solid pendulum workflow.
Development **0.1.0a3.dev6** includes native Assembly motion, a linear-static task and
a preconfigured tension example. Use a source checkout containing this guide for
the development features; they are not in the released 0.1.0a2 ZIP.

## Install the extension and engine

The qualified desktop is **FreeCAD 1.1.3 / Qt 6 on Linux x86-64**. For the
current static workflow, use the host rebuilt with both the BREP stream and
placement persistence patches in the [automated host build recipe](FREECAD_HOST.md).
The standard FreeCAD distribution does not contain these fixes. Without them,
an unchanged static result can become stale after save/reopen. The extension
ZIP does not patch FreeCAD.

From the repository root, create the extension ZIP:

```sh
python3 apps/freecad/package.py /tmp/Vinkulum-FreeCAD.zip
```

Close FreeCAD. Extract the `Vinkulum` folder into its user `Mod` directory,
replacing the previous extension folder if present. The qualified Linux location
is `~/.local/share/FreeCAD/Mod`; with a custom XDG setup it is
`$XDG_DATA_HOME/FreeCAD/Mod`. `Vinkulum/InitGui.py` must sit directly under that
folder. Restart FreeCAD. This packages Python extension files; it does not rebuild
FreeCAD or the Rust kernel.

Prepare a fresh **0.6.1.dev2** adapter environment using the
[Linux installation guide](FREECAD_ENGINE.md). Earlier engines lack the CAD
face-witness contract required by the current static worker.
Use the virtual environment's Python path printed by that installer. FreeCAD keeps
its own Python and OCCT; install no second OCCT binding into its embedded console.
The engine installer does not supply the external Gmsh or CalculiX executables.
For statics, build the pinned OCCT 8-enabled Gmsh with the
[Linux mesher recipe](../ci/build_mesher.py):

```sh
python3 ci/build_mesher.py /tmp/vinkulum-mesher --jobs 4
/tmp/vinkulum-mesher/gmsh-install/bin/gmsh -info
```

The recipe requires Git, CMake, Ninja, a C++ compiler and Linux/X11 development
headers. Install CalculiX separately; the qualified executable is `ccx` 2.21.
The [engine build record](STUDIO_CAD_MESHING.md#engines-and-local-builds) describes
the pinned sources and build dependencies. These executables are also usable by
the FreeCAD extension; the custom Studio interface is not required.

Choose either workflow below. Opening an example creates or opens editable CAD;
no calculation starts until you click its Run button.

## Rigid motion: a native double pendulum

1. Choose **Vinkulum → Open double-pendulum assembly**. The native Assembly and its
   Revolute joints are already present; this example activates FreeCAD's Assembly
   workbench so its native tools are available.
2. Choose **Vinkulum → Motion analysis…** with that Assembly selected.
3. Set density to **2700 kg/m³**, duration to **0.5 s** and native step to **0.005 s**.
4. Choose **Engine Python** and a directory for captured calculations, then click
   **Run motion**. Use the playback controls after it completes.
5. Save the FCStd document to keep the source and analysis inputs. Retain the
   captured calculation folder for reopening the motion without rerunning it.

Motion displays temporary copies. Saving or closing removes those copies and
preserves the editable source. The [Assembly guide and reference](bancs/freecad-assembly-analysis-2026/README.md)
explain the qualified mechanism and current limits.

## Linear elasticity: a bar in tension

1. Choose **Vinkulum → Open static tension example**. A new editable 120 × 20 ×
   15 mm bar is created, with a fixed end and an outward pressure of **2 MPa**.
   The Static object is selected. The current workbench stays active.
2. Choose **Vinkulum → Static analysis…**. The material, mesh size and boundary
   conditions are ready. This reference deliberately uses **Poisson ratio zero**
   and **Young modulus 210 GPa** so its displacement field is known exactly.
3. Open **Configure engines and output folder**. Choose the separate engine Python,
   an **OCCT 8-enabled Gmsh HXT** executable, a **CalculiX** executable and an
   absolute output directory. The task remembers these paths after launch.
4. Click **Run static analysis**. The result is a native FEM displacement object.
   The free-end displacement should be approximately **0.001143 mm**:
   `displacement = pressure × length / Young modulus` for this reference.
5. Save the FCStd file to retain the inputs and native result. Keep the capture
   folder as well if you want the STEP, mesh and raw solver evidence.

The task displays pressure in MPa: positive acts inward, negative pulls outward.
The example therefore contains **−2 MPa**. Selecting its row loads that value;
enter a new value and click **Update selected pressure**. Right-click →
**Edit pressure…** focuses the value field. Undo restores both the stored pressure
and the selected row's displayed value. With the boundary list focused, **Menu**
or **Shift+F10** opens the context menu for the selected condition.
To try your own part,
select a top-level solid and open Static analysis; select faces in the 3D view
before adding supports and pressure. The [static task guide](FREECAD_STATIC_TASK.md)
covers units, executable requirements, persistence and scientific limits.

Changing material or boundary inputs makes a linked result stale. On recomputation,
the old mesh is hidden and the source is restored. A geometry change requires
explicit face reselection; an old face number is not silently accepted on a new
shape. Cancel, task closure and document closure retire the static process group.

## What this example establishes

The [first-run qualification](bancs/freecad-first-run-2026/README.md) exercises the
actual installed menu, repeated opening/closing, saved inputs, a real static solve
and the reference displacement. This is a small linear-elastic reference, not a
proof for arbitrary geometry, nonlinear material, contact or general FEM accuracy.
For another physical case, include an independent reference and a refinement study.
