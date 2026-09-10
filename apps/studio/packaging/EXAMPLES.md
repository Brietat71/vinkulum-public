# Try Studio in five minutes

The files in this folder travel with the application. Work on a copy if
you want to preserve the original examples. The projects and captures are also
available at https://github.com/Brietat71/vinkulum-public.

1. **Edit a real part.** Open `parametric-plate.vinkulum.json` from the File menu.
   Select the plate and open its CAD history (Ctrl+Shift+H). Change the stock
   length from 120 to 150 mm, preview the dependent cut and fillets, then apply.
   Undo restores the previous solid. This uses the bundled OCCT 8 CAD worker.
2. **Inspect an existing CAD study.** With the plate selected, use **Run →
   Static study from CAD…**. Choose **Open CAD study…** and select
   `cad-statics/calculation/study.json`. The captured 120 mm plate, its mesh,
   clamped end and top pressure replace the window's input snapshot. Click its
   boundary faces; Ctrl-click adds or removes a face from the selection.
3. **See its displacement.** Use **Open in CalculiX**, then **Open result…**
   and select `cad-statics/calculation/result.json`. Deformation is visibly
   amplified; displacement values remain in metres. The raw input and solver
   output stay alongside the result and are rechecked when it opens.
4. **Compare with an analytic case.** In the CalculiX window, open
   `tetra-bending/result.json`. This six-element quadratic tetrahedron case has
   a published pure-bending reference. Its special polynomial solution does
   not establish general finite-element accuracy.
5. **Inspect dynamics operators.** Use **Run → Articulated operators · Pinocchio…**,
   then **Open result…** and select `pinocchio/result.json`. Inspect the mass
   matrix, state and body Jacobians of the captured two-link mechanism.

Steps 2–5 inspect saved data without executing Gmsh, CalculiX or Pinocchio.
To compute new results, install the relevant external engine: Gmsh built with
OCCT 8 or newer for CAD meshing, `ccx` for statics, or the documented separate
Pinocchio environment for articulated operators. These engines are not included
in this archive. Native Vinkulum simulations and CAD regeneration are included.

This is an alpha research preview. Saved results have explicit units and checks;
the complete application and arbitrary engineering models are not certified.
