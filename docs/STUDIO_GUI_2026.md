# Vinkulum Studio interface programme — September 2026

The objective is a major upgrade of the engineering workflow, from mechanism
authoring to inspection of traceable results. Qualification concerns real
application interactions and executable packages. Screenshots and visual style
alone do not establish parity with a mature CAD/CAE suite.

## Requirements and expected evidence

| ID | Required behaviour | Acceptance evidence |
|---|---|---|
| UX-01 | An uncluttered 3D area, rearrangeable/restorable panels, native menus and grouped commands | Model/result captures, 1280×800 and 1440×950 windows, layout restoration after closing |
| UX-02 | Searchable command palette, native shortcuts and keyboard access to essential functions | Real keyboard tests; unavailable commands cannot execute |
| UX-03 | Filterable Browser, consistent tree/scene/property selection and reversible visibility/isolation | Name/type search; hidden objects remain in the model and calculation |
| UX-04 | Structured Inspector, explicit units/frames, visible pending edits and validation without losing fields | Invalid input, selection change, cancellation and results arriving during editing |
| UX-05 | Precise 3D interaction: orientation, framing, projection, manipulation modes and navigation information | Qt/VTK tests, numerical and mouse manipulation, Linux and macOS rendering |
| UX-06 | Interactive curves, sample tables, playback controls and identified comparisons | Time synchronises scene/curve/table; each calculation's model and units remain visible |
| UX-07 | Clear document/calculation state, navigable diagnostics, accessible provenance and errors | Successful run, invalid model, failure, cancellation and preservation of previous results |
| UX-08 | Contrasting themes, display scaling, visible focus, accessible names and persistent layout | 100%/200% captures and keyboard workflows; accessibility limits documented |
| UX-09 | Responsive execution, enforced memory budgets and intact files/results | Existing and new tests, rendering and large-curve measurements |
| UX-10 | Current documentation and a genuinely executed package; local Linux builds first | Extracted/tested Linux archive, bundle report and fingerprint; Apple Silicon DMG qualified separately |

GUI-01–GUI-08 and UI-01–UI-06 of the
[engineering specification v1.1](../outputs/Cahier_des_charges_suite_ingenierie_Vinkulum.md)
also apply. An interface redesign does not change the kernel's scientific
guarantees. Passing model-input checks does not certify a trajectory. Future
CAD, collaboration and multiphysics functions need implementations and evidence
beyond a change of presentation.

## Current status

Work remains in progress. Studio 0.2.0 is the historical baseline before this
redesign, not the current desktop release. The
[public Linux 0.5.0 preview](https://github.com/Brietat71/vinkulum-public/releases/tag/studio-v0.5.0-linux)
includes CAD, native dynamics, CalculiX studies and checked result reopening.
Its extracted application was tested through both normal startup and actual
result-opening interaction.

The source subsequently adds the [Pinocchio operator workspace](PINOCCHIO_OPERATORS.md)
and [parametric CAD feature editing](STUDIO_CAD_HISTORY.md). Their documentation
distinguishes source functionality from features in published binaries. Every
requirement above must be checked against its actual evidence before declaring
the overall upgrade complete.

## Design direction and primary references

The requested direction goes beyond recolouring forms: the scene should dominate,
values should remain readable, commands should use consistent icons, and density
should suit the active workspace. Real captures help identify clipping, alignment
errors and weak visual hierarchy.

References reviewed on 9 September 2026:

- [FreeCAD 1.1 release](https://freecad.github.io/Website/download/releases/1-1/): precise transforms, navigation, search, themes and feedback inform UX-03–UX-05. These directions are not evidence that Studio implements them completely.
- [ParaView 6 customisation](https://docs.paraview.org/en/v6.0.0/ReferenceManual/customizingParaView.html): property search, basic/advanced separation, persistent preferences and restoration inform UX-01/UX-04.
- [Qt accessibility](https://doc.qt.io/qt-6/accessible.html): keyboard navigation, size adaptation, contrast and component semantics inform UX-02/UX-08.

The requested products guide specific choices. They are not integrated libraries
or a promise to reproduce every capability of those products.

| Primary reference | Principle used for Studio | Application in the redesign |
|---|---|---|
| [Plasticity interface](https://doc.plasticity.xyz/plasticity-essentials/plasticity-interface/user-interface-overview), [command palette](https://doc.plasticity.xyz/plasticity-essentials/plasticity-interface/command-palette) | A dominant scene, compact tools and direct command access | Graphite interface, vector icons and searchable palette |
| [Blender tools/workspaces](https://docs.blender.org/manual/en/4.5/interface/tool_system.html) | An explicit active tool and work context | Model/Simulate/Inspect workspaces; select/move/rotate modes |
| [Shapr3D adaptive interface](https://support.shapr3d.com/hc/en-us/articles/7873882619548-Adaptive-user-interface) | Actions relevant to the selection | Contextual tools and explicitly disabled unavailable commands |
| [Abaqus/CAE management and visualisation](https://www.3ds.com/fileadmin/Products/Simulia/PDF/datasheets/Abaqus_CAE_Datasheet.pdf) | Separate authoring, execution and examination | Read-only captures, run history and provenance |
| [Onshape tool search](https://cad.onshape.com/help/Content/Home/search_tools.htm) | Discover a command without knowing its location | Word search, visible shortcuts and keyboard activation |
| [NX command access](https://blogs.sw.siemens.com/designcenter/designcenter-x-nx-tips-and-tricks-copilot/) | Find relevant operations in a rich application | Shared command registry and contextual access; no simulated AI assistant |
| [STAR-CCM+ 2606](https://blogs.sw.siemens.com/simcenter/simcenter-star-ccm-2606-released/) | Inspect differences between simulations | Overlaid series, object identity and differences in settings |
| [SolidWorks shortcuts and contextual menus](https://blogs.solidworks.com/products/solidworks/useful-keyboard-shortcuts-workflow-customizations-solidworks/) | Reduce the distance to common tools | S palette, native menus and tools near the scene |
| [3DEXPERIENCE tree, scene and actions](https://3dswym.3dexperience.3ds.com/wiki/solidworks-news-info/getting-started-with-3dexperience-simulation-solidpractices_rFZtKhrBSdO2cIlYITOksg) | Connect selection with action context | Stable identities shared by Browser, scene and Inspector |
| [Autodesk Fusion interface](https://help.autodesk.com/view/fusion360/ENU/?contextId=LP-STEPS-P13N-SNP-GS-OTH-CRD-1) | Direct workspace and spatial navigation access | Persistent workspaces, interactive orientation and isolation |
| [Rhino Gumball](https://www.rhino3d.com/en/docs/guides/user-guide/gumball-basics/), [layouts](https://www.rhino3d.com/features/user-interface/window-layouts/) | Combine direct manipulation and numerical precision | Manipulator modes, X/Y/Z components and restorable panels |
| [Creo command search](https://support.ptc.com/help/creo/creo_pma/r12/usascii/fundamentals/fundamentals/to_search_a_command.html) | Locate commands by name and help | Accent-insensitive search through one command registry |

Established principles are not presented as inventions of 2026. These references
document interface patterns, not a comprehensive local evaluation of commercial
products or an assertion of functional parity.

## Implemented interaction baseline

- Reorganised Qt/VTK editor, native menus, workspaces and panels.
- Global Ctrl/Cmd+K and contextual S palettes; Browser filtering.
- Graphite/light themes, original vector icons and visible focus.
- Separate vector components. Compact numeric presentation retains each original
  binary value until that component is edited.
- Interactive orientation, selection/move/rotate, visibility and isolation
  without changing document or calculation inputs.
- Session history of up to eight calculations within a 128 MiB array budget.
  Older calculations are evicted at the limit; exports remain explicit.
- Curves with sample selection, zoom, pan and adjustable playback; lazy sample
  tables avoid copying every cell.
- Comparisons use UUIDs and physical channels on each run's native time grid.
  Renamed objects retain identity; different identities are not matched by array
  position. No resampling or interpolated numerical difference is implicit.
- Save/discard handling, document dirty state, retained invalid fields and
  protection of pending edits when a calculation finishes.

## Historical qualification records

The first Linux redesign baseline passed 38 tests: 29 previous cases and nine
new cases. GUI tests also watch deferred Qt exceptions because an **OK** unittest
summary alone does not detect every callback failure. A VTK-wrapper cleanup
defect discovered in the logs was corrected.

The [200% Linux report](bancs/studio-gui-2026/linux-hidpi.json) includes
[modelling](bancs/studio-gui-2026/modeling-dark.png),
[comparison](bancs/studio-gui-2026/comparison-dark.png) and a
[small light-theme window](bancs/studio-gui-2026/comparison-light-1280.png).
The large capture is 2880×1900 pixels for a logical 1440×950 window. With three
bodies and two curves, software rendering measured 36.2 ms median and 56.1 ms
at the 95th percentile in that environment. This is not a GPU measurement on
the user's Mac or a frame-rate guarantee. The recipe also checks that playback
controls remain present in a logical 1280×800 window.

Subsequent source refinements added a faded, hideable metric XY grid, less
saturated lighting, a compact orientation widget, and revolute symbols at both
actual attachment positions. Construction lines appear for the selected joint.
Perspective/orthographic projection is selected through Views. These changes
do not alter mechanical volumes or kernel inputs.

Those 38 tests, four 200% captures and two native calculations passed again.
An attempted FXAA setting produced a black main scene while leaving the
orientation widget visible; it was removed. The bundle check now samples the
reference scene image to catch that failure. It does not measure overall visual
quality. The historical Linux 0.3.0 archive preceded those final refinements;
presentation edits subsequently ran from source before milestone packaging.

## Precision and density pass — 10 September 2026

Studio 0.4.0 reduced command, Browser and property spacing. CAD commands acquired
labels, and menu indicators retained their own space. Inspector fields choose
a readable preview for the available width while storing the exact value
separately; focus and tooltips expose that value. Renaming a body does not rebuild
its rotation or round its position and mass.

The [CAD/interface qualification](bancs/studio-cad-040/README.md) passed 48 tests,
including focus, resizing and small components in scientific notation. Captures
at 1280×844 and high DPI retained complete exponents. These checks are not a
usability study or proof of parity with mature engineering applications.

Studio 0.4.1 translated the interface, and 0.4.2 qualified normal startup through
the real application entry point. Studio 0.5.0 added the CalculiX workspace,
checked saved-result reopening and framing of visible bodies/attachments.
The [Pinocchio installed-package record](bancs/studio-pinocchio-060/README.md)
covers its separate worker, captured-state editor and operator inspection.
The [CAD feature guide](STUDIO_CAD_HISTORY.md) explains subsequent parametric
previews, source-file compatibility and remaining topological-reference work.

## Still required before closing the overall objective

Final contextual-tool ergonomics, broader accessibility and display-scale
coverage, measured performance on representative models, robust layout
restoration across multiple monitors and a current Apple Silicon package remain
open. Each new source milestone needs its own extracted Linux application checks
before it becomes a downloadable standalone release.

Constrained sketches, persistent face/edge references, collaboration and broader
multiphysics workflows also remain substantive implementation work. These
interaction records do not declare the full engineering-suite objective achieved.
