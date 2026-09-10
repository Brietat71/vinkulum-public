# FreeCAD extension 0.1.0a1 — Linux qualification

The [extension](../../../apps/freecad/README.md) was installed from its ZIP into
an isolated FreeCAD user `Mod` directory and exercised in the real FreeCAD GUI.
It adds a global menu and a native task panel; selecting PartDesign keeps the
Vinkulum menu available. No custom Studio GUI is opened.

![Actual installed extension in FreeCAD](freecad-extension.png)

## Runtime and identity

The host is the previously checked official FreeCAD 1.1.3 x86-64 AppImage,
Python 3.11.14, Qt 6 and OCCT 7.8.1 on Linux/X11 under Xvfb. The separate
engine uses Python 3.14, Vinkulum 0.20, Studio's backend modules, OCCT 8.0.1
and adapted build123d. The earlier [runtime locator and checksum](../freecad-bridge-2026/README.md)
identify the unchanged FreeCAD runtime.

[precommit-extension-check.json](precommit-extension-check.json) records the
installed test before the source commit. Its `source_dirty: true` is intentional;
the packaged application files are individually fingerprinted. The final public
extension is built with `--require-clean` and qualified again after extraction;
its release report and ZIP manifest identify that clean commit. The screenshot
is an actual FreeCAD window, not a mockup or a GPU performance benchmark.

## Checks and scientific scope

The recipe checks installed-file hashes, menu availability after native
workbench activation and an actual pendulum calculation started through the
panel's Run button. It checks the 3.744 kg captured mass, explicit thread
allocation and the displayed copy's actual centre against a native pose.
The original source's shape and placement remain unchanged by playback.

A captured result reopens in a fresh task panel with an unavailable engine path.
Saving during playback removes the temporary motion object before native
`.FCStd` serialization and restores the source visibility. Closing and reopening
that actual saved design preserves admission of its captured calculation.
Changing the source placement refuses old-motion playback.

A second calculation applies a 90° rotation around world Z and a translation
of (250, -100, 500) mm to the source and pivot. The hinge axis becomes (-1, 0, 0).
The run uses one engine thread instead of two. All native positions are compared
with `Q p + t`, and all body rotations with `Q R Qᵀ`. The observed maximum
component differences are `1.6461e-13 m` and `3.0182e-13`, respectively, against
budgets of `1e-8 m` and `1e-8`. This checks rigid-frame transport for the recorded
case, not a general trajectory error bound. The original
[independent finite-section pendulum reference](../freecad-bridge-2026/README.md)
and two external-worker regression tests remain available separately.

The lifecycle checks deliberately fail process startup, attempt a second task
while one owns admission, cancel a real running process behind a filesystem
barrier, and close the source document during a job. Admission is released after
the worker settles. Cancellation and failed startup retain the previous captured
result. A closed panel cannot receive a late result.

This first host covers one top-level rigid solid, one explicit revolute joint,
rest initial velocities and gravity -Z. It does not qualify general assemblies,
FEM controls, other operating systems, arbitrary long trajectories or all desktop
lifecycle combinations. STEP capture still uses FreeCAD's GUI thread; large-solid
capture latency is not established by this small example. Geometry identity is
conservative: an edit/recompute can require a new capture even if the part looks
the same. No solver tolerance was relaxed for this interface.

## Reproduce the installed-archive check

With the qualified FreeCAD runtime, an installed Vinkulum CAD environment and
`xvfb-run`/`xauth` when no graphical session is available:

```sh
python3 ci/freecad_extension.py \
  --freecad /absolute/path/to/FreeCAD/AppRun \
  --engine-python /absolute/path/to/vinkulum-environment/bin/python \
  --archive /absolute/path/to/Vinkulum-FreeCAD-0.1.0a1.zip \
  --output /tmp/new-freecad-extension-check
```

The output directory must be new. The recipe creates isolated user configuration,
data and cache directories, installs the archive, supplies a cancellable barrier
process and runs [qualify_extension.FCMacro](../../../apps/freecad/qualify_extension.FCMacro).
It leaves the exact package, report, console log, screenshot, saved `.FCStd`
model and native calculation captures for inspection. Timeout handling terminates
only the process group started by this qualification.

To produce the release ZIP from committed source:

```sh
python3 apps/freecad/package.py /tmp/Vinkulum-FreeCAD-0.1.0a1.zip --require-clean
```

The extension ZIP includes the original example, installation instructions,
licences and packaged-file hashes. FreeCAD and the separate Vinkulum Python
runtime are not bundled. The native `.FCStd` design and captured calculation
folders remain separate, independently inspectable files.
