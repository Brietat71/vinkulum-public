# Build the patched FreeCAD host on Linux

This development recipe completed a build from a fresh upstream clone and a
full resume with the final installer. The resulting launcher passed 99 checks
across five native FreeCAD scenarios. The [retained build and qualification
record](bancs/freecad-host-build-2026/README.md) states the tested scope.

The static workflow uses two FreeCAD 1.1.3 persistence fixes: BREP stream
formatting and exact placement round trips. See the
[geometry identity contract](FREECAD_GEOMETRY_IDENTITY.md) for their scope.
The engine remains separate and uses OCCT 8; this host uses its locked OCCT 7.8.1
SDK. This command builds FreeCAD; it does not install the Vinkulum extension,
the engine, Gmsh or CalculiX.

## Build

Requirements: Linux x86-64, Python 3, Git, Pixi **0.80.0**, network access and
enough disk/memory for the complete FreeCAD SDK and source build. Use a new
user-owned directory whose absolute path will remain stable. From this checkout:

```sh
python3 ci/freecad_host.py /path/to/new/freecad-host --jobs 4
```

If Pixi is not on PATH, pass `--pixi /absolute/path/to/pixi`.
The script verifies the FreeCAD commit and SDK lock hash before applying the
reviewed patches. Patch line endings are preserved in Git because the upstream
checkout uses CRLF. It installs the locked SDK in the new directory, configures
CMake and builds with the requested concurrency (1–32 jobs; linking limited to
one job). It performs no system installation.

A successful build writes a `FreeCAD-Vinkulum` launcher. Use that launcher rather
than the raw executable so the locked native environment is selected. The
installation uses absolute paths and is not a relocatable binary bundle.

## Preparation and controlled resume

To fetch and patch sources without installing the SDK or starting compilation:

```sh
python3 ci/freecad_host.py /path/to/new/freecad-host --prepare-only
```

After preparation has succeeded, continue in the same directory:

```sh
python3 ci/freecad_host.py /path/to/new/freecad-host --resume --jobs 4
```

An existing destination is refused without `--resume`. Resume checks the pinned
commit, lock, patch identities and tracked source diff against the retained
record. Changed sources are refused. Failures retain files and logs for diagnosis;
a failure before successful preparation may require choosing a new directory.
Do not edit the source tree or patch files while building.

An exclusive lock prevents simultaneous installer invocations in one destination.
After an installer crash, resume also checks for surviving commands in its
recorded Linux process session, including the separate groups created by Ninja.
Wait for those commands to finish before resuming. An interrupted launch without
a recorded session cannot be safely resolved automatically and requires a new
destination. Do not remove the lock file to force a second invocation. Commands
that deliberately detach into another session are outside this guard's scope.

`build.json` records invocation hashes, commands, status and executable hash;
`build.log` retains tool output. A `built` status means compilation succeeded,
not that mechanics or GUI behaviour have been certified. Run the native
[qualification recipes](FREECAD_GEOMETRY_IDENTITY.md) with the resulting launcher
and the separate qualified engine before relying on the installation.

## Qualify the generated launcher

After the build finishes, run these checks from the Vinkulum checkout. Replace
the paths with the new host, a separately installed engine (adapter
`0.6.1.dev2` for this candidate), the OCCT 8 Gmsh executable and CalculiX.
The output directories must not already exist. Linux qualification uses
`xvfb-run`; no interactive desktop session is required.

```sh
qualify_host() {
    python3 ci/freecad_static.py \
        --freecad /path/to/new/freecad-host/FreeCAD-Vinkulum \
        --engine-python /path/to/engine/venv/bin/python \
        --gmsh /path/to/gmsh --ccx /path/to/ccx "$@"
}
qualify_host --recipe static-pose --output /path/to/checks/pose
qualify_host --recipe static-identity --output /path/to/checks/identity
qualify_host --recipe static-curved-result --rotated --output /path/to/checks/curved
qualify_host --recipe static-reopen \
    --result-file /path/to/checks/curved/result.FCStd \
    --output /path/to/checks/reopen
qualify_host --recipe static-task --output /path/to/checks/task
```

Run each command only after the preceding command succeeds. These checks cover
placement round trips, the declared geometry identity contract, an actual
rotated-cylinder solve, reopening its saved result in a separate process and
the static task's interaction/process lifecycle. Keep each directory's
`provenance.json`, `console.log` and `static-check.json` with the host's
`build.json` and `build.log`. Passing these bounded recipes does not establish
general solver accuracy or qualify every FreeCAD workbench.
