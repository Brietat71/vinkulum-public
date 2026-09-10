# FreeCAD engine installation — fresh Linux environment

The [installation command](../../../ci/freecad_engine.py) was exercised in a new
directory on Linux x86-64 / Ubuntu 24.04. It created its own Python environment,
installed the repository kernel and adapted CAD dependencies, checked the
dependency graph and ran the two actual external-worker regressions.

[Retained records](record.zip) · [SHA-256](SHA256SUMS) ·
[Installation guide](../../FREECAD_ENGINE.md)

The archive contains the executed installer, its source hash, installation steps,
the complete installed package-version inventory, worker-test output and the
actual FreeCAD qualification report using this newly installed engine. The
installer was run before committing this increment: `source_dirty: true` is
intentional. Its recorded base is `818bbaa`, and the executed installer hash is
`bffb09280f1462a3a0485f491940b9cc06c302d56b9ee7be75fadb234926d4b3`.

The resolved environment uses Python 3.14.7, Vinkulum 0.20.0, Studio backend
0.6.0a2.dev6, OCP binding 8.0.1.0, build123d 0.11.1+vinkulum.occt8 and
ocpsvg 0.6.0+vinkulum.occt8. The dependency check succeeds. Both worker tests
pass: a real native pendulum calculation agrees with the independent reference,
and a doubled captured mass is rejected before simulation. The recorded tests
complete in approximately five seconds; this is not an installation benchmark.

The public **FreeCAD extension 0.1.0a2** was then installed in a fresh official
FreeCAD 1.1.3 process, configured with this environment's `venv/bin/python`.
All 18 document, mechanics and lifecycle checks pass. This tests the newly
created installation rather than reusing the earlier manually prepared engine.
The GUI host retains OCCT 7.8.1; OCCT 8 stays in the separate engine process.

The installer refusal checks launch its real CLI against an existing directory
containing an environment marker and against a new destination without `uv` on
`PATH`. The first preserves file contents and timestamps; the second creates no
destination. These two checks are part of `ci/local.sh`. They do not download
packages or rebuild the kernel on each CI invocation.

This is a source-installation recipe, not a bundled binary or a globally locked
dependency graph. It requires existing build tools and network access. Primary
CAD versions and source hashes are pinned, while the resolved transitive set is
recorded. The current backend distribution still installs Qt/VTK dependencies.
Other operating systems, Linux distributions and CPU architectures are not
qualified by this record. The scientific scope remains the existing pendulum
reference and FreeCAD extension cases; no general trajectory guarantee is added.
