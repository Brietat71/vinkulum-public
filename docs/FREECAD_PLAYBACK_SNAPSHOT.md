# Assembly playback snapshot reuse

Development extension **0.1.0a3.dev9** creates Assembly playback copies from the
same source snapshot used to validate their identity. Previously, creating an
Assembly preview performed two full snapshots on the GUI thread: one to compare
the source with the saved calculation, then another to obtain the shapes.

The returned shapes now remain local to that display invocation. Recreating a
preview always validates the source, including when the document observer has
not marked it dirty. No snapshot cache persists between frames and no geometry
identity tolerance changes. Existing source-change invalidation is retained.

## Native verification

Run from the repository root with the [qualified Linux host](FREECAD_HOST.md)
and [separate engine](FREECAD_ENGINE.md):

```sh
python3 ci/freecad_extension.py --recipe assembly-analysis \
  --freecad /path/to/FreeCAD-Vinkulum \
  --engine-python /path/to/engine/bin/python \
  --output /tmp/vinkulum-playback-check
```

The installed recipe runs a real double-pendulum calculation, saves and reopens
the native document and capture, then runs a second calculation. It records
snapshot counts and elapsed times while recreating each preview with and
without a dirty flag. Each recreation must take exactly one snapshot and
produce both moving solids. Displayed centres are compared with the calculated
poses; the grounded component and source geometry must remain unchanged.
Moving the Assembly or editing a connector must still invalidate playback.
The remaining checks exercise preview deletion and worker cancellation.

The [retained records](bancs/freecad-playback-snapshot-2026/README.md) include the
previous implementation's two-call observation and the corrected native run.
The qualification scope is this Assembly workflow, not arbitrary mechanisms.

## Performance limits

The deterministic improvement is removing one complete snapshot per preview
creation. Individual timings on the small example overlap across runs and are
not a statistically established speedup. CAD capture, result admission and
preview construction still execute on the GUI thread. This change does not
establish responsiveness for large CAD models or eliminate their potential
pauses. The native source document is not moved into an unqualified background
thread.
