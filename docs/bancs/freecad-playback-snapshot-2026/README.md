# Native Assembly playback snapshot record

See [behavior, reproduction and limits](../../FREECAD_PLAYBACK_SNAPSHOT.md).

- `before`: native run with the original dev8 host; 26 checks passed. Instrumented
  preview creation called the source snapshot twice on each of two calculations.
- `dev9`: final installed extension; 30 checks passed. Four preview recreations
  each use one validated snapshot, including recreation without a dirty flag.
  Source/module hashes and the final recipe hash are recorded. Native documents,
  captured calculations, logs and the installation ZIP are retained.
- `intermediate`: two earlier passing reports retain timing variability; their
  checks preceded the final version metadata and recipe-hash recording.
- `source`: final production module, native recipe and runner.

The archive excludes runtime profiles, caches and extracted installed copies.
`sha256.json` identifies the ZIP; its internal manifest identifies every payload
file. CRC and internal hashes were verified. Paths in reports identify this
local Linux experiment, not a portable installation. Timings are observations,
not an established general speedup; the robust result is the removed duplicate
snapshot with source validation and physical playback checks retained.
