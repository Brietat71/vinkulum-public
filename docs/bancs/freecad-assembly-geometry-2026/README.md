# Native Assembly geometry admission record

See the [recipe, analytic budgets and limits](../../FREECAD_ASSEMBLY_GEOMETRY.md).
`record.zip` contains 74 files plus an internal SHA-256 manifest; `sha256.json`
identifies the complete archive. CRC and all internal hashes were checked.

- `original-failure`: native FreeCAD observation that extra faces, edges and
  vertices were accepted and discarded by the original extraction function.
- `native-geometry`: installed development 0.1.0a3.dev8, 17 passing checks for
  analytic world properties, native links, source preservation, explicit
  rejection, repeated snapshots and the full Assembly capture boundary.
  The restored capture also passed actual OCCT 8 import and motion calculation
  in the separate engine; its command, outputs and log are retained.
- `native-motion`: the existing 26-check native Assembly analysis recipe passed
  with the same corrected extraction module. This run preceded the dev8
  packaging version increment; its installation manifest records that version.
- `recompute-diagnostic`: the initially overstrong expectation of an identical
  whole-Assembly snapshot after recomputation failed. The only observed snapshot
  difference was Lower's geometry hash; its BREP location changed from zero to
  about -2.19047482849133e-29 mm. The report remains marked failed. The final
  recipe requires exact repeated snapshots without edits, restored ground
  geometry and admission of a new capture; it does not assert bit-identical
  native solver recomputation or weaken production identity checks.
- `source`: the corrected module and reproducible native recipe/runner.

The corrected extraction module's hash matches both successful installed
reports. Cache/profile/installed copies are excluded; installation ZIPs and
source manifests remain. Absolute paths identify this Linux experiment, not a
portable installation. These records do not certify general topology handling,
nested mechanisms, arbitrary constraints or full-kernel correctness.
