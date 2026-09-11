# Single-solid compound capture regression

The [recipe and scientific scope](../../FREECAD_COMPOUND_CAPTURE.md) describe the
inputs and independent references. `record.zip` retains 115 files with an
internal SHA-256 manifest; `sha256.json` identifies the archive itself.

- `original-failure`: the original native profiling macro and report. A valid
  four-hole plate is a `Part.Compound`; its static capture fails on
  `CenterOfMass`. The process exits normally, but its report records failure.
- `installed`: development 0.1.0a3.dev7 installed in an isolated FreeCAD profile;
  16 native capture/rejection checks passed. Includes source provenance and ZIP.
- `workers`: the four earlier captures were reimported by the separate OCCT 8
  engine and all four motion calculations passed. Native analytic checks also
  passed. These captures were created directly from the edited modules before
  the extension version increment; their provenance records that fact.
- `static`: the installed native task passed on a rotated parametric cylinder
  inside a `Part::Compound`. Maximum analytic displacement error was
  `3.462137686674256e-11 m`; 106 GUI timer ticks occurred while polling the job.
  Save/reopen preserves the result; a `1e-10 mm` edit invalidates it. This run
  preceded the packaging-only version increment.
- `source`: the corrected capture modules, native recipes and packaging runner.

The SHA-256 hashes of both corrected capture modules agree across the installed,
worker and static records. Profile/cache directories and installed file copies
are excluded; the installation ZIP is retained. Absolute paths refer to the
recorded Linux experiment and are not portable installation locations.

The fixture's first worker attempt used a non-UUID identifier and was rejected;
it was corrected before producing the retained successful `workers` captures.
The original capture failure is distinct from that test-fixture error.

These records establish the stated analytic cases and refusal behaviour, not
arbitrary topology compatibility, a general FEM bound or full-kernel certification.
