# Keyboard CAD qualification — 10 September 2026

Recipe and scope: [Studio keyboard workflow](../../STUDIO_KEYBOARD.md).
Source baseline: public `bd0c4a0`, with the accompanying inspector-focus fix.

| Desktop scaling | Inspector widths | Outcome | Evidence |
|---|---|---|---|
| 100% (DPR 1) | 420 / 280 logical px | Both pass | [log](100/recipe.log), [normal](100/keyboard-420.png), [narrow](100/keyboard-280.png) |
| 200% (DPR 2) | 420 / 280 logical px | Both pass | [log](200/recipe.log), [normal](200/keyboard-420.png), [narrow](200/keyboard-280.png) |

[Environment versions](versions.txt.gz). The [ordinary Studio suite](studio-suite.log.gz)
ran 137 tests successfully, with 14 skips: 13 require separately configured external
engines; one is this opt-in Openbox recipe, executed separately above. These are
actual test outcomes for this snapshot, not a claim that skipped paths were tested.

The captures use the native X11 window (including VTK). Position X is focused,
showing its editing outline and horizontally scrollable full-precision value.

The [required local CI](ci-local.log.gz) also completed successfully, including Rust,
kernel verification, archive checks and Lean proofs (`CI locale OK`).
