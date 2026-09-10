# Studio dev6: measured CAD process reuse

Qualification on 2026-09-10, Linux x86-64/X11, Python 3.14.7, an installed
Studio **0.6.0a2.dev6** wheel based on main `07cf6f7` (including dev5 background
admission). Vinkulum 0.20.0, PySide6 6.11.2, VTK 9.7.0,
`cadquery-ocp-novtk` 8.0.1.0.0 and adapted build123d 0.11.1+vinkulum.occt8.
The [implementation contract](../../CAD_PROCESS_REUSE.md) describes ownership,
CPU admission, transaction identity and failure handling.

## Actual desktop latency

The recipe opens the real editor and CAD dialog, clicks Create/Modify, waits for
the accepted document and explicitly renders its viewport. The timer includes
request capture, CPU admission, process startup when applicable, geometry,
response admission on the background thread, document commit and rendering.
It ends before the subsequent undo/redo checks and screenshot capture.

Each row below is the median of seven operations. Every operation starts from
the same plate and bore geometry. Box: 100 × 60 × 20 mm; cut: central through bore
of radius 10 mm; fillet: 1 mm on the plate. Density is 7800 kg/m³. The worker
receives two native threads; the GUI admission budget is two, BLAS is restricted
to one. No other qualification job ran during these measurements.

| DPI | Operation | Fresh process, ms | Reused process, ms | Fresh / reused |
| --- | --- | ---: | ---: | ---: |
| 100% | Box | 2638.36 | 307.59 | 8.58 |
| 100% | Through bore | 2588.90 | 289.30 | 8.95 |
| 100% | Fillet | 2694.84 | 369.79 | 7.29 |
| 200% | Box | 2969.51 | 589.99 | 5.03 |
| 200% | Through bore | 2812.82 | 578.32 | 4.86 |
| 200% | Fillet | 2922.41 | 612.61 | 4.77 |

The cached path's first box, excluded from warm medians, took **2111.78 ms** at
100% and **2266.08 ms** at 200%. These are two observations, not a startup-time
distribution. The one-shot control uses the same interpreter and geometry code,
selected through the controller's explicit executable override. Measurement
order was fresh 100%, cached 100%, cached 200%, fresh 200%; OS file caches were
not flushed. These are workload measurements, not portable performance limits.

Across the 86 captured bodies, the verifier checks independent analytic box and
through-bore masses and inertias, plus cross-mode/cross-DPI volume, position,
dimensions, full inertia and tessellation equivalence. The fillet comparison is
an equivalence check against the fresh result, not an independent analytic fillet
oracle. Each operation creates one undo transaction, and undo/redo restores the
expected bodies. Request and body files retain their SHA-256 fingerprints.

## Memory and lifecycle

| Mode | GUI initial / final RSS, MiB | Observed worker peak RSS, MiB | Retained idle worker RSS, MiB |
| --- | ---: | ---: | ---: |
| Fresh 100% | 393.04 / 419.78 | 450.71 | none |
| Cached 100% | 393.11 / 413.02 | 454.59 | 454.59 |
| Cached 200% | 566.47 / 603.13 | 454.57 | 451.96 |
| Fresh 200% | 566.66 / 609.22 | 450.80 | none |

RSS samples run every 20 ms while the CAD dialog is active, plus operation
boundaries. They can miss short peaks. Final GUI RSS is recorded after undo/redo
and screenshot capture; it can exceed the in-operation samples because those
steps allocate additional display buffers. RSS is not unique physical memory
and cannot simply be summed into a precise system memory cost.

Each cached run used **one PID for all 22 operations**, with distinct transaction
IDs and no CPU reservation between operations. After the last checks/capture,
the 15-second idle timer had 14.028 s and 12.862 s remaining. Actual expiry and
process reaping were observed after 14.092 s and 12.948 s, respectively, with zero
allocated CAD CPUs. GUI RSS stayed at its final value after worker expiry.
The recipe found **zero retained closed CAD dialogs** after each apply. An earlier
pilot exposed one additional hidden dialog per operation; the editor now copies
the admitted result before scheduling dialog destruction. These short sessions
do not establish a long-term absence of memory growth.

## Reproduction and evidence

Run from a checkout with the installed dev6 CAD environment and Linux Xvfb:

```sh
env QT_QPA_PLATFORM=xcb QT_SCALE_FACTOR=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  xvfb-run -a -s '-screen 0 1600x1100x24' "$PY" \
  ci/studio_cad_reuse_recipe.py /tmp/cad-fresh-100 --mode fresh --repeats 7
```

Use `--mode cached` for the reusable worker. At scale 2, use a 3200×2200 screen
and a separate output directory. The recipe fixes the CPU budget and the default
idle interval, captures all installed Studio module hashes and records its own
hash plus the independent-reference helper's hash. All 62 installed Studio
modules matched the implementation sources used for these measurements.

[captured-gui.zip](captured-gui.zip) contains the four complete runs: requests,
body payloads, per-operation timings, memory samples, saved final projects and
100%/200% screenshots. Its SHA-256 is
`dd711edd13477b3b2b6d4fb7cd80fc7a813eacdd819cb043091635a89d0877b2`.
[qualification.json](qualification.json) exposes the run summaries without
extracting the archive. To recheck the captured geometry without OCCT or Qt:

```sh
unzip docs/bancs/cad-reuse-dev6/captured-gui.zip -d /tmp/cad-reuse-captures
"$PY" ci/verify_cad_reuse.py /tmp/cad-reuse-captures
```

The archive regression also changes a captured box mass and checks rejection.
The full Studio suite passed **202 tests with 17 skips** before the final partial
protocol-response guard; the final **14-test lifecycle suite passed** and covers that guard as
well as reuse, concurrency, cache invalidation, idle expiry, cancellation,
timeout, crashes, native operation errors and background-validator lifetime.
The skips are 16 optional numerical Pinocchio tests in the GUI interpreter and
the separately invoked keyboard recipe; external Pinocchio worker scenarios ran.
Native Pinocchio operators are unchanged by this contribution.

`ci/local.sh --bancs` passed, including **46/46 mechanical** and **9/9 contact**
cases; three optional Exudyn comparisons were skipped. This run preceded the dev5
integration; its native/Python kernel, proof and CI inputs are unchanged by that
integration or by this Studio-only feature. See [ci-local.log.gz](ci-local.log.gz)
and [studio-tests.log.gz](studio-tests.log.gz).

The keyboard recipe exercises creation, feature editing, undo/redo and saving at
100%/200% with 420/280-pixel inspectors. Initial runs on the preceding integration
were timing-sensitive around modal activation; a baseline dev4 run passed. The
recipe now waits for native window activation before sending keys, while retaining
the original field-specific focus assertions and without moving focus itself.
The initial failure and baseline logs are retained alongside the final recipe log.
All four combinations passed in the final run. See [keyboard.log.gz](keyboard.log.gz),
[lifecycle-tests.log.gz](lifecycle-tests.log.gz), and
[archive-test.log.gz](archive-test.log.gz).

Qualification covers an installed Linux source wheel. The package's
`--cad-worker --service` dispatch is checked separately without GUI startup.
[keyboard-and-cli.zip](keyboard-and-cli.zip) retains the four keyboard projects
and captures, dependency versions, and the two-request package-dispatch probe.
Frozen Windows/macOS bundles, long edit sessions and large STEP models remain
separate qualification work. The keyboard recipe covers small-screen access;
the performance captures are not a complete visual audit of all DPI layouts.
