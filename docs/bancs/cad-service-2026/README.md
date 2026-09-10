# CAD service experiment — 10 September 2026

The retained run uses the installed Studio **0.6.0a2.dev4** wheel built from
main `efc76921d1d6b618b754bfeb581d282ae2b52208`, with kernel 0.20.0,
Python 3.14.7, OCCT 8.0.1.0.0, build123d 0.11.1+vinkulum.occt8 and NumPy 2.5.3.
All 58 installed Studio Python modules were compared byte-for-byte with that
source tree. The benchmark ran outside the checkout without `PYTHONPATH` on
Linux x86-64, after the local CI process had finished.

The [method and integration contract](../../CAD_SERVICE_EXPERIMENT.md) describe
the fixtures, SI reference equations, tolerances and limits. There are seven
paired repetitions per case, with alternating fresh/resident order, two native
OCCT threads and single-thread BLAS. The same production worker entry point is
called in both modes. This experiment does not change Studio's worker lifecycle.

| Case | Fresh median [ms] | Fresh min–max [ms] | Resident median [ms] | Resident min–max [ms] |
|---|---:|---:|---:|---:|
| Box | 2290.55 | 2274.25–2311.84 | 7.79 | 6.29–8.79 |
| Through-bore cut | 2309.02 | 2279.75–2320.85 | 23.09 | 20.79–31.79 |
| All-edge fillet | 2388.14 | 2365.82–2401.92 | 114.21 | 104.14–138.66 |

Resident rows exclude initial startup and final shutdown. The first resident
box, including startup and first imports, takes **1799.63 ms**; final graceful
shutdown takes **501.90 ms**. Fresh rows include termination and parent result
reading/validation. Fresh imports have medians near 1.69 s; observed exit waits
have medians near 0.49 s. Timing is observed through Linux process descriptors,
without timed-wait polling delays. GUI dispatch, rendering and undo/apply are
outside this experiment, so the table is not a whole-application speedup claim.

The resident process reaches **456.21 MiB** RSS and retains its native libraries
between requests. Its RSS grows from about 455.68 to 456.21 MiB across the measured
sequence. Fresh-worker observed peaks range up to 443.73 MiB for the box,
449.55 MiB for the cut and 452.35 MiB for the fillet, before termination releases
the process. Counters are sampled after operations; the short sequence cannot
rule out a longer-term leak. The parent reaches about 49.5 MiB before offline
archive verification. Raw counters, phase times, byte counts and load averages
remain in the record.

All 45 retained benchmark responses pass structural checks, captured-request
identity, native-allocation checks and numerical comparisons. Box and cut mass
and inertia agree with independent analytic references. Fillet equivalence is
checked between fresh and resident execution, including display tessellation;
no independent analytic fillet reference is claimed.

The actual invalid 1,000 mm fillet is rejected and the same process subsequently
produces the correct box. Injected cancellation, a 100 ms deadline and SIGABRT
preserve the previous accepted result; each replacement process produces the
reference box. Observed settle times are 41.77, 136.32 and 141.91 ms respectively
(the deadline case includes its 100 ms wait). Replacement startup plus its first
box takes about 1.79–1.80 s. These are tests of the experiment's supervisor;
the injected cases are not identified upstream OCCT faults or GUI tests.

## Retained files

- `qualification.json`: versions, complete Studio source fingerprints, harness
  fingerprint, request/response file hashes, all measurements and outcomes.
- `captured-results.zip`: original inputs, all worker responses and logs, plus
  the same qualification record. Approximately 1.39 MB compressed.
- `qualification.log.gz`: completed live experiment and its offline recheck.
- `ci-local.log.gz`: completed `ci/local.sh`, including native tests, Lean proofs,
  reference archives and documentation checks. Three optional Exudyn checks are
  skipped. The fast-forward over #20 changes no native, proof or CAD-worker
  files from the tree checked by that process.
- `studio-tests.log.gz`: 184 Studio tests, OK with 17 skips (16 numerical
  Pinocchio checks belonging to the separate engine environment and one opt-in
  keyboard recipe). The new archive regression passes. The external worker is
  installed for the GUI/controller checks; numerical Pinocchio tests and the
  optional keyboard recipe were not rerun for this measurement-only change.

Reproduce the live run with `ci/qualify_cad_service.py` as described in the method.
To verify this retained archive without executing CAD:

```sh
python -m zipfile -e docs/bancs/cad-service-2026/captured-results.zip /tmp/cad-service-archive
python ci/qualify_cad_service.py /tmp/cad-service-archive --verify
```

The Studio regression suite also opens this archive in an engine-free subprocess
and checks that changing a retained mass invalidates the record. The offline
check verifies file and numerical consistency; it does not reproduce timings,
rerun failure injection or authenticate the original measurements.
