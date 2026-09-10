# Native CPU execution — qualification

Kernel 0.20.0 / Studio 0.6.0a2.dev3, Linux x86-64, CPython 3.14.7.

The [architecture contract](../../NATIVE_THREADING.md) defines exactly which paths use the new executor.

## Integration measurements

One native pool is alive at a time. The process is restricted to CPUs 0–7.
Pool construction and model creation are excluded; each case has an untimed warm-up
and five measured samples. All final states equal the one-thread reference.
These are observations on this host, with no claim about another machine.

| Bodies | 1 thread | 2 threads | 4 threads | 8 threads |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.359 ms (1.00×) | 0.393 ms (0.91×) | 0.309 ms (1.16×) | 0.343 ms (1.05×) |
| 64 | 19.430 ms (1.00×) | 18.856 ms (1.03×) | 18.190 ms (1.07×) | 18.060 ms (1.08×) |
| 256 | 385.009 ms (1.00×) | 389.851 ms (0.99×) | 380.105 ms (1.01×) | 383.994 ms (1.00×) |
| 900 | 172.436 ms (1.00×) | 167.929 ms (1.03×) | 151.433 ms (1.14×) | 147.510 ms (1.17×) |

These 20-step runs show modest gains for 900 bodies and essentially none for
256 bodies. Sub-millisecond measurements for one body are dominated by noise.
The non-monotone time across model sizes also precludes inferring an asymptotic
complexity from this table. This increment establishes bounded execution;
it does not demonstrate a universal acceleration or an advantage over other solvers.

An initial protocol retained all idle pools at once. It was discarded because
their worker spin periods could interfere between samples. The committed protocol
destroys each pool before the next allocation and records all timings, solver
counters, residuals and relevant source hashes in [measurement.json](measurement.json).

## Reproduction

```bash
taskset -c 0-7 python ci/measure_native_threads.py --output report.json
python -m unittest vinkulum.test_execution
PY=/path/to/installed-studio/bin/python bash ci/studio.sh
```

The kernel tests also exercise shared pools from multiple calling threads.
The desktop tests verify queue bounds, real worker provenance, failures,
cancellation and window shutdown. Timing ratios are not CI gates.

## Installed-artifact checks

- 90 Rust tests and 4 native Python execution tests passed.
- Studio: 157 tests, 147 passed and 10 optional native Pinocchio tests skipped.
  All 10 skipped cases passed in the separate 11-test Pinocchio suite.
- 80 Qt lifecycle tests passed across 10 repetitions in one process, including
  cancellation during validation and closing from the stage-change callback.
- All 56 Studio Python modules match the wheel and both installed environments.
  Wheel hashes and source/journal hashes are recorded in
  [qualification.json](qualification.json) and [module-identity.json](module-identity.json).

The first desktop run exposed an obsolete “Stopping…” message after synchronous
queued cancellation and a deferred command-capture assumption. Both were fixed.
A new close-during-validation test also exposed a native QObject-destruction crash;
receiver-bound Qt slots and application ownership of the scheduler replace the
capturing callback topology. The final lifecycle stress and complete suite pass.
These tests do not prove absence of every possible concurrency defect.

The standalone Linux 0.6.0a1 release is unchanged. This qualification covers the
new installed source wheels; no new standalone binary was built.
