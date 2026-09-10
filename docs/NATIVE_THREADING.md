# Native multithreading

**Source: kernel 0.20.0 / Studio 0.6.0a2.dev3.** Native parallel execution is an
architecture requirement. The implementation below is the first explicit
resource-management increment. The downloadable Studio 0.6.0a1 still contains
kernel 0.19.0; it does not contain these changes.

## Kernel execution

```python
from concurrent.futures import ThreadPoolExecutor
from vinkulum import ExecutionPool, Noyau

cpu = ExecutionPool(4)
models = [Noyau(executor=cpu), Noyau(executor=cpu)]
# Populate each independent model before submitting it.
# with ThreadPoolExecutor(max_workers=2) as callers:
#     results = list(callers.map(lambda m: m.simule(1.0, .001), models))
assert models[0].execution_pool_id == models[1].execution_pool_id == cpu.pool_id
assert cpu.threads == 4
```

A pool owns a fixed number of native Rayon workers. Models retain shared
ownership; deleting the Python pool variable does not invalidate them. The
same pool can service several calling threads without constructing another
worker set. Separate explicit pools have separate budgets: the caller must
account for their sum. This bounds Rayon compute workers, not every OS thread,
calling thread or third-party library in a process. The pool ID is process-local;
it is neither a global identifier nor a performance certificate.

Explicit pool execution currently covers `assemble`, `simule`,
`simule_multirythme`, `simule_em`, `_trace_em`, `modes_complexes` and `spectre`.
These numerical operations release the Python GIL. Assembly retains its complete
rollback on failure. Unrelated API operations, including `statique` and other
modal/matrix exports, have not all been migrated to this executor yet.
`Noyau()` without `executor` preserves the historical global Rayon path;
its `execution_pool_id` is zero. Concurrent mutation of one model is unsupported;
concurrent independent models are supported.

The element assembly already used Rayon; this version makes its executor
explicit and shareable. Indexed collection and ordered scatter retain the
existing summation order. The threshold of 48 bodies/elements and the sparse
parallel threshold of 2,000 unknowns remain in place. Small cases can be slower
with additional workers. A native thread pool does not make a sequential
algorithm parallel or override dependencies between time steps.

## Studio execution

Choose **Run → Native dynamics CPU threads…** for the next native calculation.
The choice is captured when Run is pressed; changing it does not change an
active or queued job. The default per-run allocation is the application budget.
Reduce it to admit multiple independent native calculations concurrently.

The application budget defaults to `max(1, os.process_cpu_count() - 1)` to leave
headroom for interaction. `VINKULUM_STUDIO_CPUS=4 vinkulum-studio` sets a smaller
explicit budget. CPU availability is an affinity-aware estimate, not a CPU
reservation: competing processes and container quotas can reduce real capacity.
Use the explicit setting for constrained environments. Each separate Studio
application has its own budget; there is no machine-wide broker.

Within an application, FIFO admission ensures that the sum of allocated native
threads does not exceed its budget. A queued job can be cancelled before it
starts. Failed startup, cancellation and completion release its allocation.
The lease remains held while the result is checked. Queue wait does not consume
a compute worker. No GUI callback waits for CPU admission.

Each worker constructs an `ExecutionPool` from its immutable request. Nested
BLAS/OpenMP pools are limited to one thread in this native-dynamics process.
The manifest records the requested threads, application budget, actual native
pool identity and worker CPU availability. The parent rejects a result whose
allocation differs from the captured request. This is a consistency check, not
protection against a malicious worker or a guarantee that every worker was busy.

Trajectory file reading, hash checks and scientific array validation now use a
Qt worker thread. Cancellation during validation preserves the previous result
and keeps the temporary files alive until the reader exits. Widgets and VTK
rendering remain on the GUI thread. Python-only validation still obeys CPython's
GIL; this change does not claim parallel Python execution. Closing a window
joins its bounded result reader before destroying its objects.

## Remaining work required by the architecture

| Path | Current status | Next required increment |
| --- | --- | --- |
| Native dynamics | Explicit Rayon pools, CPU admission, off-GUI result admission | Tune scheduling using measured workload sizes; extend coverage to remaining heavy kernel APIs |
| OCCT 8 / build123d | Isolated CAD process | Audit native operation-level parallelism and route the allocated CPU budget into supported OCCT operations |
| Gmsh | Explicit mesher thread count | Join the common admission budget and verify effective engine concurrency |
| CalculiX | Isolated worker with one-thread settings | Qualify threaded builds and sparse-solver behavior before increasing allocation |
| Pinocchio | Isolated native worker; single captured-state analysis | Use native batch/state parallelism where independent work exists; thread count alone cannot parallelize one recursive solve |
| GUI | Qt/VTK rendering on GUI thread; native-result validation moved out | Profile and migrate remaining expensive CAD/FEM/Pinocchio result admission |

Process isolation and native multithreading solve different requirements.
Until all engine adapters join admission, this budget covers **native dynamics**
only. It must not be described as a global Studio CPU limit.

## Verification and measurement

- Rust tests broadcast into distinct OS worker threads and verify that nested
  Rayon work executes within its pool. Shared and concurrent pools retain their
  identities and capacities.
- `python -m unittest vinkulum.test_execution` checks strict arguments, retained
  ownership, concurrent independent models, recovery after numerical failure,
  all three time integrators and equality of states/iteration counts/constraints
  across 1, 2 and 4 threads. Chains of 64 and 256 bodies cross the assembly and
  sparse parallel thresholds respectively. Equality is empirical for these
  cases, not a theorem for every model or machine.
- Studio tests cover FIFO bounds, cancellation before dispatch and from a
  synchronous busy signal, failed startup with queued recovery, real workers
  with captured allocations, and a worker-thread reader with a live GUI timer.
- [Measurements](bancs/native-threads-020/measurement.json) include all timed
  samples, CPU time, wall time, constraints, solver counters, source hashes and
  numerical equivalence checks. Reproduce with
  `python ci/measure_native_threads.py --output report.json`; use an otherwise
  idle machine and a fixed CPU affinity. No noisy timing ratio is a CI gate.

Underlying contracts: [Rayon ThreadPool](https://docs.rs/rayon/latest/rayon/struct.ThreadPool.html),
[Qt thread affinity](https://doc.qt.io/qt-6/threads-qobject.html).
