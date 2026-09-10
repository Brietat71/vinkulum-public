# CAD process reuse: measurements and integration contract

Small CAD operations can spend more time importing the geometry stack than
building the requested solid. The [original CAD recipe](bancs/studio-cad-040/README.md)
identified this possibility. The [retained experiment](bancs/cad-service-2026/README.md)
compares fresh processes with sequential reuse of one process, invoking the
same production `cad_worker.main` entry point for each request.

This is the measurement deliverable from the
[responsive CAD service contributor project](CONTRIBUTOR_PROJECTS.md).
Studio continues to use one process per CAD operation. The experiment provides
a basis for reviewing a future service, including its memory and failure costs.

## Reproduce

Use the [qualified CAD environment](STUDIO_CAD.md), Python 3.14 and Studio
0.6.0a2.dev4 or later, with its native CPU-allocation contract. The harness uses
Linux `/proc`, process descriptors and POSIX signals. Run it outside the checkout
against installed packages:

```sh
env -u PYTHONPATH /path/to/cad-env/bin/python \
  /path/to/vinkulum/ci/qualify_cad_service.py \
  /tmp/cad-service-measurement --repeats 7 --threads 2

env -u PYTHONPATH /path/to/gui-env/bin/python \
  /path/to/vinkulum/ci/qualify_cad_service.py \
  /tmp/cad-service-measurement --verify
```

The output directory must not exist. The live run writes every request,
response, worker log, source fingerprint, version, phase measurement and
process-recovery outcome. Verification checks the retained files and numerical
references without importing OCCT. It does not rerun the process lifecycle
experiment or authenticate the record's author.

## Workload and numerical checks

Each pair uses identical captured inputs, two native OCCT threads and one BLAS
thread. Fresh/resident order alternates across repetitions and cases. A separate
first-use row for each resident workload retains its initial cost. The OS file
cache is not cleared: a fresh process is not a cold machine or cold filesystem.

The three original fixtures are a 100 × 60 × 20 mm box, a 10 mm radius through
bore cut from that box, and a 1 mm all-edge fillet of the box. Density is
7,800 kg/m³. Inputs contain stable original body identities and numerical
construction parameters, not third-party CAD files.

For the box, with SI dimensions x, y, z, the independent reference is
`m = ρxyz` and `I = m/12 diag(y²+z², x²+z², x²+y²)`. For the centred bore,
subtract the cylinder's mass `ρπr²z` and its centre-of-mass inertia
`m_cylinder diag((3r²+z²)/12, (3r²+z²)/12, r²/2)`.
Mass relative tolerance is 1e-10; inertia relative tolerance is 1e-9 with
absolute tolerance 1e-15 kg·m². These references do not use OCCT or a mesh.

Every result also undergoes the normal `Body`/CAD-data validation, request hash
and native-allocation checks. Fresh and resident results agree on mass, volume,
centre, dimensions, inertia, display vertices and triangle indices. Inertia
equivalence uses relative tolerance 1e-10 and absolute 1e-15 kg·m²; vertex
tolerance is 1e-12 m. The fillet has a numerical equivalence check between the
two process modes, with no independent analytic fillet reference claimed.

## What the timings mean

Both modes use the same small instrumentation wrapper and the production worker
entry point. Fresh timings include interpreter/wrapper startup, input writing,
the request/response exchange, actual process termination, output reading and
`Body` construction. Resident timings cover the same request path while the
process stays alive; their startup, first imports and final graceful shutdown
are reported separately.
Linux process descriptors observe termination without timed-wait polling delays.

Production phase timers separate imports, request handling and geometry/capture.
The outer entry-point timer also includes output serialization. Parent-side
reading and structural validation are timed separately. Analytic comparisons
and archive verification run outside the timed request. The protocol uses JSON
files plus a small acknowledgement pipe; it does not measure Qt dispatch, VTK
rendering, application of an undo transaction or interactive frame latency.

RSS and peak RSS are Linux process counters in KiB, sampled after each operation
and before process termination. A resident worker keeps its
imported native libraries and allocator state between requests. Fresh workers
release their process memory at exit. The short sequence measures this cost;
it cannot establish the absence of a leak in a long-running CAD session.
Minima, medians, maxima, input/output sizes and raw observations are retained.
These local observations do not establish a speedup for the whole application.

## Failure and integration requirements

The prototype first requests a 1,000 mm fillet that the real geometry operation
rejects. It checks that the previous accepted result stays unchanged and that a
subsequent box succeeds in the same process. Separate injected stalls exercise
cancellation and a 100 ms deadline; an injected SIGABRT exercises replacement
after abrupt process death. Each replacement produces the expected box.
These injections test the experiment's supervisor, not an identified OCCT crash
or the production GUI controller. One recoverable operation error does not
qualify arbitrary native failure recovery.

A production proposal must preserve these contracts:

- Keep native CAD execution in a separate process, admit only the existing
  bounded numerical requests and preserve the captured request identity.
- Acquire the shared CPU lease for each operation. Start its timeout after
  admission, hold its lease through validation and release it after cancellation
  has actually settled. An idle process must hold no active CPU reservation.
- Treat idle, queued, running, validating and cancelling as distinct states.
  A live resident process alone does not identify an active operation. Preserve
  the current FIFO behaviour across CAD, meshing, statics and Pinocchio.
- Keep the native pool size fixed for a process lifetime, or qualify a separate
  reconfiguration contract. A changed allocation can require a replacement
  process; this experiment uses a fixed allocation throughout.
- Give each operation its own files and identity. Reject stale responses;
  cancellation and late completion must preserve the previous document and
  remain one explicit undoable apply transaction.
- On cancellation, deadline or process failure, stop and reap the worker before
  releasing resources. Restart lazily for the next request. Restart after an
  operation error until reuse after that error class is independently qualified.
- Bound idle lifetime and the number of retained workers. Measure the complete
  Studio memory budget, shutdown, package/interpreter changes and repeated
  failures before enabling retention by default.
- Validate the actual Qt controller and keyboard workflows at normal and high
  DPI. The headless harness supplies geometry and supervision evidence, not GUI
  responsiveness evidence.

The [shared CPU contract](ENGINE_CPU_ADMISSION.md) remains authoritative for
engine admission. Process reuse must preserve it while making the measured
startup savings available to users.
