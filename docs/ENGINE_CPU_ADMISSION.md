# Studio engine CPU admission

Studio `0.6.0a2.dev4` extends the application's FIFO CPU admission to CAD,
meshing, CalculiX and Pinocchio, alongside the native dynamics worker introduced
in [0.20.0](NATIVE_THREADING.md). An immutable allocation is captured when a job
is submitted. Waiting does not start an engine or consume its execution timeout.
Cancellation, failed startup and completion release the allocation; background
result validation retains it until the reader has stopped.

| Engine | Default allocation | Native execution |
| --- | --- | --- |
| Vinkulum dynamics | Application budget; adjustable in the Run menu | Explicit Rayon `ExecutionPool` |
| OCCT 8 / build123d | Up to 4 CPUs within the application budget | OCCT's native `OSD_ThreadPool`, configured before importing build123d |
| Gmsh | Up to 4 CPUs; adjustable in Execution settings | HXT with an OpenMP build, explicit `-nt` and dimensional thread limits |
| CalculiX | Up to 4 CPUs; adjustable in Execution settings | Engine phase settings, with reported maxima retained in the result |
| Pinocchio | 1 CPU | Sequential recursive operators for one captured articulated state |

The application budget still defaults to one fewer than the affinity-aware
CPU availability, with a minimum of one. For example,
`VINKULUM_STUDIO_CPUS=4 vinkulum-studio` sets a four-CPU admission budget.
Several small jobs can run together if their allocations fit. A job requiring
the whole budget waits for the earlier jobs to finish. Each Studio process has
its own scheduler. CLI calculations run independently of this GUI scheduler.

This bounds **admitted engine allocations**, not every operating-system thread
in Studio. Rendering, external applications, container quotas and numerical
libraries loaded into the GUI process are outside that accounting. Native
engines can also choose sequential algorithms for small problems. A requested
thread count is neither a measurement of CPU use nor a speedup guarantee.

## Engine-specific controls

The CAD worker selects OCCT's own thread implementation and configures its
default pool before geometry construction. The response identifies the pool
size and captured allocation; the controller rejects inconsistent responses.
Parallel Boolean operations supplied by the adapted build123d package can use
this pool. Serial geometry operations remain serial.

New GUI meshes default to HXT (`Mesh.Algorithm3D = 10`). Gmsh must report both
`Hxt` and `OpenMP` in its build options, in addition to OCCT 8 or newer.
Legacy Delaunay remains selectable. The new `HxtMeshRequest` has an explicit
algorithm field; historical `MeshRequest` payloads and generated scripts keep
their exact structure, preserving archived mesh fingerprints.

CalculiX receives `NUMBER_OF_CPUS`, `CCX_NPROC_STIFFNESS`, `CCX_NPROC_RESULTS`,
`CCX_NPROC_EQUATION_SOLVER` and the OpenMP allocation. Nested OpenMP teams are
disabled and BLAS libraries in external workers are restricted to one thread.
The CLI accepts `--threads N`. Result schema 4 stores `execution.json` and the
solver's reported per-phase maxima. Reading the archive reconstructs that
metadata from the captured allocation and raw solver log, with strict integer
validation. Older result schemas 1–3 remain readable without this extra file.
The log records phase maxima; it does not attest actual simultaneous CPU use.

Pinocchio's current single-state analysis consumes one scheduler slot. Increasing
an OpenMP setting would not parallelise that recursive algorithm. Batch
evaluation of independent states is a separate future increment.

## GUI lifecycle and verification

New CalculiX and Pinocchio result readers use receiver-bound Qt slots and worker
threads. CalculiX and mesh reports are written to a staged file off the GUI
thread, then published by rename after final cancellation checks. Closing the
statics or articulated-analysis window during validation defers destruction
until the reader finishes. The previous accepted result survives cancellation
and validation failure.

Studio `0.6.0a2.dev5` also moves CAD response admission, attachment-preserving
document preparation and interactive CalculiX/Pinocchio archive reads off the
GUI thread. Archive reads request one slot from the same scheduler; CAD checks
retain their operation's existing allocation. See the
[lifecycle contract and tests](STUDIO_BACKGROUND_ADMISSION.md).
CAD request serialization, ordinary project file I/O, STEP export writing and
VTK scene updates still have GUI-thread work. GUI numerical backends are not yet
fully accounted for by the scheduler.

[The engine tests](../apps/studio/tests/test_engine_threads.py) exercise real
queued external processes, cancellation before launch, failure cleanup,
off-GUI readers, deferred window destruction and strict archive metadata.
An actual CalculiX traction problem is compared at 1, 2 and 4 threads for
displacements, stresses, reactions and energy. Performance ratios are not CI
pass/fail criteria.

[Retained Linux probes](bancs/engine-threads-dev4/README.md) of an 81-hole Boolean plate recorded 1, 2 and 4 active
OCCT threads, with the same volume matching the analytic reference. Small HXT
box meshes reopened successfully at all three allocations; their Delaunay
insertion phases selected one thread and their connectivity could differ.
These are bounded experiments, not general solver performance qualifications.
