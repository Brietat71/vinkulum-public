# Engine allocation probes — Studio 0.6.0a2.dev4

Linux x86-64, Python 3.14.7, OCCT 8.0.1, adapted build123d 0.11.1 and
Gmsh 5.0.0-git-91b4154 built with OCCT 8 and OpenMP. The source and installed
thread-policy hashes, package versions and CPU affinity are in each JSON.
The [table](table.md) reports exploratory observations, not a performance gate.

The [recipe](../../../ci/probe_engine_threads.py) runs in a fresh process for
each engine/allocation. It requires Studio's CAD extra and `psutil` for OCCT
thread accounting. For example, on a machine allowing CPUs 0–7:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python ci/probe_engine_threads.py occt 4 /tmp/occt-new.json --affinity 0,1,2,3,4,5,6,7
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python ci/probe_engine_threads.py gmsh-hxt 4 /tmp/hxt-new.json --affinity 0,1,2,3,4,5,6,7 --gmsh /path/to/gmsh
```

OCCT cuts 81 through-holes of radius 2 mm in a 100 × 100 × 10 mm plate,
three times per process. The measured interval includes the cuts, validity
checks and volume queries, excluding imports and tool construction. All final
volumes are 89821.23980236918 mm³, against the independent formula
`100000 − 81π × 2² × 10 = 89821.23980236908 mm³`.
Per-thread CPU-time increments show 1, 2 and 4 active OS threads respectively
(activity threshold 0.01 s). Timings use one aggregate interval per allocation:
the observed ratio is not a statistically established or universal speedup.

HXT meshes a 40 × 30 × 20 mm CAD box with 5 mm target size and quadratic
tetrahedra. Times include the complete adapter and archive admission, but not
CAD construction. The ZIPs retain the original inputs, mesh, result, engine
information and logs. Extract one to a fresh directory and call
`vinkulum_studio.meshing.load_mesh(directory)` to recheck it without Gmsh.
The startup logs report the configured 1, 2 and 4 thread limits. Delaunay
insertion selected one thread on these small meshes; connectivity varies.
No 3D HXT parallel speedup is inferred from these observations.

The retained logs record 176 Studio tests with no failures (11 optional skips),
11 successful separate-environment Pinocchio tests, and keyboard CAD recipes
at both 100%/200% scaling and 420/280 px inspector widths. The 10 skipped
Pinocchio unit tests are covered by that separate environment; the keyboard
test is opt-in and runs in its own X11/Openbox sessions.
