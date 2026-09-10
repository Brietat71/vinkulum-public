"""Real chain scaling: identical models, native pools, wall/CPU time, state checks.

Run in an otherwise idle environment; measurements are observations, not gates.
Example: taskset -c 0-7 python ci/measure_native_threads.py --output report.json
"""

import argparse
import hashlib
import json
import os
import platform
import statistics
import time
from pathlib import Path

import numpy as np
import vinkulum
from vinkulum import ExecutionPool
from vinkulum.test_execution import chain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bodies", type=int, nargs="+", default=[1, 64, 256, 900])
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if min(args.bodies + args.threads + [args.repeats]) < 1 or 1 not in args.threads:
        parser.error("Positive counts and a one-thread reference are required.")
    counts = sorted(set(args.threads))
    rows = []
    for count in args.bodies:
        samples = {n: [] for n in counts}
        reference = None
        # Rotate execution order to reduce systematic warm-up/drift bias.
        for repeat in range(-1, args.repeats):
            order = list(counts)
            if repeat >= 0:
                shift = repeat % len(order)
                order = order[shift:] + order[:shift]
            for threads in order:
                pool = ExecutionPool(threads)
                model = chain(count, pool)
                cpu0, wall0 = time.process_time(), time.perf_counter()
                model.simule(.02, .001, tous=10**9)
                wall, cpu = time.perf_counter() - wall0, time.process_time() - cpu0
                state = model.etat()
                if reference is None:
                    reference = state
                same_state = state == reference
                if not same_state:
                    raise RuntimeError(f"State differs for {count} bodies / {threads} threads.")
                if repeat >= 0:
                    samples[threads].append({
                        "wall_s": wall, "cpu_s": cpu, "cpu_over_wall": cpu / wall,
                        "chronos_s": model.chronos(), "stats": model.stats(),
                        "constraint_l2": float(np.linalg.norm(model.phi())),
                        "identical_to_one_thread": same_state,
                    })
                del model, pool  # No other idle pool can spin during the next sample.
        baseline = statistics.median(s["wall_s"] for s in samples[1])
        for threads, sample in samples.items():
            median = statistics.median(s["wall_s"] for s in sample)
            rows.append({"bodies": count, "unknowns": 9 * count, "threads": threads,
                         "wall_median_s": median, "speedup": baseline / median,
                         "samples": sample})
            print(f"{count:4} bodies  {threads} threads  {median:.6f}s  {baseline / median:.2f}x", flush=True)
    root = Path(__file__).resolve().parents[1]
    report = {
        "schema": "vinkulum.native-thread-measurement.1",
        "kernel_version": vinkulum.__version__, "python": platform.python_version(),
        "platform": platform.platform(), "available_cpus": os.process_cpu_count(),
        "affinity": sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else None,
        "steps": 20, "step_s": .001, "repeats": args.repeats,
        "scope": "Integration only; pool construction and Python model creation excluded. One untimed warm-up per case/thread count. Pools are destroyed between samples, with at most the requested native workers alive. Ordered state equality is empirical, not a general roundoff theorem.",
        "source_sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                          for name in ("src/lib.rs", "src/execution.rs", "src/tangent.rs", "Cargo.lock", "ci/measure_native_threads.py")},
        "measurements": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
