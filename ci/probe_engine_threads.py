"""Exploratory OCCT thread activity and HXT archive probes, outside timed CI gates."""

import argparse
import hashlib
import json
import os
import time
from importlib.metadata import version
from pathlib import Path


def occt(count):
    from vinkulum_studio.engine_threads import configure_occt_threads

    allocation = configure_occt_threads(count)
    import psutil
    from build123d import Box, Cylinder, Pos

    p = psutil.Process()
    block = Box(100, 100, 10)
    holes = [
        Pos(x, y, 0) * Cylinder(2, 20)
        for x in range(-40, 41, 10)
        for y in range(-40, 41, 10)
    ]
    before = {t.id: t.user_time + t.system_time for t in p.threads()}
    start = time.perf_counter()
    volumes = []
    for _ in range(3):
        shape = block.cut(*holes)
        if not shape.is_valid or len(shape.solids()) != 1:
            raise ValueError("Boolean returned an invalid solid")
        volumes.append(shape.volume)
    elapsed = time.perf_counter() - start
    cpu = {
        str(t.id): t.user_time + t.system_time - before.get(t.id, 0)
        for t in p.threads()
    }
    expected = 100000 - 81 * 3.141592653589793 * 4 * 10
    if any(abs(v - expected) > 1e-8 * expected for v in volumes):
        raise ValueError("Perforated plate volume disagrees with analytic reference")
    return {
        "engine": "occt",
        "allocation": allocation,
        "elapsed_s": elapsed,
        "iterations": 3,
        "cpu_by_os_thread_s": cpu,
        "active_os_threads": sum(v > 0.01 for v in cpu.values()),
        "volume_mm3": volumes,
        "expected_volume_mm3": expected,
    }


def mesh(count, output, gmsh):
    from vinkulum_studio.engine_threads import configure_occt_threads

    configure_occt_threads(1)
    from vinkulum_studio.cad import execute
    from vinkulum_studio.document import Body
    from vinkulum_studio.meshing import HxtMeshRequest, load_mesh, run_mesh

    body = Body.from_dict(
        execute(
            {
                "operation": "box",
                "name": "HXT reference",
                "dimensions_mm": [40, 30, 20],
                "density": 7800,
            }
        )["body"]
    )
    start = time.perf_counter()
    result = run_mesh(HxtMeshRequest(body, 5.0, 2, count), output, executable=gmsh)
    elapsed = time.perf_counter() - start
    assert load_mesh(output)[0] == result
    return {
        "engine": "gmsh-hxt",
        "threads": count,
        "elapsed_s": elapsed,
        "nodes": len(result.mesh.nodes),
        "elements": len(result.mesh.elements),
        "archive": str(output),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("engine", choices=["occt", "gmsh-hxt"])
    parser.add_argument("threads", type=int)
    parser.add_argument("output", type=Path)
    parser.add_argument("--gmsh", default="gmsh")
    parser.add_argument(
        "--affinity",
        help="Comma-separated Linux CPU IDs; otherwise retain caller affinity",
    )
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix("").exists():
        raise SystemExit("Choose a new output path; prior evidence is preserved")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.affinity:
        os.sched_setaffinity(0, [int(v) for v in args.affinity.split(",")])
    report = (
        occt(args.threads)
        if args.engine == "occt"
        else mesh(args.threads, args.output.with_suffix(""), args.gmsh)
    )
    import vinkulum_studio.engine_threads as policy

    report["provenance"] = {
        "studio": version("vinkulum-studio"),
        "occt": version("cadquery-ocp-novtk"),
        "build123d": version("build123d"),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "thread_policy_sha256": hashlib.sha256(
            Path(policy.__file__).read_bytes()
        ).hexdigest(),
        "affinity": sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else None,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)
