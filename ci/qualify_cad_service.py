#!/usr/bin/env python3
"""Linux CAD worker experiment: fresh processes versus sequential process reuse.

This harness invokes the production cad_worker.main entry point for every
request. It is an experiment, not a persistent service used by Studio.
"""

import argparse
import hashlib
import json
import math
import os
import platform
import resource
import select
import signal
import statistics
import subprocess
import sys
import time
import uuid
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def memory():
    status = Path("/proc/self/status").read_text().splitlines()
    rss = next(int(line.split()[1]) for line in status if line.startswith("VmRSS:"))
    return {
        "rss_kib": rss,
        "peak_rss_kib": next(
            int(line.split()[1]) for line in status if line.startswith("VmHWM:")
        ),
    }


def wait_exit(process, timeout=5):
    """Observe actual termination without Popen's timed-wait polling delays."""
    if process.returncode is not None:
        return process.returncode
    descriptor = os.pidfd_open(process.pid)
    try:
        if not select.select([descriptor], [], [], timeout)[0]:
            raise TimeoutError("CAD experiment process did not terminate")
        return process.wait()
    finally:
        os.close(descriptor)


def worker(threads):
    # Keep acknowledgements separate from Python and native-library stdout.
    channel = os.fdopen(os.dup(sys.stdout.fileno()), "w", buffering=1)
    os.dup2(sys.stderr.fileno(), sys.stdout.fileno())
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    start = time.perf_counter()
    from vinkulum_studio import cad_worker

    channel.write(
        json.dumps(
            {
                "ready": True,
                "setup_ms": 1000 * (time.perf_counter() - start),
                **memory(),
            }
        )
        + "\n"
    )
    for line in sys.stdin:
        command = json.loads(line)
        if command.get("fault"):
            channel.write(json.dumps({"fault_entered": command["fault"]}) + "\n")
            if command["fault"] == "abort":
                os.kill(os.getpid(), signal.SIGABRT)
            # A controlled stall exercises supervisor cancellation and deadlines.
            # It is not a claim to reproduce an OCCT hang.
            time.sleep(120)
            raise RuntimeError("Supervisor did not stop the injected stall")
        sys.argv = ["cad_worker", command["input"], command["output"]]
        start = time.perf_counter()
        code = cad_worker.main()
        entrypoint_ms = 1000 * (time.perf_counter() - start)
        from OCP.OSD import OSD_ThreadPool

        pool = OSD_ThreadPool.DefaultPool_s()
        assert pool.NbThreads() == threads
        channel.write(
            json.dumps(
                {
                    "code": code,
                    "entrypoint_ms": entrypoint_ms,
                    "native_threads": pool.NbThreads(),
                    "native_default_threads": pool.NbDefaultThreadsToLaunch(),
                    **memory(),
                }
            )
            + "\n"
        )
        if command.get("exit"):
            break


class Session:
    def __init__(self, root, threads):
        root.mkdir()
        self.root, self.threads = root, threads
        self.log = (root / "worker.log").open("wb")
        env = dict(os.environ)
        for key in (
            "OPENBLAS_NUM_THREADS",
            "GOTO_NUM_THREADS",
            "MKL_NUM_THREADS",
            "BLIS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            env[key] = "1"
        env.update(
            OMP_NUM_THREADS=str(threads),
            OMP_THREAD_LIMIT=str(threads),
            OMP_MAX_ACTIVE_LEVELS="1",
            OMP_DYNAMIC="FALSE",
            PYTHONUNBUFFERED="1",
        )
        self.start = time.perf_counter()
        self.process = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "--worker", str(threads)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=self.log,
            cwd=root,
            env=env,
            bufsize=0,
        )
        self.buffer = b""
        try:
            self.ready = self.receive(30)
            assert self.ready.get("ready")
        except BaseException:
            self.close()
            raise
        self.startup_ms = 1000 * (time.perf_counter() - self.start)
        self.count, self.previous = 0, None

    def send(self, value):
        self.process.stdin.write((json.dumps(value) + "\n").encode())

    def receive(self, timeout):
        deadline = time.monotonic() + timeout
        while b"\n" not in self.buffer:
            remaining = deadline - time.monotonic()
            if (
                remaining <= 0
                or not select.select([self.process.stdout], [], [], remaining)[0]
            ):
                raise TimeoutError("CAD experiment response deadline expired")
            data = os.read(self.process.stdout.fileno(), 65536)
            if not data:
                raise EOFError("CAD experiment worker terminated")
            self.buffer += data
        line, self.buffer = self.buffer.split(b"\n", 1)
        return json.loads(line)

    def evaluate(self, request, *, fresh=False, expect_failure=False):
        from vinkulum_studio.document import Body

        self.count += 1
        source = self.root / f"{self.count:03}-input.json"
        target = self.root / f"{self.count:03}-output.json"
        request = {
            **request,
            "execution": {"threads": self.threads, "budget": self.threads},
        }
        start = time.perf_counter()
        save(source, request)
        written = time.perf_counter()
        self.send({"input": str(source), "output": str(target), "exit": fresh})
        ack = self.receive(60)
        answered = time.perf_counter()
        if fresh:
            assert wait_exit(self.process) == 0
        exited = time.perf_counter()
        response = json.loads(target.read_text())
        if expect_failure:
            assert ack["code"] == 1 and response["status"] == "failed", response
            body = None
        else:
            assert ack["code"] == 0 and response["status"] == "completed"
            assert (
                response["request_sha256"]
                == hashlib.sha256(source.read_bytes()).hexdigest()
            )
            body = Body.from_dict(response["body"])
            assert (
                ack["native_threads"] == ack["native_default_threads"] == self.threads
            )
            execution = response["execution"]
            assert execution["allocation"] == request["execution"]
            assert execution["threads"] == execution["default_threads"] == self.threads
            assert math.isclose(
                body.mass, body.cad.volume_m3 * request["density"], rel_tol=1e-12
            )
            self.previous = body
        ended = time.perf_counter()
        return body, {
            "input": str(source.relative_to(self.root.parent)),
            "output": str(target.relative_to(self.root.parent)),
            "request_write_ms": 1000 * (written - start),
            "round_trip_ms": 1000 * (answered - written),
            "exit_wait_ms": 1000 * (exited - answered),
            "response_read_validation_ms": 1000 * (ended - exited),
            "end_to_end_ms": 1000 * (ended - start) + (self.startup_ms if fresh else 0),
            "request_bytes": source.stat().st_size,
            "response_bytes": target.stat().st_size,
            "production_timings_ms": response.get("timings_ms"),
            **ack,
        }

    def close(self):
        if self.process.poll() is None:
            self.process.kill()
        wait_exit(self.process)
        self.process.stdin.close()
        self.process.stdout.close()
        self.log.close()


def fixtures():
    from vinkulum_studio.document import Body

    def part(name, shape, dimensions):
        return asdict(
            Body(
                str(uuid.uuid5(uuid.NAMESPACE_URL, "vinkulum/cad-service/" + name)),
                name,
                shape,
                dimensions,
            )
        )

    box = part("Plate", "box", (0.1, 0.06, 0.02))
    tool = part("Bore", "cylinder", (0.01, 0.04))
    return {
        "box": {
            "operation": "box",
            "name": "Plate",
            "dimensions_mm": [100, 60, 20],
            "density": 7800,
        },
        "cut": {"operation": "cut", "a": box, "b": tool, "density": 7800},
        "fillet": {"operation": "fillet", "a": box, "radius_mm": 1, "density": 7800},
    }


def equivalent(actual, expected):
    import numpy as np

    np.testing.assert_allclose(actual.mass, expected.mass, rtol=1e-11)
    np.testing.assert_allclose(actual.cad.volume_m3, expected.cad.volume_m3, rtol=1e-11)
    np.testing.assert_allclose(
        actual.position, expected.position, atol=1e-12, rtol=1e-11
    )
    np.testing.assert_allclose(
        actual.dimensions, expected.dimensions, atol=1e-12, rtol=1e-11
    )
    np.testing.assert_allclose(
        actual.inertia(), expected.inertia(), atol=1e-15, rtol=1e-10
    )
    # Same OCCT version/tolerances: also retain display tessellation equivalence.
    np.testing.assert_allclose(
        actual.cad.vertices_m, expected.cad.vertices_m, atol=1e-12, rtol=1e-11
    )
    np.testing.assert_array_equal(actual.cad.triangles, expected.cad.triangles)


def analytic(name, body):
    import numpy as np

    if name == "fillet":
        return
    x, y, z, radius, rho = 0.1, 0.06, 0.02, 0.01, 7800
    mass = x * y * z * rho
    inertia = mass / 12 * np.array([y * y + z * z, x * x + z * z, x * x + y * y])
    if name == "cut":
        removed = math.pi * radius**2 * z * rho
        mass -= removed
        inertia -= removed * np.array(
            [(3 * radius**2 + z * z) / 12] * 2 + [radius**2 / 2]
        )
    np.testing.assert_allclose(body.mass, mass, rtol=1e-10)
    np.testing.assert_allclose(
        body.inertia(), np.diag(inertia).ravel(), rtol=1e-9, atol=1e-15
    )


def faults(root, threads, requests):
    outcomes = {}
    resident = Session(root / "operation-failure", threads)
    try:
        previous, _ = resident.evaluate(requests["box"])
        bad = {**requests["fillet"], "radius_mm": 1000}
        body, failed = resident.evaluate(bad, expect_failure=True)
        assert body is None and resident.previous is previous
        recovered, _ = resident.evaluate(requests["box"])
        equivalent(recovered, previous)
        outcomes["occt_operation_error"] = {
            "same_process_recovery": True,
            "previous_preserved": True,
            "response": failed["output"],
        }
    finally:
        resident.close()
    for fault in ("cancel", "timeout", "abort"):
        resident = Session(root / fault, threads)
        try:
            previous, _ = resident.evaluate(requests["box"])
            resident.send({"fault": fault})
            assert resident.receive(5) == {"fault_entered": fault}
            start = time.perf_counter()
            if fault == "timeout":
                try:
                    resident.receive(0.1)
                except TimeoutError:
                    pass
                else:
                    raise AssertionError("Stall escaped its deadline")
            if fault == "abort":
                assert wait_exit(resident.process) == -signal.SIGABRT
            else:
                resident.process.kill()
                wait_exit(resident.process)
            assert resident.previous is previous
            outcomes[fault] = {
                "injected": True,
                "settle_ms": 1000 * (time.perf_counter() - start),
                "previous_preserved": True,
            }
        finally:
            resident.close()
        replacement = Session(root / (fault + "-replacement"), threads)
        try:
            body, row = replacement.evaluate(requests["box"])
            equivalent(body, previous)
            outcomes[fault]["replacement_ready_and_box_ms"] = (
                replacement.startup_ms + row["end_to_end_ms"]
            )
        finally:
            replacement.close()
    return outcomes


def verify_archive(root):
    """Recheck retained geometry and file identity without importing OCCT."""
    from vinkulum_studio.document import Body

    root = root.resolve()
    report = json.loads((root / "qualification.json").read_text())
    assert (
        report["format"] == "vinkulum-cad-service-experiment" and report["schema"] == 1
    )
    for name, digest in report["files_sha256"].items():
        path = (root / name).resolve()
        assert path.is_relative_to(root)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
    references, counts = {}, {}
    for row in report["rows"]:
        assert row["input"] in report["files_sha256"]
        assert row["output"] in report["files_sha256"]
        source = root / row["input"]
        request = json.loads(source.read_text())
        result = json.loads((root / row["output"]).read_text())
        assert (
            result["request_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
        )
        assert result["status"] == "completed"
        assert result["execution"]["allocation"] == request["execution"]
        body = Body.from_dict(result["body"])
        if request.get("a"):
            assert body.id == request["a"]["id"]
        analytic(row["case"], body)
        if row["mode"] == "resident_first_use":
            assert row["case"] not in references
            references[row["case"]] = body
        else:
            equivalent(body, references[row["case"]])
        key = row["mode"], row["case"]
        counts[key] = counts.get(key, 0) + 1
    assert set(references) == {"box", "cut", "fillet"}
    for name in references:
        for mode in ("fresh", "resident"):
            assert counts[mode, name] == report["repeats"]
    assert not any(
        name == "OCP" or name.startswith("OCP.") or name == "build123d"
        for name in sys.modules
    )
    print(
        f"Verified {len(report['rows'])} retained CAD results without importing OCCT."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, nargs="?")
    parser.add_argument("--worker", type=int)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=7)
    args = parser.parse_args()
    if sys.platform != "linux":
        parser.error(
            "This experiment records Linux /proc RSS and POSIX process signals"
        )
    if args.worker:
        worker(args.worker)
        return
    if args.verify:
        if args.output is None:
            parser.error("Supply the extracted archive directory")
        verify_archive(args.output)
        return
    if (
        args.output is None
        or not 1 <= args.threads <= 32
        or not 3 <= args.repeats <= 100
    ):
        parser.error("Supply an output directory, 1–32 threads and 3–100 repeats")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    requests = fixtures()
    save(root / "fixtures.json", requests)
    import vinkulum_studio

    package = Path(vinkulum_studio.__file__).parent
    report = {
        "format": "vinkulum-cad-service-experiment",
        "schema": 1,
        "platform": platform.platform(),
        "load_average_start": os.getloadavg(),
        "parent_memory_start": memory(),
        "studio_source_version": vinkulum_studio.__version__,
        "python": sys.version,
        "versions": {
            name: version(name)
            for name in (
                "vinkulum",
                "vinkulum-studio",
                "build123d",
                "cadquery-ocp-novtk",
                "numpy",
            )
        },
        "source_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(package.glob("*.py"))
        },
        "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "threads": args.threads,
        "repeats": args.repeats,
        "rows": [],
    }
    references = {}
    resident = Session(root / "resident", args.threads)
    report["resident_startup"] = {"wall_ms": resident.startup_ms, **resident.ready}
    try:
        # Retain first-use cost; steady-state rows do not hide it in an average.
        for name, request in requests.items():
            body, row = resident.evaluate(request)
            analytic(name, body)
            references[name] = body
            report["rows"].append({"mode": "resident_first_use", "case": name, **row})
        for repeat in range(args.repeats):
            for index, (name, request) in enumerate(requests.items()):
                # Alternate paired order to reduce systematic drift bias.
                order = (
                    ("fresh", "resident")
                    if (repeat + index) % 2 == 0
                    else ("resident", "fresh")
                )
                for mode in order:
                    session = (
                        resident
                        if mode == "resident"
                        else Session(root / f"fresh-{repeat}-{name}", args.threads)
                    )
                    try:
                        body, row = session.evaluate(request, fresh=mode == "fresh")
                        analytic(name, body)
                        equivalent(body, references[name])
                        report["rows"].append(
                            {"mode": mode, "case": name, "repeat": repeat, **row}
                        )
                    finally:
                        if mode == "fresh":
                            session.close()
            save(root / "partial.json", report)
        stopping = time.perf_counter()
        resident.process.stdin.close()
        assert wait_exit(resident.process) == 0
        report["resident_shutdown_ms"] = 1000 * (time.perf_counter() - stopping)
    finally:
        resident.close()
    report["faults"] = faults(root, args.threads, requests)
    report["load_average_end"] = os.getloadavg()
    report["parent_memory_end"] = memory()
    report["summary"] = {}
    for name in requests:
        report["summary"][name] = {}
        for mode in ("fresh", "resident"):
            rows = [
                row
                for row in report["rows"]
                if row["case"] == name and row["mode"] == mode
            ]
            times = [row["end_to_end_ms"] for row in rows]
            report["summary"][name][mode] = {
                "median_ms": statistics.median(times),
                "min_ms": min(times),
                "max_ms": max(times),
                "max_peak_rss_kib": max(row["peak_rss_kib"] for row in rows),
                "rss_first_kib": rows[0]["rss_kib"],
                "rss_last_kib": rows[-1]["rss_kib"],
            }
    report["files_sha256"] = {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name not in ("partial.json", "qualification.json")
    }
    save(root / "qualification.json", report)
    verify_archive(root)
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
