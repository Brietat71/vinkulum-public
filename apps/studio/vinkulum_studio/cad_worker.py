"""One CAD operation per process, with bounded requests and no script evaluation."""

import hashlib
import sys
import time
from dataclasses import asdict
from pathlib import Path

from .document import MAX_PROJECT_BYTES
from .model import read_json, write_json


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: cad_worker request.json response.json")
    try:
        start = time.perf_counter()
        from .engine_threads import configure_occt_threads
        from .execution import ExecutionPlan

        request = read_json(sys.argv[1], MAX_PROJECT_BYTES)
        plan = ExecutionPlan.from_dict(request["execution"])
        runtime = configure_occt_threads(plan.threads)
        from .cad import execute

        imported = time.perf_counter()
        request_sha256 = hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()
        loaded = time.perf_counter()
        result = execute(request)
        result["execution"] = {"allocation": asdict(plan), **runtime}
        result["request_sha256"] = request_sha256
        result["timings_ms"] = {
            "import": 1000 * (imported - start),
            "request_read": 1000 * (loaded - imported),
            "geometry_and_capture": 1000 * (time.perf_counter() - loaded),
        }
        write_json(sys.argv[2], {"status": "completed", **result})
        return 0
    except Exception as error:  # noqa: BLE001 -- native failure is contained at the process boundary.
        write_json(sys.argv[2], {"status": "failed", "error": str(error)[:2048]})
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
