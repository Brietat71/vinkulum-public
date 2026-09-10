"""One CAD operation per process, with bounded requests and no script evaluation."""

import hashlib
import sys
import time
from pathlib import Path

from .document import MAX_PROJECT_BYTES
from .model import read_json, write_json


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: cad_worker request.json response.json")
    try:
        start = time.perf_counter()
        from .cad import execute

        imported = time.perf_counter()
        request = read_json(sys.argv[1], MAX_PROJECT_BYTES)
        request_sha256 = hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()
        loaded = time.perf_counter()
        result = execute(request)
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
