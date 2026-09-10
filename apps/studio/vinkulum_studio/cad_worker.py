"""Bounded CAD requests, as a one-shot worker or a supervised sequential service."""

import hashlib
import json
import os
import sys
import time
import uuid
from dataclasses import asdict
from pathlib import Path

from .document import MAX_PROJECT_BYTES
from .model import read_json, write_json


def evaluate(request_path, response_path, *, request_id=None, native_threads=None):
    try:
        start = time.perf_counter()
        from .engine_threads import configure_occt_threads
        from .execution import ExecutionPlan

        request = read_json(request_path, MAX_PROJECT_BYTES)
        plan = ExecutionPlan.from_dict(request["execution"])
        if native_threads is not None and plan.threads != native_threads:
            raise ValueError(
                "A resident CAD process retains its initial native pool size."
            )
        runtime = configure_occt_threads(plan.threads)
        from .cad import execute

        imported = time.perf_counter()
        request_sha256 = hashlib.sha256(Path(request_path).read_bytes()).hexdigest()
        loaded = time.perf_counter()
        result = execute(request)
        result["execution"] = {"allocation": asdict(plan), **runtime}
        result["request_sha256"] = request_sha256
        if request_id is not None:
            result["service_request_id"] = request_id
        result["timings_ms"] = {
            "import": 1000 * (imported - start),
            "request_read": 1000 * (loaded - imported),
            "geometry_and_capture": 1000 * (time.perf_counter() - loaded),
        }
        write_json(response_path, {"status": "completed", **result})
        return 0
    except (
        Exception
    ) as error:  # noqa: BLE001 -- native failure is contained at the process boundary.
        write_json(response_path, {"status": "failed", "error": str(error)[:2048]})
        return 1


def service():
    # Python and native diagnostic output goes to stderr; stdout is only the
    # small acknowledgement protocol. The full result stays in its own file.
    channel = os.fdopen(os.dup(sys.stdout.fileno()), "w", buffering=1)
    os.dup2(sys.stderr.fileno(), sys.stdout.fileno())
    native_threads = None
    while line := sys.stdin.buffer.readline(16385):
        if len(line) > 16384 or not line.endswith(b"\n"):
            raise ValueError("Oversized or incomplete CAD service command.")
        command = json.loads(line)
        if not isinstance(command, dict) or set(command) != {
            "id",
            "input",
            "output",
            "threads",
        }:
            raise ValueError("Invalid CAD service command.")
        identifier = command["id"]
        if not isinstance(identifier, str) or uuid.UUID(identifier).hex != identifier:
            raise ValueError("Invalid CAD transaction identity.")
        if type(command["threads"]) is not int or command["threads"] < 1:
            raise ValueError("Invalid CAD native pool size.")
        if any(not isinstance(command[key], str) for key in ("input", "output")):
            raise ValueError("CAD service file paths are required.")
        if Path(command["output"]).exists():
            raise ValueError("A CAD transaction cannot overwrite an existing response.")
        if native_threads is None:
            native_threads = command["threads"]
        if command["threads"] != native_threads:
            raise ValueError("A resident CAD process cannot change native pool size.")
        code = evaluate(
            command["input"],
            command["output"],
            request_id=identifier,
            native_threads=native_threads,
        )
        channel.write(
            json.dumps({"protocol": 1, "id": identifier, "code": code}) + "\n"
        )
        if code:
            return code  # Retire after every operation error, including native errors.
    return 0


def main():
    if sys.argv[1:] == ["--service"]:
        return service()
    if len(sys.argv) != 3:
        raise SystemExit("Usage: cad_worker request.json response.json, or --service")
    return evaluate(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    raise SystemExit(main())
