"""Immutable CPU allocation recorded with each native calculation (no Qt)."""

import os
from dataclasses import asdict, dataclass


def available_cpus():
    """Affinity-aware estimate; an explicit budget can further constrain it."""
    return max(1, os.process_cpu_count() or 1)


def default_budget():
    available = available_cpus()
    raw = os.environ.get("VINKULUM_STUDIO_CPUS")
    if raw is None:
        # Leave room for the event loop on workstations, without pinning threads.
        return max(1, available - 1)
    try:
        value = int(raw)
    except ValueError as error:
        raise ValueError("VINKULUM_STUDIO_CPUS must be a positive integer.") from error
    if not 1 <= value <= available:
        raise ValueError(f"VINKULUM_STUDIO_CPUS must be between 1 and {available}.")
    return value


@dataclass(frozen=True)
class ExecutionPlan:
    threads: int
    budget: int

    def __post_init__(self):
        if any(type(v) is not int for v in (self.threads, self.budget)) or not (
            1 <= self.threads <= self.budget
        ):
            raise ValueError("Native execution requires 1 <= threads <= CPU budget.")

    @classmethod
    def from_dict(cls, data):
        if not isinstance(data, dict) or set(data) != {"threads", "budget"}:
            raise ValueError("Invalid native CPU allocation.")
        return cls(**data)

    def executor(self):
        from vinkulum import ExecutionPool

        if self.threads > available_cpus():
            raise ValueError("CPU allocation exceeds this worker's available CPUs.")
        return ExecutionPool(self.threads)

    def manifest(self, solver):
        if solver.execution_threads != self.threads or solver.execution_pool_id <= 0:
            raise ValueError("Native executor does not match its captured allocation.")
        return asdict(self) | {
            "backend": "rayon",
            "pool_id": solver.execution_pool_id,
            "worker_available_cpus": available_cpus(),
        }

    def check_manifest(self, data):
        if not isinstance(data, dict) or any(
            data.get(key) != value or type(data.get(key)) is not int
            for key, value in asdict(self).items()
        ):
            raise ValueError("Result CPU allocation differs from the captured request.")
        if (
            data.get("backend") != "rayon"
            or type(data.get("pool_id")) is not int
            or data["pool_id"] <= 0
            or type(data.get("worker_available_cpus")) is not int
            or data["worker_available_cpus"] < self.threads
        ):
            raise ValueError("Invalid native executor provenance.")
