"""Bounded session results and identity-based comparison, independent of Qt."""

from dataclasses import fields

ARRAYS = (
    "time",
    "position",
    "rotation",
    "velocity",
    "angular_velocity",
    "joint_coordinate",
)


def result_bytes(result):
    return sum(getattr(result, name).nbytes for name in ARRAYS)


class RunArchive:
    def __init__(self, max_runs=8, max_bytes=128 * 1024 * 1024):
        self.max_runs = max_runs
        self.max_bytes = max_bytes
        self.results = []

    def add(self, result):
        if result_bytes(result) > self.max_bytes:
            raise ValueError("Result exceeds the session history budget.")
        self.results = [r for r in self.results if r.run_id != result.run_id]
        self.results.append(result)
        removed = []
        while (
            len(self.results) > self.max_runs
            or sum(map(result_bytes, self.results)) > self.max_bytes
        ):
            removed.append(self.results.pop(0).run_id)
        return removed

    def find(self, run_id):
        return next((r for r in self.results if r.run_id == run_id), None)


def series_keys(project):
    """Stable physical identities, in MechanicalResult.series() order."""
    keys = [
        (body.id, quantity, axis)
        for body in project.bodies
        for axis in range(3)
        for quantity in ("position", "velocity")
    ]
    keys += [
        (joint.id, joint.kind, 0)
        for joint in project.joints
        if joint.kind in {"pivot", "glissiere"}
    ]
    return keys


def model_differences(a, b):
    """Human-readable summary; full captured inputs remain in provenance."""
    changes = []
    for key, label in (
        ("id", "project identity"),
        ("gravity", "gravity"),
        ("duration", "duration"),
        ("step", "time step"),
    ):
        if getattr(a, key) != getattr(b, key):
            changes.append(label)
    for category, label in (
        ("bodies", "bodies"),
        ("joints", "joints"),
        ("loads", "loads"),
    ):
        left = {o.id: o for o in getattr(a, category)}
        right = {o.id: o for o in getattr(b, category)}
        if left.keys() != right.keys():
            changes.append(f"set of {label}")
        for key in left.keys() & right.keys():
            different = [
                f.name
                for f in fields(left[key])
                if getattr(left[key], f.name) != getattr(right[key], f.name)
            ]
            if different:
                names = {
                    "name": "name",
                    "mass": "mass",
                    "dimensions": "dimensions",
                    "position": "position",
                    "orientation": "orientation",
                    "inertia_mode": "inertia mode",
                    "explicit_inertia": "inertia",
                    "motion": "prescribed motion",
                    "force": "force",
                    "moment": "moment",
                    "point": "application point",
                }
                changes.append(
                    f"{left[key].name} : "
                    + ", ".join(names.get(f, f) for f in different)
                )
    return changes
