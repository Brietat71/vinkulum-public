"""One immutable input, one fresh native solver, one atomic result file."""

from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import math
from pathlib import Path
import platform
import sys

from . import __version__
from .model import (
    G,
    INERTIA,
    RHO,
    PARAMETER_UNITS,
    UNITS,
    Parameters,
    read_json,
    write_json,
)


def simulate(parameters, run_id):
    # Import the native library exclusively in the compute process.
    import vinkulum
    import vinkulum._vinkulum as native

    p = parameters
    angle = math.radians(p.angle_deg)
    x, z = p.length * math.sin(angle), -p.length * math.cos(angle)
    solver = vinkulum.Noyau(g=[0.0, 0.0, -G])
    body = solver.corps(
        "masse",
        p.mass,
        [INERTIA, 0.0, 0.0, 0.0, INERTIA, 0.0, 0.0, 0.0, INERTIA],
        [x, 0.0, z],
    )
    solver.liaison("rotule", None, body, bloque_t=[0, 1, 2], bloque_r=[])
    trajectory = solver.simule(p.duration, p.step, rho=RHO)
    samples = [[0.0, x, z, angle]]
    for frame in trajectory:
        t, positions = frame[:2]
        x, y, z = positions[0]
        if abs(y) > 1e-7 * max(1, p.length):
            raise ValueError("Le mouvement sort du plan du cas G0.")
        samples.append([t, x, z, math.atan2(x, -z)])
    return {
        "schema_version": 1,
        "run_id": run_id,
        "parameters": asdict(p),
        "status": "completed",
        "samples": samples,
        "manifest": {
            "app_version": __version__,
            "kernel_version": vinkulum.__version__,
            "kernel_sha256": hashlib.sha256(
                Path(native.__file__).read_bytes()
            ).hexdigest(),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "machine": platform.machine(),
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "method": "generalized-alpha / SO(3)",
            "rho_infinity": RHO,
            "gravity_m_s2": G,
            "body_inertia_kg_m2": INERTIA,
            "initial_angular_velocity_rad_s": 0.0,
            "threads": 2,
            "units": UNITS,
            "parameter_units": PARAMETER_UNITS,
            "sampling": "native solver steps plus explicit initial state; no interpolation",
            "scientific_status": "NotAssessed",
            "limitations": "No trajectory error bound; no checkpoint; G0 planar pendulum only.",
        },
    }


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: python -m vinkulum_studio.worker input.json output.json"
        )
    try:
        from .document import MAX_PROJECT_BYTES

        data = read_json(sys.argv[1], MAX_PROJECT_BYTES)
        if data.get("kind") == "mechanism":
            from .document import Project
            from .mechanism import simulate_project

            result = simulate_project(
                Project.from_dict(data["parameters"]),
                data["run_id"],
                Path(sys.argv[2]).parent,
            )
        else:
            result = simulate(Parameters.from_dict(data["parameters"]), data["run_id"])
        write_json(sys.argv[2], result)
    except Exception as exc:
        write_json(
            sys.argv[2],
            {
                "status": "failed",
                "error_type": type(exc).__name__,
                "message": str(exc)[:2000],
            },
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
