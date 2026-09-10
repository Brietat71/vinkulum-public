"""Optional CalculiX linear statics adapter; SI, affine C3D8 meshes only.

CalculiX is an external executable. This module does not import its libraries,
modify the multibody document or claim a bound on the finite-element error.
Node IDs are one-based; constrained axes are 1=X, 2=Y, 3=Z in the world frame.
"""

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .model import finite_number as finite
from .model import read_json, write_json

MAX_NODES = 6000
MAX_ELEMENTS = 5000
MAX_OUTPUT_BYTES = 32 * 1024 * 1024
FACES = (
    (0, 1, 2, 3),
    (4, 5, 6, 7),
    (0, 1, 5, 4),
    (1, 2, 6, 5),
    (2, 3, 7, 6),
    (3, 0, 4, 7),
)


@dataclass(frozen=True)
class StaticStudy:
    nodes: tuple
    elements: tuple
    fixed_dofs: tuple
    forces: tuple
    young_pa: float
    poisson: float

    def __post_init__(self):
        for key in ("nodes", "elements", "fixed_dofs", "forces"):
            rows = getattr(self, key)
            if not isinstance(rows, (list, tuple)) or not all(
                isinstance(r, (list, tuple)) for r in rows
            ):
                raise ValueError(f"{key}: expected rows of numbers.")
            object.__setattr__(self, key, tuple(tuple(r) for r in rows))
        if not finite(self.young_pa) or self.young_pa <= 0:
            raise ValueError("Young modulus must be positive and finite, in Pa.")
        if not finite(self.poisson) or not -1 < self.poisson < 0.5:
            raise ValueError("Isotropic elasticity requires -1 < Poisson ratio < 0.5.")
        if (
            not 8 <= len(self.nodes) <= MAX_NODES
            or not 1 <= len(self.elements) <= MAX_ELEMENTS
        ):
            raise ValueError("Static mesh budget exceeded or mesh empty.")
        if any(len(r) != 3 or not all(finite(v) for v in r) for r in self.nodes):
            raise ValueError(
                "Nodes must have three finite world coordinates in metres."
            )
        n = len(self.nodes)
        valid_node = lambda i: type(i) is int and 1 <= i <= n
        if any(
            len(e) != 8 or len(set(e)) != 8 or not all(valid_node(i) for i in e)
            for e in self.elements
        ):
            raise ValueError("Each C3D8 element needs eight distinct valid node IDs.")
        if len({tuple(sorted(e)) for e in self.elements}) != len(self.elements):
            raise ValueError("Duplicate elements are not allowed.")
        if {i for e in self.elements for i in e} != set(range(1, n + 1)):
            raise ValueError("Every node must belong to an element.")
        if not self.fixed_dofs or len(set(self.fixed_dofs)) != len(self.fixed_dofs):
            raise ValueError("Supports must be nonempty and contain no duplicate DOFs.")
        if any(
            len(r) != 2
            or not valid_node(r[0])
            or type(r[1]) is not int
            or r[1] not in (1, 2, 3)
            for r in self.fixed_dofs
        ):
            raise ValueError("Supports require a valid node ID and axis 1, 2 or 3.")
        if any(
            len(r) != 4 or not valid_node(r[0]) or not all(finite(v) for v in r[1:])
            for r in self.forces
        ):
            raise ValueError(
                "Loads require a node ID and three finite world forces in N."
            )
        if len({r[0] for r in self.forces}) != len(self.forces):
            raise ValueError("Combine forces applied to the same node before solving.")
        points = np.array(self.nodes, dtype=float)
        faces, neighbours = {}, [set() for _ in self.elements]
        for index, element in enumerate(self.elements):
            p = points[np.array(element) - 1]
            a, b, c = p[1] - p[0], p[3] - p[0], p[4] - p[0]
            scale = max(np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c))
            expected = np.array(
                (
                    (0, 0, 0),
                    (1, 0, 0),
                    (1, 1, 0),
                    (0, 1, 0),
                    (0, 0, 1),
                    (1, 0, 1),
                    (1, 1, 1),
                    (0, 1, 1),
                )
            )
            jac = np.column_stack((a, b, c)) / scale if scale > 0 else np.zeros((3, 3))
            if not np.isfinite(jac).all() or np.linalg.det(jac) <= 1e-12:
                raise ValueError(
                    "An element is inverted, degenerate or too poorly scaled."
                )
            # Measure distortion in all element directions. A tolerance based
            # on the longest edge alone can admit inversion of a thin element.
            local_points = np.linalg.solve(jac, ((p - p[0]) / scale).T).T
            if (
                not np.isfinite(local_points).all()
                or np.max(abs(local_points - expected)) > 1e-10
            ):
                raise ValueError(
                    "The first adapter supports affine C3D8 elements only."
                )
            for face in FACES:
                key = tuple(sorted(element[i] for i in face))
                owners = faces.setdefault(key, [])
                owners.append(index)
                if len(owners) > 2:
                    raise ValueError("A mesh face has more than two owners.")
                if len(owners) == 2:
                    left, right = owners
                    neighbours[left].add(right)
                    neighbours[right].add(left)
        seen, pending = set(), [0]
        while pending:
            index = pending.pop()
            if index not in seen:
                seen.add(index)
                pending.extend(neighbours[index] - seen)
        if len(seen) != len(self.elements):
            raise ValueError(
                "Elements must form one connected mesh through shared faces."
            )
        # Check removal of the six global rigid modes with dimensionless rows.
        centred = (points - points.mean(axis=0)) / np.ptp(points, axis=0).max()
        rows = []
        for node, axis in self.fixed_dofs:
            x, y, z = centred[node - 1]
            modes = np.column_stack((np.eye(3), ((0, z, -y), (-z, 0, x), (y, -x, 0))))
            rows.append(modes[axis - 1])
        if np.linalg.matrix_rank(np.array(rows), tol=1e-10) != 6:
            raise ValueError(
                "Supports leave a rigid-body mode or are numerically ambiguous."
            )

    def input_deck(self):
        lines = ["*HEADING", "Vinkulum linear statics, SI: m N Pa", "*NODE,NSET=ALLN"]
        lines += [
            f"{i}," + ",".join(format(float(v), ".17g") for v in p)
            for i, p in enumerate(self.nodes, 1)
        ]
        lines += ["*ELEMENT,TYPE=C3D8,ELSET=ALLE"]
        lines += [
            f"{i}," + ",".join(map(str, e)) for i, e in enumerate(self.elements, 1)
        ]
        lines += [
            "*MATERIAL,NAME=ISOTROPIC",
            "*ELASTIC",
            f"{self.young_pa:.17g},{self.poisson:.17g}",
            "*SOLID SECTION,ELSET=ALLE,MATERIAL=ISOTROPIC",
            "*BOUNDARY",
        ]
        lines += [f"{node},{axis},{axis},0" for node, axis in self.fixed_dofs]
        lines += ["*STEP", "*STATIC"]
        if any(value != 0 for row in self.forces for value in row[1:]):
            lines += ["*CLOAD"]
            lines += [
                f"{row[0]},{axis},{value:.17g}"
                for row in self.forces
                for axis, value in enumerate(row[1:], 1)
                if value != 0
            ]
        lines += [
            "*NODE PRINT,NSET=ALLN",
            "U,RF",
            "*EL PRINT,ELSET=ALLE",
            "S,ENER",
            "*END STEP",
        ]
        return "\n".join(lines) + "\n"


def load_study(path):
    data = read_json(path, max_bytes=4 * 1024 * 1024)
    if (
        not isinstance(data, dict)
        or set(data) != {"format", "schema", "study"}
        or data["format"] != "vinkulum-static-study"
        or type(data["schema"]) is not int
        or data["schema"] != 1
        or not isinstance(data["study"], dict)
    ):
        raise ValueError("Unsupported or incomplete static-study document.")
    try:
        return StaticStudy(**data["study"])
    except TypeError as exc:
        raise ValueError("Invalid or unknown static-study fields.") from exc


def parse_dat(text, study):
    """Accept exactly one complete final state; preserve integration-point stress."""
    tables = {}
    specifications = (
        ("displacements", "displacements ", "ALLN", len(study.nodes), 4),
        ("external_forces", "forces ", "ALLN", len(study.nodes), 4),
        ("stress", "stresses ", "ALLE", len(study.elements) * 8, 8),
        (
            "energy_density",
            "internal energy density ",
            "ALLE",
            len(study.elements) * 8,
            3,
        ),
    )
    for name, label, node_set, row_count, columns in specifications:
        pattern = re.compile(
            r"^ "
            + re.escape(label)
            + r"[^\n]* for set "
            + node_set
            + r" and time\s+(\S+)\s*\n\s*\n([^\n]+(?:\n[^\n]+)*)",
            re.MULTILINE,
        )
        matches = list(pattern.finditer(text))
        if len(matches) != 1:
            raise ValueError(f"Expected one {name} table in the final static state.")
        match = matches[0]
        if float(match[1].replace("D", "E")) != 1.0:
            raise ValueError("Incomplete static load step.")
        try:
            rows = np.array(
                [
                    [float(v.replace("D", "E")) for v in line.split()]
                    for line in match[2].splitlines()
                ]
            )
        except ValueError as exc:
            raise ValueError(f"Malformed {name} table.") from exc
        if rows.shape != (row_count, columns) or not np.isfinite(rows).all():
            raise ValueError(f"Incomplete or nonfinite {name} table.")
        identifiers = rows[:, :1] if columns == 4 else rows[:, :2]
        expected = (
            np.arange(1, row_count + 1)[:, None]
            if columns == 4
            else np.array(
                [(e, p) for e in range(1, len(study.elements) + 1) for p in range(1, 9)]
            )
        )
        if not np.array_equal(identifiers, expected):
            raise ValueError(f"Unexpected or duplicate identities in {name}.")
        tables[name] = rows[:, identifiers.shape[1] :].copy()
    applied = np.zeros_like(tables["external_forces"])
    for node, x, y, z in study.forces:
        applied[node - 1] = (x, y, z)
    # RF is the sum of reactions AND applied nodal/distributed loads in CalculiX.
    tables["reactions"] = tables["external_forces"] - applied
    fixed = np.zeros_like(applied, dtype=bool)
    for node, axis in study.fixed_dofs:
        fixed[node - 1, axis - 1] = True
    if np.any(tables["displacements"][fixed] != 0):
        raise ValueError("Calculated displacements violate a zero support.")
    force_scale = np.abs(applied).sum() + np.abs(tables["external_forces"]).sum()
    if not math.isfinite(force_scale):
        raise ValueError("Calculated force scale is not finite.")
    force_tolerance = max(2e-6 * force_scale, np.finfo(float).tiny)
    if np.any(abs(tables["reactions"][~fixed]) > force_tolerance):
        raise ValueError("Calculated forces contain an unsupported reaction.")
    if np.any(abs(tables["external_forces"].sum(axis=0)) > force_tolerance):
        raise ValueError("Calculated external forces do not balance.")
    points = np.array(study.nodes, dtype=float)
    # Use undeformed coordinates for linear statics, centred and scaled to
    # avoid making the moment tolerance depend on the world-frame origin.
    levers = (points - points.mean(axis=0)) / np.ptp(points, axis=0).max()
    moment = np.cross(levers, tables["external_forces"]).sum(axis=0)
    moment_tolerance = force_tolerance * np.linalg.norm(levers, axis=1).max()
    if not np.isfinite(moment).all() or np.any(abs(moment) > moment_tolerance):
        raise ValueError("Calculated external moments do not balance.")
    if np.any(tables["energy_density"] < 0):
        raise ValueError("Negative elastic energy density.")
    volumes = []
    for element in study.elements:
        p = points[np.array(element) - 1]
        volumes.append(
            np.linalg.det(np.column_stack((p[1] - p[0], p[3] - p[0], p[4] - p[0])))
        )
    internal = float(
        np.dot(volumes, tables["energy_density"].reshape(-1, 8).mean(axis=1))
    )
    work = float(0.5 * np.sum(applied * tables["displacements"]))
    energy_scale = max(abs(internal), abs(work), np.finfo(float).tiny)
    if (
        not math.isfinite(internal)
        or not math.isfinite(work)
        or abs(internal - work) > 5e-6 * energy_scale
    ):
        raise ValueError("Elastic energy and external work are inconsistent.")
    tables["strain_energy_J"] = np.array(internal)
    return {name: value.tolist() for name, value in tables.items()}


def run_static(study, directory, *, executable="ccx", timeout=60):
    """Run in a new directory. Keep the input, raw output and diagnostics on failure."""
    if not isinstance(study, StaticStudy):
        raise TypeError("A validated StaticStudy is required.")
    if not finite(timeout) or not 0 < timeout <= 300:
        raise ValueError("Timeout must be between 0 and 300 seconds.")
    found = shutil.which(str(executable))
    if not found:
        raise FileNotFoundError(
            "CalculiX executable unavailable; install ccx separately."
        )
    engine = Path(found).resolve()
    engine_hash = hashlib.sha256(engine.read_bytes()).hexdigest()
    version_run = subprocess.run(
        [str(engine), "-v"], capture_output=True, text=True, timeout=10, check=False
    )
    version = re.search(r"Version\s+([0-9.]+)", version_run.stdout)
    if not version:
        raise ValueError(
            "The executable did not identify itself as a supported CalculiX CLI."
        )
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    deck = study.input_deck()
    (root / "study.inp").write_text(deck, encoding="ascii")
    write_json(
        root / "study.json",
        {"format": "vinkulum-static-study", "schema": 1, "study": asdict(study)},
    )
    with (root / "solver.log").open("wb") as log:
        result = subprocess.run(
            [str(engine), "-i", "study"],
            cwd=root,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
            env={
                **os.environ,
                "OMP_NUM_THREADS": "1",
                "CCX_NPROC_RESULTS": "1",
                "CCX_NPROC_EQUATION_SOLVER": "1",
            },
        )
    log_path, dat_path = root / "solver.log", root / "study.dat"
    if log_path.stat().st_size > MAX_OUTPUT_BYTES:
        raise ValueError("CalculiX log exceeded the output budget.")
    log_text = log_path.read_text(errors="replace")
    if (
        result.returncode != 0
        or "Job finished" not in log_text
        or "*ERROR" in log_text.upper()
    ):
        raise RuntimeError(f"CalculiX did not complete; inspect {log_path}.")
    if not dat_path.is_file() or dat_path.stat().st_size > MAX_OUTPUT_BYTES:
        raise ValueError("Missing or oversized CalculiX result.")
    raw = dat_path.read_text(encoding="ascii")
    values = parse_dat(raw, study)
    if hashlib.sha256(engine.read_bytes()).hexdigest() != engine_hash:
        raise RuntimeError("CalculiX executable changed during the calculation.")
    report = {
        "format": "vinkulum-calculix-static",
        "schema": 1,
        "engine": "CalculiX",
        "engine_version": version[1],
        "adapter_version": "0.1.0",
        "engine_sha256": engine_hash,
        "input_sha256": hashlib.sha256(deck.encode("ascii")).hexdigest(),
        "dat_sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "units": {
            "displacements": "m",
            "external_forces": "N",
            "reactions": "N",
            "stress": "Pa",
            "energy_density": "J/m3",
            "strain_energy_J": "J",
        },
        "frame": "world",
        "stress_components": ["xx", "yy", "zz", "xy", "xz", "yz"],
        "stress_location": "element integration points 1..8, element-major order",
        "load_factor": 1.0,
        "requested_threads": 1,
        "scientific_status": "NotAssessed",
        "limitations": "Linear isotropic elasticity, affine C3D8, zero supports, nodal loads. No general mesh or discretisation-error certification. ASCII output is rounded.",
        **values,
    }
    write_json(root / "result.json", report)
    return report


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path, help="Static-study JSON document")
    parser.add_argument(
        "--output", type=Path, required=True, help="New calculation directory"
    )
    parser.add_argument("--ccx", default="ccx", help="Installed CalculiX executable")
    parser.add_argument(
        "--timeout",
        type=float,
        default=60,
        help="Solver timeout in seconds; version probe has a separate 10 s limit",
    )
    args = parser.parse_args()
    try:
        result = run_static(
            load_study(args.study),
            args.output,
            executable=args.ccx,
            timeout=args.timeout,
        )
    except (
        ValueError,
        TypeError,
        OSError,
        RuntimeError,
        subprocess.SubprocessError,
    ) as exc:
        parser.exit(2, f"CalculiX study failed: {exc}\n")
    print(
        json.dumps(
            {
                "status": "completed",
                "engine_version": result["engine_version"],
                "result": str(args.output / "result.json"),
                "strain_energy_J": result["strain_energy_J"],
                "scientific_status": result["scientific_status"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
