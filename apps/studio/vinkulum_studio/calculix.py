"""Optional CalculiX linear statics adapter; SI, C3D4/C3D10 and affine C3D8.

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

from .finite_elements import (
    NODE_COUNTS,
    POINT_COUNTS,
    TET_EDGES,
    face_nodes,
    integration_weights,
)
from .model import finite_number as finite
from .model import read_json, write_json

MAX_NODES = 6000
MAX_ELEMENTS = 5000
MAX_OUTPUT_BYTES = 32 * 1024 * 1024


@dataclass(frozen=True)
class StaticStudy:
    nodes: tuple
    elements: tuple
    fixed_dofs: tuple
    forces: tuple
    young_pa: float
    poisson: float
    element_type: str = "C3D8"

    def __post_init__(self):
        for key in ("nodes", "elements", "fixed_dofs", "forces"):
            rows = getattr(self, key)
            if not isinstance(rows, (list, tuple)) or not all(
                isinstance(r, (list, tuple)) for r in rows
            ):
                raise ValueError(f"{key}: expected rows of numbers.")
            object.__setattr__(self, key, tuple(tuple(r) for r in rows))
        if (
            not isinstance(self.element_type, str)
            or self.element_type not in NODE_COUNTS
        ):
            raise ValueError("Supported solid elements are C3D4, C3D8 and C3D10.")
        count = NODE_COUNTS[self.element_type]
        if not finite(self.young_pa) or self.young_pa <= 0:
            raise ValueError("Young modulus must be positive and finite, in Pa.")
        if not finite(self.poisson) or not -1 < self.poisson < 0.5:
            raise ValueError("Isotropic elasticity requires -1 < Poisson ratio < 0.5.")
        if (
            not count <= len(self.nodes) <= MAX_NODES
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
            len(e) != count or len(set(e)) != count or not all(valid_node(i) for i in e)
            for e in self.elements
        ):
            raise ValueError(
                f"Each {self.element_type} element needs {count} distinct valid node IDs."
            )
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
        if self.element_type == "C3D10":
            corner_ids = {i for e in self.elements for i in e[:4]}
            mid_ids = {i for e in self.elements for i in e[4:]}
            if corner_ids & mid_ids:
                raise ValueError(
                    "A quadratic node cannot be both a corner and an edge node."
                )
            edges, midpoints = {}, {}
            for element in self.elements:
                for mid, (a, b) in zip(element[4:], TET_EDGES):
                    key = tuple(sorted((element[a], element[b])))
                    if (
                        edges.setdefault(key, mid) != mid
                        or midpoints.setdefault(mid, key) != key
                    ):
                        raise ValueError(
                            "Nonconforming quadratic edge-node identities."
                        )
        for index, element in enumerate(self.elements):
            p = points[np.array(element) - 1]
            integration_weights(
                self.element_type, tuple(tuple(float(v) for v in row) for row in p)
            )
            for face in face_nodes(self.element_type):
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
        lines += [f"*ELEMENT,TYPE={self.element_type},ELSET=ALLE"]
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


def study_document(study):
    return {"format": "vinkulum-static-study", "schema": 2, "study": asdict(study)}


def load_study(path):
    data = read_json(path, max_bytes=4 * 1024 * 1024)
    if (
        not isinstance(data, dict)
        or set(data) != {"format", "schema", "study"}
        or data["format"] != "vinkulum-static-study"
        or type(data["schema"]) is not int
        or data["schema"] not in (1, 2)
        or not isinstance(data["study"], dict)
    ):
        raise ValueError("Unsupported or incomplete static-study document.")
    if data["schema"] == 1 and data["study"].get("element_type", "C3D8") != "C3D8":
        raise ValueError("Tetrahedral elements require static-study schema 2.")
    if data["schema"] == 2 and "element_type" not in data["study"]:
        raise ValueError("Static-study schema 2 requires an explicit element type.")
    try:
        return StaticStudy(**data["study"])
    except TypeError as exc:
        raise ValueError("Invalid or unknown static-study fields.") from exc


def parse_dat(text, study):
    """Accept exactly one complete final state; preserve integration-point stress."""
    tables = {}
    points_per_element = POINT_COUNTS[study.element_type]
    specifications = (
        ("displacements", "displacements ", "ALLN", len(study.nodes), 4),
        ("external_forces", "forces ", "ALLN", len(study.nodes), 4),
        ("stress", "stresses ", "ALLE", len(study.elements) * points_per_element, 8),
        (
            "energy_density",
            "internal energy density ",
            "ALLE",
            len(study.elements) * points_per_element,
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
                [
                    (e, p)
                    for e in range(1, len(study.elements) + 1)
                    for p in range(1, points_per_element + 1)
                ]
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
    weights = [
        integration_weights(
            study.element_type,
            tuple(
                tuple(float(v) for v in row) for row in points[np.array(element) - 1]
            ),
        )
        for element in study.elements
    ]
    internal = float(
        np.sum(
            np.array(weights) * tables["energy_density"].reshape(-1, points_per_element)
        )
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


def identify_engine(executable):
    found = shutil.which(str(executable))
    if not found:
        raise FileNotFoundError(
            "CalculiX executable unavailable; install ccx separately."
        )
    engine = Path(found).resolve()
    engine_hash = hashlib.sha256(engine.read_bytes()).hexdigest()
    return engine, engine_hash


def parse_version(text):
    version = re.search(r"Version\s+([0-9.]+)", text)
    if not version:
        raise ValueError(
            "The executable did not identify itself as a supported CalculiX CLI."
        )
    return version[1]


@dataclass(frozen=True)
class PreparedRun:
    study: StaticStudy
    root: Path
    engine: Path
    engine_hash: str
    version: str
    deck: str


def prepare_run(study, directory, engine, engine_hash, version):
    """Capture inputs once, shared by synchronous CLI and direct Qt supervision."""
    if not isinstance(study, StaticStudy):
        raise TypeError("A validated StaticStudy is required.")
    root = Path(directory).resolve()
    root.mkdir(parents=True, exist_ok=False)
    deck = study.input_deck()
    (root / "study.inp").write_text(deck, encoding="ascii")
    write_json(
        root / "study.json",
        study_document(study),
    )
    return PreparedRun(study, root, engine, engine_hash, version, deck)


def _read_bounded(path, limit=MAX_OUTPUT_BYTES):
    with Path(path).open("rb") as stream:
        raw = stream.read(limit + 1)
    if len(raw) > limit:
        raise ValueError(f"CalculiX file exceeded its size budget: {Path(path).name}")
    return raw


def _result_contract(study=None):
    contract = {
        "format": "vinkulum-calculix-static",
        "schema": 1,
        "engine": "CalculiX",
        "adapter_version": "0.1.0",
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
    }

    if study is not None:
        contract.update(
            schema=2,
            adapter_version="0.2.0",
            element_type=study.element_type,
            integration_points_per_element=POINT_COUNTS[study.element_type],
            stress_location=f"element integration points 1..{POINT_COUNTS[study.element_type]}, element-major order",
            limitations="Linear isotropic elasticity, C3D4/C3D10 or affine C3D8, zero supports, nodal loads. Local element checks do not certify global mesh injectivity or discretisation error. ASCII output is rounded.",
        )
    return contract


def _successful_log(raw, returncode=0):
    text = raw.decode("utf-8", errors="replace")
    return returncode == 0 and "Job finished" in text and "*ERROR" not in text.upper()


def finish_run(run, returncode):
    """Validate raw files and publish a result only after successful completion."""
    study, root, engine = run.study, run.root, run.engine
    if _read_bounded(root / "study.inp", 4 * 1024 * 1024) != run.deck.encode("ascii"):
        raise ValueError("CalculiX input deck changed during the calculation.")
    log_path, dat_path = root / "solver.log", root / "study.dat"
    if not _successful_log(_read_bounded(log_path), returncode):
        raise RuntimeError(f"CalculiX did not complete; inspect {log_path}.")
    # Hash and parse the same captured bytes, even if a file changes afterwards.
    raw = _read_bounded(dat_path)
    values = parse_dat(raw.decode("ascii"), study)
    if hashlib.sha256(engine.read_bytes()).hexdigest() != run.engine_hash:
        raise RuntimeError("CalculiX executable changed during the calculation.")
    report = {
        **_result_contract(study),
        "engine_version": run.version,
        "engine_sha256": run.engine_hash,
        "input_sha256": hashlib.sha256(run.deck.encode("ascii")).hexdigest(),
        "dat_sha256": hashlib.sha256(raw).hexdigest(),
        **values,
    }
    write_json(root / "result.json", report)
    return report


def load_static_result(directory):
    """Recheck an existing run without executing a solver or modifying its files.

    Hashes and raw-data comparisons establish internal consistency, not the
    authorship of an archive or the historical identity of its executable.
    """
    root = Path(directory).resolve()
    if root.is_file() and root.name == "result.json":
        root = root.parent
    report = read_json(root / "result.json", MAX_OUTPUT_BYTES)
    study = load_study(root / "study.json")
    if (
        not isinstance(report, dict)
        or type(report.get("schema")) is not int
        or report["schema"] not in (1, 2)
    ):
        raise ValueError("Unsupported CalculiX result schema.")
    if report["schema"] == 1 and study.element_type != "C3D8":
        raise ValueError("Legacy CalculiX results require C3D8.")
    contract = _result_contract(study if report["schema"] == 2 else None)
    value_keys = {
        "displacements",
        "external_forces",
        "reactions",
        "stress",
        "energy_density",
        "strain_energy_J",
    }
    hash_keys = {"engine_sha256", "input_sha256", "dat_sha256"}
    if not isinstance(report, dict) or set(report) != set(
        contract
    ) | value_keys | hash_keys | {"engine_version"}:
        raise ValueError("Unsupported or incomplete CalculiX result document.")
    for key, expected in contract.items():
        if type(report[key]) is not type(expected) or report[key] != expected:
            raise ValueError(f"Unsupported CalculiX result contract: {key}.")
    if (
        not isinstance(report["engine_version"], str)
        or re.fullmatch(r"\d+(?:\.\d+)+", report["engine_version"]) is None
    ):
        raise ValueError("Invalid captured CalculiX version.")
    for key in hash_keys:
        if (
            not isinstance(report[key], str)
            or re.fullmatch(r"[0-9a-f]{64}", report[key]) is None
        ):
            raise ValueError(f"Invalid archive fingerprint: {key}.")
    deck = _read_bounded(root / "study.inp", 4 * 1024 * 1024)
    if hashlib.sha256(deck).hexdigest() != report[
        "input_sha256"
    ] or deck != study.input_deck().encode("ascii"):
        raise ValueError("Archived input deck and captured study do not match.")
    if not _successful_log(_read_bounded(root / "solver.log")):
        raise ValueError("Archived log does not describe a completed calculation.")
    raw = _read_bounded(root / "study.dat")
    if hashlib.sha256(raw).hexdigest() != report["dat_sha256"]:
        raise ValueError("Archived raw CalculiX output fingerprint does not match.")
    values = parse_dat(raw.decode("ascii"), study)
    for key, expected in values.items():
        actual = report[key]
        if key == "strain_energy_J":
            if not finite(actual) or not math.isclose(
                actual, expected, rel_tol=1e-12, abs_tol=0
            ):
                raise ValueError("Archived total energy disagrees with raw output.")
        else:
            if not isinstance(actual, list) or any(
                not isinstance(row, list) or not all(finite(v) for v in row)
                for row in actual
            ):
                raise ValueError(f"Invalid numeric archive channel: {key}.")
            if not np.array_equal(actual, expected):
                raise ValueError(f"Archived {key} values disagree with raw output.")
    return study, report, root


def run_static(study, directory, *, executable="ccx", timeout=60):
    """Run in a new directory. Keep the input, raw output and diagnostics on failure."""
    if not isinstance(study, StaticStudy):
        raise TypeError("A validated StaticStudy is required.")
    if not finite(timeout) or not 0 < timeout <= 300:
        raise ValueError("Timeout must be between 0 and 300 seconds.")
    engine, engine_hash = identify_engine(executable)
    version_run = subprocess.run(
        [str(engine), "-v"], capture_output=True, text=True, timeout=10, check=False
    )
    run = prepare_run(
        study, directory, engine, engine_hash, parse_version(version_run.stdout)
    )
    with (run.root / "solver.log").open("wb") as log:
        result = subprocess.run(
            [str(engine), "-i", "study"],
            cwd=run.root,
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
    return finish_run(run, result.returncode)


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
