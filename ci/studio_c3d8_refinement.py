"""Reproduce a C3D8 cantilever study and compare it with explicit beam estimates."""

import argparse
import json
import platform
import re
import tempfile
import time
import zipfile
from pathlib import Path

import numpy as np
from vinkulum_studio.calculix import StaticStudy, load_static_result, run_static

PROBLEM = {
    "length_m": 0.12,
    "width_m": 0.012,
    "height_m": 0.006,
    "young_pa": 210e9,
    "poisson": 0.3,
    "force_z_N": -1.0,
    "shear_factor": 5 / 6,
}
MESHES = ((20, 2, 1), (40, 4, 2), (80, 8, 4))
FILES = ("study.json", "study.inp", "study.dat", "solver.log", "result.json")


def beam_reference():
    p = PROBLEM
    length, width, height = (p[k] for k in ("length_m", "width_m", "height_m"))
    area, second_moment = width * height, width * height**3 / 12
    shear_modulus = p["young_pa"] / (2 * (1 + p["poisson"]))
    bending = p["force_z_N"] * length**3 / (3 * p["young_pa"] * second_moment)
    shear = p["force_z_N"] * length / (p["shear_factor"] * shear_modulus * area)
    return {
        "area_m2": area,
        "second_moment_y_m4": second_moment,
        "euler_bernoulli_uz_m": bending,
        "timoshenko_uz_m": bending + shear,
        "euler_bernoulli_energy_J": p["force_z_N"] * bending / 2,
        "timoshenko_energy_J": p["force_z_N"] * (bending + shear) / 2,
    }


def cantilever(cells):
    nx, ny, nz = cells
    length, width, height = (PROBLEM[k] for k in ("length_m", "width_m", "height_m"))

    def node(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    nodes = tuple(
        (length * i / nx, width * (j / ny - 0.5), height * (k / nz - 0.5))
        for k in range(nz + 1)
        for j in range(ny + 1)
        for i in range(nx + 1)
    )
    corners = (
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    )
    elements = tuple(
        tuple(node(i + a, j + b, k + c) for a, b, c in corners)
        for k in range(nz)
        for j in range(ny)
        for i in range(nx)
    )
    fixed = tuple(
        (node(0, j, k), axis)
        for k in range(nz + 1)
        for j in range(ny + 1)
        for axis in (1, 2, 3)
    )
    # Integral of a bilinear face shape function: one quarter of its face area.
    loads = {}
    for k in range(nz):
        for j in range(ny):
            for dj, dk in ((0, 0), (1, 0), (1, 1), (0, 1)):
                identifier = node(nx, j + dj, k + dk)
                loads[identifier] = loads.get(identifier, 0.0) + PROBLEM[
                    "force_z_N"
                ] / (4 * ny * nz)
    return StaticStudy(
        nodes,
        elements,
        fixed,
        tuple((i, 0.0, 0.0, f) for i, f in sorted(loads.items())),
        PROBLEM["young_pa"],
        PROBLEM["poisson"],
    )


def measurement(study, result, log, cells):
    actual = study.solver_study
    nodes, displacement = np.asarray(actual.nodes), np.asarray(result["displacements"])
    applied = np.zeros_like(nodes)
    for node, *force in actual.forces:
        applied[node - 1] = force
    # Area average over the loaded face: these weights also produce the traction.
    average = float(np.dot(applied[:, 2] / PROBLEM["force_z_N"], displacement[:, 2]))
    external = applied + np.asarray(result["reactions"])
    timer = re.findall(r"Total CalculiX Time:\s+([0-9.Ee+-]+)", log)
    if len(timer) != 1:
        raise ValueError("A unique solver-reported elapsed time is required.")
    beam = beam_reference()
    return {
        "cells": list(cells),
        "nodes": len(nodes),
        "elements": len(study.elements),
        "total_dofs": 3 * len(nodes),
        "free_dofs": 3 * len(nodes) - len(study.fixed_dofs),
        "loaded_end_uz_m": average,
        "strain_energy_J": result["strain_energy_J"],
        "half_external_work_J": float(0.5 * np.sum(applied * displacement)),
        "force_residual_N": external.sum(axis=0).tolist(),
        "moment_residual_Nm": np.cross(nodes, external).sum(axis=0).tolist(),
        "gap_to_euler_bernoulli_percent": 100
        * (average / beam["euler_bernoulli_uz_m"] - 1),
        "gap_to_timoshenko_percent": 100 * (average / beam["timoshenko_uz_m"] - 1),
        "solver_elapsed_s": float(timer[0]),
        "engine_version": result["engine_version"],
        "requested_threads": result["requested_threads"],
        "scientific_status": result["scientific_status"],
    }


def table(report):
    lines = [
        "| Cells X×Y×Z | Free DOFs | Mean end uz (µm) | Strain energy (J) | Gap to EB (%) | Gap to Timoshenko (%) | Force residual (N) | Moment residual (N m) | Solver time (s) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in report["cases"]:
        lines.append(
            f"| {'×'.join(map(str, row['cells']))} | {row['free_dofs']} | {row['loaded_end_uz_m'] * 1e6:.8g} | "
            f"{row['strain_energy_J']:.8g} | {row['gap_to_euler_bernoulli_percent']:.5f} | "
            f"{row['gap_to_timoshenko_percent']:.5f} | {np.linalg.norm(row['force_residual_N']):.3g} | "
            f"{np.linalg.norm(row['moment_residual_Nm']):.3g} | {row['solver_elapsed_s']:.6g} |"
        )
    return "\n".join(lines) + "\n"


def run(output, executable):
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "problem": PROBLEM,
        "beam_reference": beam_reference(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cases": [],
    }
    for cells in MESHES:
        case = "mesh-" + "-".join(map(str, cells))
        study = cantilever(cells)
        # Keep the original calculation directory, including diagnostics on failure.
        root = output / case
        start = time.perf_counter()
        result = run_static(study, root, executable=executable, timeout=120)
        elapsed = time.perf_counter() - start
        row = measurement(study, result, (root / "solver.log").read_text(), cells)
        row.update(archive=case + ".zip", adapter_elapsed_s=elapsed)
        with zipfile.ZipFile(
            output / row["archive"],
            "w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
        ) as archive:
            for filename in FILES:
                archive.write(root / filename, filename)
        report["cases"].append(row)
        print(json.dumps(row), flush=True)
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    (output / "table.md").write_text(table(report))
    return report


def verify(directory):
    report = json.loads((directory / "summary.json").read_text())
    if report["problem"] != PROBLEM or report["beam_reference"] != beam_reference():
        raise ValueError(
            "The retained study or beam reference differs from this recipe."
        )
    if [row["cells"] for row in report["cases"]] != list(map(list, MESHES)):
        raise ValueError(
            "The three specified systematic mesh refinements are required."
        )
    for row, cells in zip(report["cases"], MESHES):
        if not np.isfinite(row["adapter_elapsed_s"]) or row["adapter_elapsed_s"] <= 0:
            raise ValueError("A finite positive adapter elapsed time is required.")
        archive_path = directory / ("mesh-" + "-".join(map(str, cells)) + ".zip")
        if row["archive"] != archive_path.name:
            raise ValueError("Unexpected archive identity.")
        with (
            tempfile.TemporaryDirectory() as temporary,
            zipfile.ZipFile(archive_path) as archive,
        ):
            if sorted(archive.namelist()) != sorted(FILES):
                raise ValueError(
                    "Archive contents must match the five original calculation files."
                )
            root = Path(temporary)
            for name in FILES:
                if archive.getinfo(name).file_size > 32 * 1024 * 1024:
                    raise ValueError("Archived calculation exceeds the adapter budget.")
                (root / name).write_bytes(archive.read(name))
            study, result, _ = load_static_result(root)
            if study != cantilever(cells):
                raise ValueError(
                    "The retained mesh, material, supports or traction changed."
                )
            measured = measurement(
                study, result, (root / "solver.log").read_text(), cells
            )
            for key, value in measured.items():
                if isinstance(value, str) or key in (
                    "cells",
                    "nodes",
                    "elements",
                    "total_dofs",
                    "free_dofs",
                    "requested_threads",
                ):
                    if value != row[key]:
                        raise ValueError(f"{key}: inconsistent summary")
                else:
                    # Re-aggregation tolerance only; not a discretisation-error bound.
                    scale = {"force_residual_N": 1.0, "moment_residual_Nm": 0.12}.get(
                        key, abs(float(value)) if np.ndim(value) == 0 else 1.0
                    )
                    np.testing.assert_allclose(
                        value, row[key], rtol=1e-12, atol=1e-12 * scale, err_msg=key
                    )
    if (directory / "table.md").read_text() != table(report):
        raise ValueError("The displayed table differs from the retained measurements.")
    return report


def plot(directory):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report = verify(directory)
    counts = [row["elements"] for row in report["cases"]]
    displacement = [-1e6 * row["loaded_end_uz_m"] for row in report["cases"]]
    fig, axis = plt.subplots(figsize=(8.5, 5), layout="constrained")
    axis.semilogx(counts, displacement, "o-", label="C3D8 · area-mean end displacement")
    for label, key, style in (
        ("Euler–Bernoulli", "euler_bernoulli_uz_m", "--"),
        ("Timoshenko (κ = 5/6)", "timoshenko_uz_m", ":"),
    ):
        value = -1e6 * report["beam_reference"][key]
        axis.axhline(
            value, linestyle=style, color="0.4", label=f"{label}: {value:.5f} µm"
        )
    for index, (x, y, row) in enumerate(zip(counts, displacement, report["cases"])):
        axis.annotate(
            f"{row['gap_to_euler_bernoulli_percent']:.2f}% vs EB",
            (x, y),
            xytext=(0, -19),
            textcoords="offset points",
            ha=(
                "left"
                if index == 0
                else "right" if index == len(counts) - 1 else "center"
            ),
        )
    axis.set(
        xticks=counts,
        xticklabels=list(map(str, counts)),
        xlabel="C3D8 elements · all three mesh directions refined",
        ylabel="Downward displacement magnitude (µm)",
        title="Clamped 120 × 12 × 6 mm cantilever · uniform end traction, total 1 N",
    )
    axis.set_ylim(7.3, 14.5)
    axis.grid(alpha=0.25)
    axis.legend(loc="upper left", fontsize=9)
    fig.text(
        0.5,
        -0.015,
        "Beam estimates are asymptotic comparisons, not the exact 3D clamped solution.",
        ha="center",
        fontsize=9,
    )
    fig.savefig(directory / "refinement.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("run", "verify", "plot"))
    parser.add_argument("directory", type=Path)
    parser.add_argument("--ccx", default="ccx")
    args = parser.parse_args()
    if args.action == "run":
        run(args.directory, args.ccx)
    elif args.action == "verify":
        print(table(verify(args.directory)), end="")
    else:
        plot(args.directory)
