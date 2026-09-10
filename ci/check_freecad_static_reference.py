"""NumPy-only affine-tension reference for the two retained FreeCAD FEM cases.

SPDX-License-Identifier: Apache-2.0
Checks captured preview numbers, not CAD face correspondence or general FE error.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
from zipfile import ZipFile

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check_case(request, data, case):
    require(case in (0, 1), "Unknown affine reference case")
    require(
        request["young_pa"] == 210e9 and request["poisson"] == 0,
        "Reference material changed",
    )
    require(request["density_kg_m3"] == 7800, "Reference density changed")
    conditions = request["conditions"]
    require(len(conditions) == 2, "Reference boundary conditions changed")
    support, pressure = conditions
    require(
        support["kind"] == "support"
        and support["faces"] == [0]
        and support["axes"] == [1, 2, 3],
        "Reference support changed",
    )
    require(
        pressure["kind"] == "pressure"
        and pressure["faces"] == [1]
        and pressure["values"] == [-2e6],
        "Reference pressure changed",
    )
    frame, shift = np.eye(3), np.zeros(3)
    if case:
        axis = np.array((1.0, 2.0, 3.0)) / math.sqrt(14)
        x, y, z = axis
        cross = np.array(((0, -z, y), (z, 0, -x), (-y, x, 0)))
        angle = math.radians(37)
        frame = (
            math.cos(angle) * np.eye(3)
            + (1 - math.cos(angle)) * np.outer(axis, axis)
            + math.sin(angle) * cross
        )
        shift = np.array((0.25, -0.1, 0.5))
    dimensions = np.array((0.12, 0.02, 0.015))
    volume = float(np.prod(dimensions))
    mass = 7800 * volume
    inertia = (
        frame @ np.diag(mass * (sum(dimensions**2) - dimensions**2) / 12) @ frame.T
    )
    expected_properties = {
        "volume_m3": volume,
        "mass_kg": mass,
        "centre_m": frame @ (dimensions / 2) + shift,
        "inertia_kg_m2": inertia.ravel(),
    }
    length = float(np.linalg.norm(dimensions))
    scales = (volume, mass, length, mass * length**2)
    for (name, expected), scale in zip(expected_properties.items(), scales):
        require(
            np.allclose(
                request["properties_si"][name], expected, atol=1e-8 * scale, rtol=0
            ),
            "Captured property disagrees with analytic box: " + name,
        )
    nodes, displacements, reactions = (
        np.array(data[key], dtype=float)
        for key in ("nodes_m", "displacements_m", "reactions_N")
    )
    require(
        nodes.ndim == 2 and nodes.shape[1] == 3 and len(nodes) >= 4,
        "Invalid node coordinates",
    )
    require(
        displacements.shape == reactions.shape == nodes.shape,
        "Invalid result dimensions",
    )
    require(
        all(np.isfinite(a).all() for a in (nodes, displacements, reactions)),
        "Nonfinite reference data",
    )
    local = (nodes - shift) @ frame
    require(
        np.all(local >= -1e-9) and np.all(local <= dimensions + 1e-9),
        "Nodes leave the reference box",
    )
    expected_displacement = (2e6 / 210e9) * local[:, 0, None] * frame[:, 0]
    displacement_error = float(
        np.max(np.linalg.norm(displacements - expected_displacement, axis=1))
    )
    expected_energy = 2e6**2 * volume / (2 * 210e9)
    energy_error = abs(data["strain_energy_J"] / expected_energy - 1)
    force = 2e6 * dimensions[1] * dimensions[2]
    reaction_error = float(np.linalg.norm(reactions.sum(axis=0) + force * frame[:, 0]))
    require(displacement_error < 1e-10, "Displacement differs from affine tension")
    require(energy_error < 5e-6, "Energy differs from affine tension")
    require(reaction_error < 5e-6 * force, "Net reaction does not balance the traction")
    return {
        "case": case,
        "nodes": len(nodes),
        "max_displacement_error_m": displacement_error,
        "relative_energy_error": energy_error,
        "net_reaction_error_N": reaction_error,
    }


def check_archive(path):
    results = []
    with ZipFile(path) as archive:
        for case in (0, 1):
            prefix = f"case-{case}/"
            raw = archive.read(prefix + "request.json")
            request = json.loads(raw)
            data = json.loads(archive.read(prefix + "preview.json"))
            require(
                request["format"] == "vinkulum-freecad-static-1"
                and data["format"] == "vinkulum-freecad-static-preview-1",
                "Unsupported capture format",
            )
            require(
                data["run_id"] == request["run_id"]
                and data["request_sha256"] == hashlib.sha256(raw).hexdigest(),
                "Preview identifies a different capture",
            )
            files = [("part.step", request["step_sha256"])] + [
                (f"face-{i:03}.step", face["step_sha256"])
                for i, face in enumerate(request["faces"])
            ]
            for name, digest in files:
                require(
                    hashlib.sha256(archive.read(prefix + name)).hexdigest() == digest,
                    "Changed STEP bytes: " + name,
                )
            results.append(check_case(request, data, case))
    return {
        "status": "passed",
        "reference": "Affine tension with independent Rodrigues frame transport, box properties, energy and reaction balance",
        "cases": results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    print(json.dumps(check_archive(parser.parse_args().archive), indent=2))
