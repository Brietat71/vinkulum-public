"""Independent shape/volume references, curved traction patch and real CalculiX."""

import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

import numpy as np
from vinkulum_studio.calculix import (
    StaticStudy,
    load_static_result,
    load_study,
    run_static,
    study_document,
)
from vinkulum_studio.finite_elements import (
    TET_EDGES,
    TET_FACES,
    certify_tetra10,
    integration_weights,
    tetra_quadrature,
    tetra_shape,
)


def reference_tetra():
    p = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    return p + [tuple((x + y) / 2 for x, y in zip(p[a], p[b])) for a, b in TET_EDGES]


# An independent polynomial boundary integral. This does not call the element
# module's shape gradients, Jacobian evaluator or volume quadrature.
def add(*polynomials):
    out = {}
    for poly in polynomials:
        for key, value in poly.items():
            out[key] = out.get(key, 0) + value
    return out


def times(a, b):
    out = {}
    for (i, j), v in a.items():
        for (k, l), w in b.items():
            key = (i + k, j + l)
            out[key] = out.get(key, 0) + v * w
    return out


def constant(value):
    return {(0, 0): value}


def derivative(poly, axis):
    out = {}
    for key, value in poly.items():
        if key[axis]:
            power = list(key)
            power[axis] -= 1
            out[tuple(power)] = value * key[axis]
    return out


def triangle_forces(points, quadratic, stress, bending=False):
    bary = ({(0, 0): 1, (1, 0): -1, (0, 1): -1}, {(1, 0): 1}, {(0, 1): 1})
    if quadratic:
        basis = [times(b, add(times(constant(2), b), constant(-1))) for b in bary]
        basis += [
            times(constant(4), times(bary[a], bary[b]))
            for a, b in ((0, 1), (1, 2), (2, 0))
        ]
    else:
        basis = list(bary)
    coords = [
        add(*(times(constant(p[k]), n) for p, n in zip(points, basis)))
        for k in range(3)
    ]
    dr, ds = [[derivative(p, k) for p in coords] for k in (0, 1)]
    normal = [
        add(times(dr[a], ds[b]), times(constant(-1), times(dr[b], ds[a])))
        for a, b in ((1, 2), (2, 0), (0, 1))
    ]
    traction = [
        add(*(times(constant(stress[i, j]), normal[j]) for j in range(3)))
        for i in range(3)
    ]
    if bending:
        sigma_xx = times(constant(-1e6), add(coords[1], constant(-0.1)))
        traction = [times(sigma_xx, normal[0]), constant(0), constant(0)]
    return np.array(
        [
            [
                sum(
                    v
                    * math.factorial(i)
                    * math.factorial(j)
                    / math.factorial(i + j + 2)
                    for (i, j), v in times(n, t).items()
                )
                for t in traction
            ]
            for n in basis
        ]
    )


def tetra_specimen(kind="C3D10", *, curved=False, cells=(1, 1, 1), bending=False):
    nx, ny, nz = cells
    nodes = [
        (i / nx, 0.2 * j / ny, 0.3 * k / nz)
        for k in range(nz + 1)
        for j in range(ny + 1)
        for i in range(nx + 1)
    ]

    def node(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)

    elements = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                cube = [
                    node(i, j, k),
                    node(i + 1, j, k),
                    node(i + 1, j + 1, k),
                    node(i, j + 1, k),
                    node(i, j, k + 1),
                    node(i + 1, j, k + 1),
                    node(i + 1, j + 1, k + 1),
                    node(i, j + 1, k + 1),
                ]
                for row in (
                    (0, 1, 2, 6),
                    (0, 2, 3, 6),
                    (0, 3, 7, 6),
                    (0, 7, 4, 6),
                    (0, 4, 5, 6),
                    (0, 5, 1, 6),
                ):
                    elements.append(tuple(cube[a] for a in row))
    if kind == "C3D10":
        mids, enriched = {}, []
        for element in elements:
            row = list(element)
            for a, b in TET_EDGES:
                edge = tuple(sorted((element[a], element[b])))
                if edge not in mids:
                    mids[edge] = len(nodes)
                    nodes.append(
                        tuple(
                            (x + y) / 2 for x, y in zip(nodes[edge[0]], nodes[edge[1]])
                        )
                    )
                row.append(mids[edge])
            enriched.append(tuple(row))
        elements = enriched
    if curved:
        nodes = [(x, y, z * (1 + 0.25 * x)) for x, y, z in nodes]
    fixed = [(i + 1, 1) for i, p in enumerate(nodes) if p[0] == 0]
    fixed += [(1, 2), (1, 3), (node(0, ny, 0) + 1, 3)]
    faces = {}
    for element in elements:
        for face in TET_FACES:
            corners = tuple(element[i] for i in face)
            row = list(corners)
            if kind == "C3D10":
                lookup = {
                    frozenset((element[a], element[b])): element[4 + i]
                    for i, (a, b) in enumerate(TET_EDGES)
                }
                row += [
                    lookup[frozenset((corners[a], corners[b]))]
                    for a, b in ((0, 1), (1, 2), (2, 0))
                ]
            faces.setdefault(tuple(sorted(corners)), []).append(row)
    forces = np.zeros((len(nodes), 3))
    for owners in faces.values():
        if len(owners) != 1:
            continue
        row = owners[0]
        stress = np.diag((1e5, 0.0, 0.0))
        forces[row] += triangle_forces(
            np.array(nodes)[row], kind == "C3D10", stress, bending=bending
        )
    return StaticStudy(
        nodes,
        [tuple(i + 1 for i in e) for e in elements],
        fixed,
        [
            (i + 1, *map(float, f))
            for i, f in enumerate(forces)
            if np.max(abs(f)) > 1e-10
        ],
        210e9,
        0.3,
        element_type=kind,
    )


class TetraGeometry(unittest.TestCase):
    def test_saved_public_bending_opens_without_gui_or_solver(self):
        archive = (
            Path(__file__).resolve().parents[3]
            / "docs/bancs/studio-tetrahedra-060/calculation"
        )
        code = """
import sys
from vinkulum_studio.calculix import load_static_result
study, report, _ = load_static_result(sys.argv[1])
assert study.element_type == "C3D10" and len(study.elements) == 6
assert abs(report["strain_energy_J"] / (1e12*(.3*.2**3/12)/(2*210e9))-1) < 6e-7
assert not any(name in sys.modules for name in ("PySide6", "vtkmodules", "OCP", "build123d", "pinocchio"))
print("Saved quadratic bending archive checked without GUI or solver")
"""
        checked = subprocess.run(
            [sys.executable, "-c", code, str(archive)],
            env={**os.environ, "PATH": ""},
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        self.assertEqual(checked.returncode, 0, checked.stderr)

    def test_positive_curved_element_needs_subdivision_and_budgets_refuse(self):
        nodes = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.34375, 0.125, 0.0),
            (0.46875, 0.59375, -0.09375),
            (-0.0625, 0.4375, 0.03125),
            (0.15625, -0.125, 0.46875),
            (0.375, -0.15625, 0.53125),
            (0.09375, 0.53125, 0.625),
        ]
        self.assertEqual(certify_tetra10(nodes, max_depth=0).status, "unresolved")
        self.assertEqual(certify_tetra10(nodes, max_cells=1).status, "unresolved")
        certificate = certify_tetra10(nodes)
        self.assertEqual(certificate.status, "positive")
        self.assertGreater(certificate.cells_examined, 1)
        self.assertEqual(certificate.lower_bound, Fraction(13, 768))

    def test_basis_interpolation_partition_and_gradient(self):
        nodes = reference_tetra()
        np.testing.assert_allclose(
            [tetra_shape(p)[0] for p in nodes], np.eye(10), atol=1e-15
        )
        p = np.array((0.23, 0.19, 0.11))
        values, gradients = tetra_shape(p)
        self.assertAlmostEqual(sum(values), 1)
        np.testing.assert_allclose(gradients.sum(axis=0), 0, atol=1e-15)
        for k in range(3):
            delta = np.eye(3)[k] * 1e-5
            np.testing.assert_allclose(
                (tetra_shape(p + delta)[0] - tetra_shape(p - delta)[0]) / 2e-5,
                gradients[:, k],
                atol=2e-11,
            )

    def test_exact_bernstein_bound_and_curved_quadrature(self):
        nodes = reference_tetra()
        self.assertEqual(certify_tetra10(nodes).lower_bound, Fraction(1))
        # x=r+r²/4, y=s, z=t. J=1+r/2, so J>=1 and volume=3/16.
        curved = [(x + x * x / 4, y, z) for x, y, z in nodes]
        certificate = certify_tetra10(curved)
        self.assertEqual(certificate.status, "positive")
        self.assertEqual(certificate.lower_bound, Fraction(1))
        weights = integration_weights("C3D10", tuple(curved))
        self.assertAlmostEqual(sum(weights), 3 / 16, places=14)
        for w, (q, _) in zip(weights, tetra_quadrature("C3D10")):
            self.assertAlmostEqual(w, (1 + q[0] / 2) / 24, places=15)
        shifted = [(x + 2**20, y - 2**20, z + 2**20) for x, y, z in curved]
        self.assertEqual(certify_tetra10(shifted).lower_bound, certificate.lower_bound)

    def test_inversion_hidden_from_all_gauss_points_is_rejected(self):
        p = reference_tetra()
        p[4] = (0.24, 0.0, 0.0)
        determinants = [
            np.linalg.det(np.array(p).T @ tetra_shape(q)[1])
            for q, _ in tetra_quadrature("C3D10")
        ]
        self.assertGreater(min(determinants), 0)
        certificate = certify_tetra10(p)
        self.assertEqual(certificate.status, "nonpositive")
        self.assertEqual(certificate.witness, (0.0, 0.0, 0.0))
        with self.assertRaisesRegex(ValueError, "nonpositive"):
            integration_weights("C3D10", tuple(p))

    def test_schema_types_quadratic_connectivity_and_legacy(self):
        study = tetra_specimen()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "study.json"
            data = study_document(study)
            path.write_text(json.dumps(data))
            self.assertEqual(load_study(path), study)
            path.write_text(json.dumps({**data, "schema": 1}))
            with self.assertRaisesRegex(ValueError, "schema 2"):
                load_study(path)
        with self.assertRaises(ValueError):
            replace(study, element_type="C3D8R")
        with self.assertRaises(ValueError):
            replace(study, element_type="C3D4")
        element = list(study.elements[1])
        shared = next(
            (a, b) for a, b in enumerate(element[4:], 4) if b in study.elements[0][4:]
        )
        element[shared[0]] = len(study.nodes) + 1
        bad_nodes = study.nodes + (study.nodes[shared[1] - 1],)
        with self.assertRaisesRegex(ValueError, "Nonconforming"):
            replace(
                study,
                nodes=bad_nodes,
                elements=(study.elements[0], tuple(element), *study.elements[2:]),
            )


@unittest.skipUnless(
    shutil.which("ccx"), "Install CalculiX for independent physical references"
)
class TetraCalculix(unittest.TestCase):
    def test_pure_bending_quadratic_exactness_and_linear_refinement(self):
        errors = []
        exact_energy = 0.3 * 0.2**3 / 12 * 1e12 / (2 * 210e9)
        for kind, cells in (
            ("C3D4", (5, 1, 1)),
            ("C3D4", (10, 2, 2)),
            ("C3D4", (20, 4, 4)),
            ("C3D10", (1, 1, 1)),
        ):
            with (
                self.subTest(kind=kind, cells=cells),
                tempfile.TemporaryDirectory() as directory,
            ):
                study = tetra_specimen(kind, cells=cells, bending=True)
                report = run_static(study, Path(directory) / "run")
                relative = abs(report["strain_energy_J"] / exact_energy - 1)
                if kind == "C3D4":
                    self.assertLess(report["strain_energy_J"], exact_energy)
                    errors.append(relative)
                else:
                    self.assertLess(relative, 6e-7)
                    x, y, z = np.array(study.nodes).T
                    kappa = 1e6 / study.young_pa
                    expected = kappa * np.column_stack(
                        (
                            -x * (y - 0.1),
                            0.5 * (x * x + 0.3 * ((y - 0.1) ** 2 - z * z - 0.1**2)),
                            0.3 * (y - 0.1) * z,
                        )
                    )
                    np.testing.assert_allclose(
                        report["displacements"], expected, rtol=6e-7, atol=1e-14
                    )
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        # This is a convergence check, not a promise of 20% accuracy. Even
        # 1920 linear tetrahedra retain about 20.5% energy error in this case.
        self.assertLess(errors[2], errors[0] / 2)

    def test_curved_bending_preserves_physical_quadrature_weights(self):
        study = tetra_specimen(curved=True, bending=True)
        with tempfile.TemporaryDirectory() as directory:
            report = run_static(study, Path(directory) / "run")
            weights = np.array(
                [
                    integration_weights("C3D10", tuple(study.nodes[i - 1] for i in e))
                    for e in study.elements
                ]
            )
            density = np.array(report["energy_density"]).reshape(-1, 4)
            actual = np.sum(weights * density)
            swapped = np.sum(weights[:, [1, 0, 2, 3]] * density)
            self.assertAlmostEqual(actual / report["strain_energy_J"], 1, places=13)
            self.assertGreater(abs(swapped / actual - 1), 1e-4)

    def test_linear_and_quadratic_straight_and_curved_traction_patches(self):
        for kind, curved in (("C3D4", False), ("C3D10", False), ("C3D10", True)):
            with (
                self.subTest(kind=kind, curved=curved),
                tempfile.TemporaryDirectory() as directory,
            ):
                study = tetra_specimen(kind, curved=curved)
                root = Path(directory) / "run"
                report = run_static(study, root)
                strain = 1e5 / study.young_pa
                expected = np.array(study.nodes) * np.array(
                    (strain, -0.3 * strain, -0.3 * strain)
                )
                np.testing.assert_allclose(
                    report["displacements"], expected, rtol=6e-7, atol=1e-14
                )
                np.testing.assert_allclose(
                    report["stress"],
                    np.tile((1e5, 0, 0, 0, 0, 0), (len(report["stress"]), 1)),
                    rtol=6e-7,
                    atol=2e-7,
                )
                volume = 0.06 * (1.125 if curved else 1)
                self.assertAlmostEqual(
                    report["strain_energy_J"] / (0.5 * 1e10 / study.young_pa * volume),
                    1,
                    delta=6e-7,
                )
                reopened, checked, _ = load_static_result(root)
                self.assertEqual(reopened, study)
                self.assertEqual(checked, report)
                self.assertEqual(
                    report["integration_points_per_element"], 1 if kind == "C3D4" else 4
                )


if __name__ == "__main__":
    unittest.main()
