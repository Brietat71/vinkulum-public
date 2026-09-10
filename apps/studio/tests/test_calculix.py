"""Independent affine patch tests and failure paths for the optional FEM adapter."""

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
from vinkulum_studio.calculix import StaticStudy, load_study, parse_dat, run_static

CCX = shutil.which("ccx")


def specimen(
    length=1.0,
    width=0.1,
    height=0.1,
    young=210e9,
    poisson=0.3,
    force=1000.0,
    cells=(1, 1, 1),
):
    nx, ny, nz = cells
    nodes = tuple(
        (length * i / nx, width * j / ny, height * k / nz)
        for k in range(nz + 1)
        for j in range(ny + 1)
        for i in range(nx + 1)
    )

    def node(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    elements = tuple(
        (
            node(i, j, k),
            node(i + 1, j, k),
            node(i + 1, j + 1, k),
            node(i, j + 1, k),
            node(i, j, k + 1),
            node(i + 1, j, k + 1),
            node(i + 1, j + 1, k + 1),
            node(i, j + 1, k + 1),
        )
        for k in range(nz)
        for j in range(ny)
        for i in range(nx)
    )
    fixed = tuple((node(0, j, k), 1) for k in range(nz + 1) for j in range(ny + 1))
    fixed += ((node(0, 0, 0), 2), (node(0, 0, 0), 3), (node(0, ny, 0), 3))
    # Consistent nodal forces for a uniform traction on each bilinear end face.
    loads = {}
    for k in range(nz):
        for j in range(ny):
            for identifier in (
                node(nx, j, k),
                node(nx, j + 1, k),
                node(nx, j, k + 1),
                node(nx, j + 1, k + 1),
            ):
                loads[identifier] = loads.get(identifier, 0.0) + force / (4 * ny * nz)
    return StaticStudy(
        nodes,
        elements,
        fixed,
        tuple((i, f, 0.0, 0.0) for i, f in sorted(loads.items())),
        young,
        poisson,
    )


class StaticContracts(unittest.TestCase):
    def test_thin_element_cannot_hide_warping_below_long_edge_tolerance(self):
        study = specimen(width=1.0, height=1.0)
        points = [(x, y, z * 2e-12) for x, y, z in study.nodes]
        points[7] = (1.0, 1.0, -3e-12)
        regular = np.array(study.nodes) * (1.0, 1.0, 2e-12)
        self.assertLess(np.max(abs(np.array(points) - regular)), 1e-10)
        # Q1 interpolation of thickness at the far Gauss point is negative.
        far_weight = ((1 + 1 / np.sqrt(3)) / 2) ** 2
        self.assertLess(2e-12 - 5e-12 * far_weight, 0)
        # Clamp all nodes so the geometry rule is exercised independently of
        # rigid-mode admission; the nominal basis has positive determinant.
        with self.assertRaisesRegex(ValueError, "affine"):
            replace(
                study,
                nodes=tuple(points),
                fixed_dofs=tuple((i, a) for i in range(1, 9) for a in (1, 2, 3)),
            )

    def test_study_schema_rejects_other_documents_and_unknown_fields(self):
        from dataclasses import asdict

        data = {
            "format": "vinkulum-static-study",
            "schema": 1,
            "study": asdict(specimen()),
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "study.json"
            path.write_text(json.dumps(data))
            self.assertEqual(load_study(path), specimen())
            for bad in (
                {**data, "format": "vinkulum-studio-project"},
                {**data, "schema": True},
                {**data, "study": {**data["study"], "nonlinear": True}},
            ):
                path.write_text(json.dumps(bad))
                with self.assertRaises(ValueError):
                    load_study(path)

    def test_reject_material_geometry_connectivity_and_free_rigid_modes(self):
        study = specimen()
        invalid = (
            {"young_pa": 0},
            {"poisson": 0.5},
            {"poisson": -1},
            {"young_pa": True},
            {"forces": ((99, 1, 0, 0),)},
            {"forces": ((2, float("nan"), 0, 0),)},
            {"fixed_dofs": ((1, 1), (1, 2), (1, 3))},
            {"fixed_dofs": study.fixed_dofs + study.fixed_dofs[:1]},
            {"elements": study.elements * 2},
            {"elements": ((1, 4, 3, 2, 5, 8, 7, 6),)},
            {"nodes": study.nodes[:6] + ((1, 0.11, 0.1),) + study.nodes[7:]},
        )
        for values in invalid:
            with self.subTest(values=values), self.assertRaises(ValueError):
                replace(study, **values)
        with self.assertRaisesRegex(ValueError, "connected"):
            replace(
                study,
                nodes=study.nodes + tuple((x + 2, y, z) for x, y, z in study.nodes),
                elements=study.elements + (tuple(i + 8 for i in study.elements[0]),),
            )

    def test_zero_forces_do_not_generate_an_empty_cload_card(self):
        self.assertNotIn("*CLOAD", specimen(force=0).input_deck())

    def test_missing_error_and_timeout_preserve_diagnostics(self):
        study = specimen()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaises(FileNotFoundError):
                run_static(study, root / "missing", executable=root / "absent")
            self.assertFalse((root / "missing").exists())
            for name, code, exception in (
                ("error", "print('*ERROR: deliberate test failure')", RuntimeError),
                ("timeout", "import time; time.sleep(10)", subprocess.TimeoutExpired),
            ):
                executable = root / name
                executable.write_text(
                    f'#!{sys.executable}\nimport sys\nif "-v" in sys.argv:\n print("Version 2.21")\nelse:\n {code}\n'
                )
                executable.chmod(0o755)
                output = root / (name + "-output")
                with self.assertRaises(exception):
                    run_static(study, output, executable=executable, timeout=0.1)
                self.assertTrue((output / "study.inp").is_file())
                self.assertTrue((output / "solver.log").is_file())
                self.assertFalse((output / "result.json").exists())
                marker = output / "previous-result"
                marker.write_text("preserve")
                with self.assertRaises(FileExistsError):
                    run_static(study, output, executable=executable)
                self.assertEqual(marker.read_text(), "preserve")


@unittest.skipUnless(CCX, "Install the separate CalculiX ccx executable")
class CalculixReferences(unittest.TestCase):
    def test_world_axis_permutation_and_translation_preserve_physics(self):
        study = specimen(cells=(2, 2, 2))
        rotation = np.array(((0.0, -1.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)))
        translation = np.array((3.0, -2.0, 0.7))
        nodes = np.array(study.nodes) @ rotation.T + translation
        moved = replace(
            study,
            nodes=tuple(tuple(float(v) for v in p) for p in nodes),
            fixed_dofs=tuple(
                (node, {1: 2, 2: 1, 3: 3}[axis]) for node, axis in study.fixed_dofs
            ),
            forces=tuple(
                (row[0], *(float(v) for v in rotation @ row[1:]))
                for row in study.forces
            ),
        )
        with tempfile.TemporaryDirectory() as directory:
            result = run_static(moved, Path(directory) / "run", executable=CCX)
        epsilon = 100000.0 / 210e9
        expected_u = (
            np.array(study.nodes) * (epsilon, -0.3 * epsilon, -0.3 * epsilon)
        ) @ rotation.T
        np.testing.assert_allclose(
            result["displacements"], expected_u, rtol=5e-7, atol=1e-16
        )
        expected_s = np.tile(
            (0.0, 100000.0, 0.0, 0.0, 0.0, 0.0), (8 * len(study.elements), 1)
        )
        np.testing.assert_allclose(result["stress"], expected_s, rtol=5e-7, atol=1e-5)
        np.testing.assert_allclose(
            np.sum(result["reactions"], axis=0),
            (0.0, -1000.0, 0.0),
            rtol=5e-7,
            atol=1e-7,
        )

    def test_uniform_strain_displacement_stress_energy_and_balance(self):
        cases = (
            (1.0, 0.1, 0.1, 210e9, 0.3, 1000.0, (1, 1, 1)),
            (0.01, 0.005, 0.003, 70e9, 0.33, -25.0, (2, 2, 2)),
            (2.0, 1.0, 0.4, 3e6, 0.0, 100.0, (3, 2, 1)),
            (1.0, 0.2, 0.3, 12e9, -0.2, 20.0, (2, 3, 2)),
            (1.0, 0.1, 0.1, 210e9, 0.3, 0.0, (1, 1, 1)),
        )
        for length, width, height, young, poisson, force, cells in cases:
            with (
                self.subTest(cells=cells, force=force, poisson=poisson),
                tempfile.TemporaryDirectory() as directory,
            ):
                study = specimen(length, width, height, young, poisson, force, cells)
                output = Path(directory) / "run"
                result = run_static(study, output, executable=CCX)
                sigma = force / (width * height)
                epsilon = sigma / young
                expected_u = np.array(study.nodes) * np.array(
                    (epsilon, -poisson * epsilon, -poisson * epsilon)
                )
                expected_s = np.tile(
                    (sigma, 0.0, 0.0, 0.0, 0.0, 0.0), (8 * len(study.elements), 1)
                )
                np.testing.assert_allclose(
                    result["displacements"],
                    expected_u,
                    rtol=5e-7,
                    atol=max(np.max(abs(expected_u)) * 1e-10, 1e-30),
                )
                np.testing.assert_allclose(
                    result["stress"],
                    expected_s,
                    rtol=5e-7,
                    atol=max(abs(sigma) * 1e-10, 1e-30),
                )
                expected_energy = force * force * length / (2 * young * width * height)
                np.testing.assert_allclose(
                    result["strain_energy_J"], expected_energy, rtol=5e-7, atol=1e-30
                )
                np.testing.assert_allclose(
                    np.sum(result["reactions"], axis=0),
                    (-force, 0, 0),
                    rtol=5e-7,
                    atol=max(abs(force) * 1e-10, 1e-30),
                )
                saved = json.loads((output / "result.json").read_text())
                self.assertEqual(load_study(output / "study.json"), study)
                self.assertEqual(saved["scientific_status"], "NotAssessed")
                self.assertEqual(saved["input_sha256"], result["input_sha256"])
                self.assertEqual(saved["load_factor"], 1.0)

    def test_parser_rejects_missing_duplicate_nonfinite_and_unsupported_reactions(self):
        study = specimen()
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "run"
            run_static(study, output, executable=CCX)
            raw = (output / "study.dat").read_text()

        def change_value(table, node, column, value, text=raw):
            lines = text.splitlines()
            start = next(
                i for i, line in enumerate(lines) if line.startswith(" " + table)
            )
            row = next(
                i
                for i in range(start + 1, len(lines))
                if lines[i].split()[:1] == [str(node)]
            )
            fields = lines[row].split()
            fields[column] = value
            lines[row] = " " + " ".join(fields)
            return "\n".join(lines) + "\n"

        corruptions = (
            raw.replace("displacements ", "missing ", 1),
            raw + raw,
            change_value("displacements ", 2, 1, "NaN"),
            change_value("displacements ", 1, 1, "1.0"),
            change_value("displacements ", 1, 0, "2"),
            raw.replace("and time  0.1000000E+01", "and time  0.5000000E+00", 1),
            change_value("forces ", 2, 1, "1.0E6"),
        )
        for index, text in enumerate(corruptions):
            self.assertNotEqual(text, raw)
            with self.subTest(corruption=index), self.assertRaises(ValueError):
                parse_dat(text, study)
        # Two opposite perturbations at fixed X DOFs preserve resultant force
        # and work, but introduce a 10 N m couple about Z.
        couple = change_value("forces ", 1, 1, "-150.0")
        couple = change_value("forces ", 3, 1, "-350.0", couple)
        with self.assertRaisesRegex(ValueError, "moments"):
            parse_dat(couple, study)
