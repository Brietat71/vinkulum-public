"""Physical load invariants and the retained three-mesh bending experiment."""

import runpy
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RECIPE = runpy.run_path(ROOT / "ci/studio_c3d8_refinement.py")


class CantileverRefinement(unittest.TestCase):
    def test_uniform_end_traction_and_clamp_are_preserved_under_refinement(self):
        for cells in RECIPE["MESHES"]:
            with self.subTest(cells=cells):
                study = RECIPE["cantilever"](cells)
                nodes = np.asarray(study.nodes)
                forces = np.zeros_like(nodes)
                for node, *force in study.forces:
                    self.assertEqual(nodes[node - 1, 0], 0.12)
                    forces[node - 1] = force
                np.testing.assert_allclose(forces.sum(axis=0), (0, 0, -1), atol=1e-14)
                np.testing.assert_allclose(
                    np.cross(nodes, forces).sum(axis=0),
                    (0, 0.12, 0),
                    rtol=0,
                    atol=1e-14,
                )
                expected_clamp = {
                    (i + 1, axis)
                    for i, xyz in enumerate(nodes)
                    if xyz[0] == 0
                    for axis in (1, 2, 3)
                }
                self.assertEqual(set(study.fixed_dofs), expected_clamp)
                # A bilinear virtual displacement has the same exact face integral
                # for all meshes: the odd y, z and yz terms integrate to zero.
                virtual_uz = (
                    2
                    + 3 * nodes[:, 1]
                    - 7 * nodes[:, 2]
                    + 5 * nodes[:, 1] * nodes[:, 2]
                )
                self.assertAlmostEqual(
                    float(np.dot(forces[:, 2], virtual_uz)), -2, places=13
                )

    def test_retained_results_recompute_without_executing_the_solver(self):
        with patch(
            "subprocess.run",
            side_effect=AssertionError("Archive verification must not launch a solver"),
        ):
            report = RECIPE["verify"](ROOT / "docs/bancs/studio-c3d8-refinement-2026")
        self.assertEqual(len(report["cases"]), 3)
        self.assertEqual(
            [row["free_dofs"] for row in report["cases"]], [360, 1800, 10800]
        )
        self.assertEqual(
            {row["scientific_status"] for row in report["cases"]}, {"NotAssessed"}
        )


if __name__ == "__main__":
    unittest.main()
