"""Independent archived FEM reference and physically wrong finite counterexamples."""

import copy
import json
import unittest
from pathlib import Path
from zipfile import ZipFile

from check_freecad_static_reference import check_archive, check_case

ARCHIVE = (
    Path(__file__).resolve().parents[1] / "docs/bancs/freecad-static-2026/record.zip"
)


class AffineReference(unittest.TestCase):
    def setUp(self):
        with ZipFile(ARCHIVE) as archive:
            self.request = json.loads(archive.read("case-1/request.json"))
            self.preview = json.loads(archive.read("case-1/preview.json"))

    def test_retained_rotated_and_original_cases(self):
        result = check_archive(ARCHIVE)
        self.assertEqual(result["status"], "passed")
        self.assertEqual(len(result["cases"]), 2)

    def test_wrong_displacement_direction_is_not_an_affine_solution(self):
        changed = copy.deepcopy(self.preview)
        changed["displacements_m"] = [
            [-v for v in row] for row in changed["displacements_m"]
        ]
        with self.assertRaisesRegex(ValueError, "Displacement"):
            check_case(self.request, changed, 1)

    def test_correct_displacements_do_not_excuse_wrong_energy(self):
        changed = {
            **self.preview,
            "strain_energy_J": 1.01 * self.preview["strain_energy_J"],
        }
        with self.assertRaisesRegex(ValueError, "Energy"):
            check_case(self.request, changed, 1)

    def test_wrong_reaction_sign_is_not_equilibrium(self):
        changed = copy.deepcopy(self.preview)
        changed["reactions_N"] = [[-v for v in row] for row in changed["reactions_N"]]
        with self.assertRaisesRegex(ValueError, "reaction"):
            check_case(self.request, changed, 1)


if __name__ == "__main__":
    unittest.main()
