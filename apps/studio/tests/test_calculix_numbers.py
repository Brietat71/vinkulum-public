"""Exercise real F20.0 input, geometry admission and archive compatibility."""

import math
import random
import shutil
import struct
import subprocess
import tempfile
import unittest
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

from test_calculix import CCX, specimen
from vinkulum_studio.calculix import (
    StaticStudy,
    load_static_result,
    prepare_run,
    run_static,
)
from vinkulum_studio.calculix_numbers import number_field, transported_number


class DecimalTransport(unittest.TestCase):
    def test_known_truncation_and_full_binary64_range(self):
        rng = random.Random(472391)
        values = [
            -3.552713678800501e-18,
            2.0205265095961919e-15,
            -0.040617999075839344,
            0.0,
            -0.0,
            5e-324,
            -5e-324,
            1.7976931348623157e308,
            -1.7976931348623157e308,
        ]
        for _ in range(10000):
            value = struct.unpack(">d", rng.getrandbits(64).to_bytes(8, "big"))[0]
            if math.isfinite(value):
                values.append(value)
        for value in values:
            field = number_field(value)
            actual = transported_number(value)
            self.assertLessEqual(len(field), 20)
            self.assertTrue(math.isfinite(actual))
            self.assertLessEqual(
                abs(Fraction(value) - Fraction(actual)),
                64 * Fraction(math.ulp(value)),
            )
        # These failures can be repaired without changing their binary64 value.
        for value in values[:3]:
            self.assertEqual(transported_number(value), value)
        self.assertEqual(math.copysign(1.0, transported_number(-0.0)), -1.0)
        self.assertEqual(math.copysign(1.0, transported_number(0.0)), 1.0)

    @unittest.skipUnless(
        shutil.which("gfortran"), "Independent Fortran reader requires gfortran"
    )
    def test_actual_fortran_reader_matches_binary64_values(self):
        # Reproduce the upstream reader's exact field slice and format, using
        # a compiled Fortran runtime rather than our own decoder as the oracle.
        reader = """program reader
use iso_fortran_env, only: real64, int64
implicit none
character(len=132) :: line
real(real64) :: value
integer :: status
do
  read(*, '(A)', iostat=status) line
  if(status /= 0) exit
  read(line(1:20), '(f20.0)', iostat=status) value
  if(status /= 0) stop 1
  write(*, '(Z16.16)') transfer(value, 0_int64)
end do
end program
"""
        rng = random.Random(9371)
        values = [
            -3.552713678800501e-18,
            2.0205265095961919e-15,
            0.0,
            -0.0,
            5e-324,
            -5e-324,
            1.7976931348623157e308,
            -1.7976931348623157e308,
        ]
        for _ in range(10000):
            value = struct.unpack(">d", rng.getrandbits(64).to_bytes(8, "big"))[0]
            if math.isfinite(value):
                values.append(value)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, executable = root / "reader.f90", root / "reader"
            source.write_text(reader)
            subprocess.run(
                [shutil.which("gfortran"), str(source), "-o", str(executable)],
                capture_output=True,
                check=True,
                timeout=30,
            )
            result = subprocess.run(
                [str(executable)],
                input="\n".join(map(number_field, values)) + "\n",
                text=True,
                capture_output=True,
                check=True,
                timeout=30,
            )
        actual = [
            struct.unpack(">d", bytes.fromhex(row))[0]
            for row in result.stdout.splitlines()
        ]
        self.assertEqual(actual, list(map(transported_number, values)))
        self.assertEqual(
            [struct.pack(">d", v) for v in actual],
            [struct.pack(">d", transported_number(v)) for v in values],
        )

    def test_rounding_cannot_bypass_element_admission(self):
        left = -1.2345678901234568e-5
        right = math.nextafter(left, math.inf)
        self.assertEqual(transported_number(left), transported_number(right))
        # Valid, thin tetrahedron whose two distinct x coordinates collapse
        # in F20.0. Reject before creating a run, not after a solver failure.
        study = StaticStudy(
            (
                (left, 0.0, 0.0),
                (right, 0.0, 0.0),
                (left, 1e-10, 0.0),
                (left, 0.0, 1e-10),
            ),
            ((1, 2, 3, 4),),
            tuple((i, a) for i in range(1, 5) for a in (1, 2, 3)),
            (),
            210e9,
            0.3,
            "C3D4",
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "must-not-exist"
            with self.assertRaisesRegex(ValueError, "degenerate"):
                prepare_run(study, root, Path("unused"), "unused", "unused")
            self.assertFalse(root.exists())

    @unittest.skipUnless(CCX, "Requires the separate CalculiX executable")
    def test_real_solver_and_requested_input_archive(self):
        study = specimen(
            length=0.12345678901234567,
            width=0.02,
            height=0.03,
            force=123.45678901234567,
            cells=(2, 2, 2),
        )
        study = replace(
            study,
            nodes=tuple(
                (
                    x - 0.06345678901234567,
                    y - 0.010123456789012345,
                    z - 0.015123456789012345,
                )
                for x, y, z in study.nodes
            ),
        )
        self.assertGreater(study.input_transport["coordinates_m"]["changed_values"], 0)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "run"
            report = run_static(study, root, executable=CCX)
            actual = study.solver_study
            length = max(p[0] for p in actual.nodes) - min(p[0] for p in actual.nodes)
            area = (
                max(p[1] for p in actual.nodes) - min(p[1] for p in actual.nodes)
            ) * (max(p[2] for p in actual.nodes) - min(p[2] for p in actual.nodes))
            force = sum(row[1] for row in actual.forces)
            reference = force**2 * length / (2 * actual.young_pa * area)
            self.assertLess(abs(report["strain_energy_J"] / reference - 1), 5e-6)
            requested, reopened, _ = load_static_result(root)
            self.assertEqual(requested, study)
            self.assertEqual(reopened, report)
            self.assertEqual(report["input_transport"], study.input_transport)

    def test_published_legacy_archives_keep_their_original_deck_contract(self):
        root = Path(__file__).resolve().parents[3]
        for folder in (
            root / "docs/bancs/studio-tetrahedra-060/calculation",
            root / "docs/bancs/studio-050/calculation",
        ):
            if not folder.exists():
                continue
            study, report, _ = load_static_result(folder)
            self.assertIn(report["schema"], (1, 2))
            self.assertEqual(
                (folder / "study.inp").read_text(), study.input_deck(legacy=True)
            )


if __name__ == "__main__":
    unittest.main()
