from fractions import Fraction as F
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from thermal_table import render, source_cell, source_rows


class ThermalTableTests(unittest.TestCase):
    rows = [(F(t), F(2), F(1, 10), F(3), F(1, 5)) for t in (1, 2, 3)]

    def test_boundary_and_interior_knots(self):
        for t, index in [(1, 0), (F(3, 2), 0), (2, 1), (3, 1)]:
            self.assertEqual(source_cell(self.rows, F(t))[0], index)

    def test_no_extrapolation(self):
        for t in [F(99, 100), F(301, 100)]:
            with self.assertRaises(ValueError):
                source_cell(self.rows, t)

    def test_source_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "wrong.dat"
            path.write_bytes(b"1 2 0.1 3 0.2\n")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                source_rows(path)

    def test_physical_assumptions_are_printed(self):
        source = render(self.rows, F(2))
        self.assertIn("h : SaikawaShiraiDoF inputCell gρ gs", source)
        self.assertIn("(19/10)", source)
        self.assertIn("(16/5)", source)


if __name__ == "__main__":
    unittest.main()
