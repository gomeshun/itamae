import sys
import tempfile
import unittest
from fractions import Fraction
from math import pi
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from spectrum_table import evaluate, load_table, number


class SpectrumTables(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = Path(tmp.name) / "input.dat"
        self.path.write_text(".04\n.27\n.5\n.95\n.82\n2 1\n3 2\n")

    def test_exact_h_and_columns(self):
        rows, metadata = load_table(self.path, "wmap7", None)
        self.assertEqual(rows, [(Fraction(".5"), 16), (1, 24)])
        self.assertEqual(metadata["header"]["h"], "1/2")
        with self.assertRaises(ValueError):
            load_table(self.path, "wmap7", ".7")

    def test_invalid_units_shape_and_numbers(self):
        for layout, h in [("canonical", Fraction(1)), ("other", None), ("camb", None)]:
            with self.subTest(layout=layout), self.assertRaises(ValueError):
                load_table(self.path, layout, h)
        for value in [False, 0.07, "1/0", "nan"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                number(value)

    def test_no_extrapolation_and_no_silent_wdm_convention(self):
        args = dict(
            layout="wmap7", k=0, interpolation="linear", placement="continuous", tolerance=0
        )
        for model, parameter in [("cdm", None), ("wdm", "0.07"), ("wdm-q5", None)]:
            with self.subTest(model=model), self.assertRaises(ValueError):
                evaluate(self.path, transfer=model, parameter=parameter, **args)

    def test_sigma_normalization_on_separate_band(self):
        self.path.write_text("1 3\n2 3\n3 3\n")
        args = dict(
            layout="canonical",
            k=3,
            tolerance="1/5",
            interpolation="linear",
            transfer="cdm",
            placement="continuous",
            sigma8="1/2",
            anchor_radius=0,
            anchor_k_min=1,
            anchor_k_max=2,
            anchor_depth=4,
            bits=32,
            degree=8,
        )
        result = evaluate(self.path, replay=True, **args)
        self.assertEqual(result["status"], "enclosed")
        self.assertTrue(result["kernel_replayed"])
        self.assertEqual(result["normalization"]["band_Mpc_inverse"], ["1", "2"])
        self.assertFalse(result["normalization"]["full_table"])
        # V = integral_1^2 3 k²/(2 pi²) dk = 7/(2 pi²).
        expected = 3 * pi**2 / 14
        self.assertLessEqual(float(Fraction(result["lower"])), expected)
        self.assertLessEqual(expected, float(Fraction(result["upper"])))
        self.assertLessEqual(Fraction(result["width"]), Fraction(1, 5))
        self.path.write_text("1 0\n2 0\n3 3\n")
        self.assertEqual(evaluate(self.path, **args)["status"], "unresolved")
