import sys
import unittest
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import sashimi


class ScalarApiTests(unittest.TestCase):
    def test_exact_input(self):
        for value in ["0.315", "63/200", Decimal("0.315"), Fraction(63, 200)]:
            self.assertEqual(sashimi.rational(value), Fraction(63, 200))
        for value in [0.315, True, Decimal("NaN"), "1/0", "1); #eval IO.println 1"]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                sashimi.rational(value)

    def test_reject_invalid_request_before_lean(self):
        for name, params, kwargs in [
            ("hubble", {}, {}),
            ("hubble", {"z": 2, "omega_k": "1/100"}, {}),
            ("growth", {"z": 2}, {"profile": "flrw"}),
            ("wdm-q5", {"alpha": 1, "k": 2}, {"profile": "native-flat"}),
            ("hubble", {"z": 2}, {"bits": True}),
            ("hubble", {"z": 2}, {"tol": "-1"}),
            ("nfw-invert-mean", {"target": 1, "lo": 0, "hi": 5}, {}),
        ]:
            with self.subTest(name=name, params=params), self.assertRaises(ValueError):
                sashimi.evaluate(name, params, **kwargs)

    def test_w_legacy_profiles_use_source_specific_polynomials(self):
        for name, builder in [
            ("hubble", "wLegacyHz"),
            ("rho-crit", "wLegacyHz"),
            ("e2", "wLegacy "),
            ("omega-m", "wLegacyOmegaz"),
        ]:
            op = sashimi.OPERATIONS[name]
            args, profile = sashimi.prepare(op, {"z": 1}, "w-legacy")
            self.assertIn(builder, sashimi.lean_call(op, args, profile, "1", 32, 16, 50))
        values, _ = sashimi.prepare(sashimi.OPERATIONS["hubble"], {"z": 2}, "w-consistent")
        self.assertEqual(values["omega_m"], Fraction("0.27"))

    @patch("sashimi.run", side_effect=["", "RESULT none\n"])
    def test_unresolved_has_no_bounds(self, _run):
        result = sashimi.evaluate("hubble", {"z": 2}, tol=0)
        self.assertEqual(result["status"], "unresolved")
        self.assertIsNone(result["interval"])
        self.assertFalse(result["kernel_replayed"])

    @patch("sashimi.run", side_effect=["", "RESULT 2 1\n"])
    def test_rejects_malformed_lean_output(self, _run):
        with self.assertRaises(RuntimeError):
            sashimi.evaluate("nfw-f", {"x": "0.25"})

    def test_catalogue_has_alpha_units_on_both_sides(self):
        info = sashimi.catalogue()
        self.assertEqual(info["wdm-alpha"]["unit"], "Mpc/h")
        self.assertEqual(info["wdm-q5"]["inputs"]["alpha"], "Mpc")


if __name__ == "__main__":
    unittest.main()
