import copy
import sys
import unittest
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from hz_cert import make_certificate, verify_certificate, render_lean


class CertificateTests(unittest.TestCase):
    def test_z1(self):
        c = make_certificate(Fraction(70), Fraction(3, 10), Fraction(1), 8)
        self.assertEqual(c["h_squared"], "15190")
        self.assertEqual(c["lower"], "12324771803/100000000")
        self.assertEqual(c["upper"], "3081192951/25000000")
        self.assertTrue(verify_certificate(c))

    def test_z0_exact(self):
        c = make_certificate(Fraction(70), Fraction(3, 10), Fraction(0), 10)
        self.assertEqual(c["lower"], "70")
        self.assertEqual(c["upper"], "70")

    def test_lambda_only(self):
        c = make_certificate(Fraction(70), Fraction(0), Fraction(11), 6)
        self.assertEqual(c["lower"], "70")
        self.assertEqual(c["upper"], "70")

    def test_nonzero_radicand(self):
        c = make_certificate(Fraction('67.4'), Fraction('0.315'), Fraction('2.5'), 9)
        self.assertTrue(verify_certificate(c))

    def test_reject_bad_input(self):
        for h, om, z in [(-1, Fraction(3,10), 1), (70, -1, 1), (70, 2, 1), (70, 0, -1)]:
            with self.subTest(h=h, om=om, z=z):
                with self.assertRaises(ValueError):
                    make_certificate(Fraction(h), Fraction(om), Fraction(z), 8)

    def test_tamper_upper_detected(self):
        c = make_certificate(Fraction(70), Fraction(3,10), Fraction(1), 8)
        broken = copy.deepcopy(c)
        broken["upper"] = broken["lower"]
        with self.assertRaises(ValueError):
            verify_certificate(broken)

    def test_tamper_input_detected(self):
        c = make_certificate(Fraction(70), Fraction(3,10), Fraction(1), 8)
        broken = copy.deepcopy(c)
        broken["z"] = "2"
        with self.assertRaises(ValueError):
            verify_certificate(broken)

    def test_lean_contains_physlib_bridge(self):
        c = make_certificate(Fraction(70), Fraction(3,10), Fraction(1), 8)
        source = render_lean(c)
        self.assertIn("import HzCertificate", source)
        self.assertIn("HzCertificate.modelH_bounds", source)
        self.assertIn("#print axioms certified_interval", source)

if __name__ == "__main__":
    unittest.main()
