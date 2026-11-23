from fractions import Fraction
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from spectrum_table import evaluate
from table_replay import check_interval


TABLE = Path(__file__).parent / "data/wmap7-two-knots.dat"


class TableReplay(unittest.TestCase):
    def evaluate(self, **changes):
        params = dict(
            layout="wmap7",
            k=".49",
            tolerance=".01",
            interpolation="linear",
            transfer="cdm",
            placement="continuous",
            bits=40,
            degree=12,
            depth=0,
            replay=True,
        )
        return evaluate(TABLE, **(params | changes))

    def test_source_knots_lookup(self):
        result = self.evaluate()
        self.assertEqual(result["status"], "enclosed")
        self.assertTrue(result["kernel_replayed"])
        # Independent exact affine interpolation in the file's (P,k) units.
        h = Fraction(".7")
        x = Fraction(".49") / h
        a, b = map(Fraction, (".69183055", ".70794535"))
        left, right = map(Fraction, ("133.25732", "126.50208"))
        expected = (left + (right - left) * (x - a) / (b - a)) / h**3
        self.assertLessEqual(Fraction(result["lower"]), expected)
        self.assertLessEqual(expected, Fraction(result["upper"]))

    def test_one_cell_moment(self):
        result = self.evaluate(k=".485", k_max=".49", bits=32, degree=8, tolerance=".1")
        self.assertEqual(result["status"], "enclosed")
        self.assertTrue(result["kernel_replayed"])
        self.assertGreater(Fraction(result["lower"]), 0)

    def test_normalized_point_anchor_replayed(self):
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory) / "canonical.dat"
            table.write_text("1 2\n2 2\n")
            result = evaluate(
                table,
                layout="canonical",
                k="1.5",
                tolerance=8,
                interpolation="linear",
                transfer="cdm",
                placement="continuous",
                sigma8=".8",
                anchor_radius=0,
                bits=32,
                degree=12,
                anchor_depth=1,
                replay=True,
            )
        self.assertTrue(result["kernel_replayed"])
        self.assertGreater(Fraction(result["lower"]), 0)
        self.assertTrue(result["normalization"]["full_table"])

    def test_mismatched_native_endpoints_rejected(self):
        run = subprocess.run

        def corrupt_native_result(command, **kwargs):
            output = run(command, **kwargs)
            if "--run" in command:
                self.assertEqual(output.returncode, 0)
                output.stdout = "enclosed 0 0\n"
            return output

        with patch("spectrum_table.subprocess.run", side_effect=corrupt_native_result):
            with self.assertRaisesRegex(RuntimeError, "Concrete kernel replay failed"):
                self.evaluate()

    def test_configurable_replay_limits(self):
        settings = {
            "PHYSICSLEAN_REPLAY_TIMEOUT": "301",
            "PHYSICSLEAN_REPLAY_MEMORY_MB": "6000",
            "PHYSICSLEAN_REPLAY_HEARTBEATS": "12345678",
        }
        with patch.dict(os.environ, settings), patch("table_replay.check_proof") as runner:
            check_interval("", Fraction(0), Fraction(1), TABLE.parent)
            self.assertIn("maxHeartbeats 12345678", runner.call_args.args[0])
            self.assertEqual(
                runner.call_args.kwargs, {"max_memory_mb": 6000, "timeout_seconds": 301}
            )
        with patch.dict(os.environ, {"PHYSICSLEAN_REPLAY_TIMEOUT": "0"}):
            with self.assertRaisesRegex(ValueError, "positive integer"):
                check_interval("", Fraction(0), Fraction(1), TABLE.parent)

    def test_unresolved_has_no_concrete_replay(self):
        result = self.evaluate(k=".485", k_max=".49", tolerance=0)
        self.assertEqual(result["status"], "unresolved")
        self.assertFalse(result["kernel_replayed"])

    def test_wdm_particle_and_units(self):
        particle = dict(transfer="wdm-q5", wdm_mass_kev="1.5", omega_wdm=".2231")
        result = self.evaluate(**particle)
        self.assertTrue(result["kernel_replayed"])
        self.assertEqual(result["wdm_particle"]["h"], "7/10")
        cdm = self.evaluate(replay=False)
        self.assertLess(Fraction(result["upper"]), Fraction(cdm["lower"]))
        for changes in (
            dict(parameter=".04"),
            dict(transfer="fdm"),
            dict(thermal_h=".71"),
            dict(omega_wdm=None),
            dict(wdm_mass_kev=0),
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.evaluate(**(particle | changes))


if __name__ == "__main__":
    unittest.main()
