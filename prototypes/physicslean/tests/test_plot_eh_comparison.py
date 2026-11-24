"""Check the plotted series and both residuals against the published physical inputs."""

from contextlib import contextmanager
from fractions import Fraction
from importlib import import_module
import json
from pathlib import Path
from shutil import copyfile
import sys
from tempfile import TemporaryDirectory
import unittest

try:
    import matplotlib
except ModuleNotFoundError as exc:
    if exc.name != "matplotlib":
        raise
    raise unittest.SkipTest("Matplotlib is required for the comparison figure") from exc

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))
plots = import_module("plot_eh_match")


MATCHED = PROJECT / "results" / "eps-eh"
BASELINE = PROJECT / "results" / "eps-plots"


def document(directory, name):
    return json.loads((directory / f"{name}.json").read_text())


def by_coordinate(rows, key):
    return {Fraction(str(row[key])): row for row in rows}


def midpoint(box):
    return (Fraction(box["lo"]) + Fraction(box["hi"])) / 2


@contextmanager
def without_astropy_growth_at(coordinates):
    with TemporaryDirectory() as tmp:
        baseline = Path(tmp)
        for infile in BASELINE.glob("*.json"):
            if infile.name != "external.json":
                copyfile(infile, baseline / infile.name)
        external = document(BASELINE, "external")
        for row in external["astropy"]["background"]:
            if Fraction(str(row["z"])) in coordinates:
                row.pop("D")
        (baseline / "external.json").write_text(json.dumps(external))
        yield baseline


class MatchedComparisonTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = plots.comparison_rows(MATCHED, BASELINE)
        cls.native = document(BASELINE, "sashimi")
        cls.astropy = by_coordinate(document(BASELINE, "external")["astropy"]["background"], "z")
        cls.colossus = document(MATCHED, "colossus")
        cls.growth = by_coordinate(document(MATCHED, "growth")["rows"], "z")
        cls.variance = by_coordinate(document(MATCHED, "lean")["mass_rows"], "mass_msun")
        cls.downstream = document(MATCHED, "downstream")

    def source_values(self, q):
        """Use certified rational midpoints and the exact archived binary64 values."""
        key = q.key
        group = q.group
        coordinate = "z" if group == "background" else "mass_msun"
        old = by_coordinate(
            document(BASELINE, "kernel" if group == "conditional" else group)["rows"], coordinate
        )
        native = by_coordinate(self.native[group], coordinate)
        colossus = by_coordinate(self.colossus[group], coordinate)
        if q in plots.HOST or group == "conditional":
            histories = {
                method: by_coordinate(self.downstream[method][group], coordinate)
                for method in ("matched", "spline", "direct")
            }
        expected = {}
        for x in sorted(old):
            reference = midpoint(old[x][key])
            values = {"legacy": reference, "sashimi": Fraction(native[x][key])}
            if q in plots.HOST or group == "conditional":
                values["matched"] = Fraction(histories["matched"][x][key])
            elif key == "H":
                values["matched"] = reference
            elif group == "background":
                values["matched"] = midpoint(self.growth[x][key])
            else:
                values["matched"] = Fraction(self.variance[x][key])
            for method in ("direct", "spline"):
                if q in plots.HOST:
                    value = histories[method][x][key]
                elif group == "conditional":
                    value = colossus[x][f"F_{method}"]
                else:
                    value = colossus[x][method][key]
                values["colossus_" + method] = Fraction(value)
            if group == "conditional":
                values["fixed_band"] = Fraction(colossus[x]["F_fixed_band"])
            elif "fixed_band" in colossus[x]:
                values["fixed_band"] = Fraction(colossus[x]["fixed_band"][key])
            if group == "background" and key in self.astropy[x]:
                values["astropy_growth"] = Fraction(self.astropy[x][key])
            expected[x] = values
        return expected

    def inspect_panel(self, q, rows, expected):
        selected = [row for row in rows if row["quantity"] == q.key]
        self.assertEqual({Fraction(row["coordinate"]) for row in selected}, set(expected))
        fig, axes = plt.subplots(3, sharex=True)
        try:
            plots.curves_panel(axes, q, selected)
            active = set(next(iter(expected.values())))
            xs = [float(x) for x in expected]
            styles = {}
            for row_number, (ax, base) in enumerate(zip(axes, (None, "legacy", "matched"))):
                with self.subTest(quantity=q.key, row=row_number):
                    model_lines = [line for line in ax.lines if line.get_gid() is not None]
                    self.assertCountEqual([line.get_gid() for line in model_lines], active)
                    for line in model_lines:
                        source = line.get_gid()
                        self.assertTrue(line.get_visible())
                        self.assertNotEqual(line.get_alpha(), 0)
                        style = (line.get_color(), line.get_marker(), line.get_linestyle())
                        self.assertEqual(style, styles.setdefault(source, style))
                        np.testing.assert_array_equal(line.get_xdata(), xs)
                        target = [
                            float(v[source]) * q.factor
                            if base is None
                            else float(100 * (v[source] - v[base]) / v[base])
                            for v in expected.values()
                        ]
                        # Percentage scaling may be rounded before or after float conversion.
                        # Zero absolute tolerance preserves sub-ULP cancellations at the reference.
                        np.testing.assert_allclose(
                            line.get_ydata(),
                            target,
                            rtol=4 * np.finfo(float).eps,
                            atol=0,
                            err_msg=f"{q.key}: {source} against {base}",
                        )
        finally:
            plt.close(fig)
        for row in selected:
            values = expected[Fraction(row["coordinate"])]
            for base in ("legacy", "matched"):
                self.assertEqual(Fraction(row[base + "_reference_exact"]), values[base])
                for source in active:
                    expected_difference = float((values[source] - values[base]) / values[base])
                    self.assertEqual(row[f"rel_{source}_vs_{base}"], expected_difference)

    def test_published_series_and_residuals(self):
        for q in (*plots.MASS, *plots.BACKGROUND[1:], *plots.HOST):
            with self.subTest(quantity=q.key):
                self.inspect_panel(q, self.records, self.source_values(q))

    def test_wide_eh_enclosure_does_not_move_the_central_prediction(self):
        q = next(q for q in plots.MASS if q.key == "F")
        values = self.source_values(q)
        x = min(values)
        bound = by_coordinate(document(MATCHED, "kernel-bounds")["rows"], "mass_msun")[x]["F"]
        row = next(
            r for r in self.records if r["quantity"] == q.key and Fraction(r["coordinate"]) == x
        )
        self.assertEqual(
            (row["matched_lo_exact"], row["matched_hi_exact"]), (bound["lo"], bound["hi"])
        )
        self.assertNotEqual(midpoint(bound), values[x]["matched"])
        self.assertEqual(Fraction(row["matched_reference_exact"]), values[x]["matched"])
        legacy = by_coordinate(document(BASELINE, "kernel")["rows"], "mass_msun")[x]["F"]
        fig, axes = plt.subplots(3, sharex=True)
        try:
            plots.curves_panel(axes, q, [r for r in self.records if r["quantity"] == q.key])
            whiskers = [
                tuple(segment[:, 1])
                for collection in axes[0].collections
                for segment in collection.get_segments()
                if segment.shape == (2, 2) and np.all(segment[:, 0] == float(x))
            ]
            for interval in (legacy, bound):
                self.assertIn(
                    tuple(float(Fraction(interval[end])) for end in ("lo", "hi")), whiskers
                )
        finally:
            plt.close(fig)

    def test_coincident_hubble_references_preserve_binary64_residual(self):
        q = next(q for q in plots.BACKGROUND if q.key == "H")
        sources = self.source_values(q)
        zero = sources[0]
        self.assertEqual(zero["legacy"], zero["matched"])
        self.assertEqual(float(zero["matched"]), float(zero["colossus_spline"]))
        self.assertNotEqual(zero["matched"], zero["colossus_spline"])
        # The reference lines coincide, but the binary64 Colossus residual is nonzero.
        self.inspect_panel(q, self.records, sources)

    def test_astropy_wholly_unavailable_for_a_quantity_is_omitted(self):
        q = next(q for q in plots.BACKGROUND if q.key == "D")
        expected = self.source_values(q)
        for values in expected.values():
            del values["astropy_growth"]
        with without_astropy_growth_at(set(expected)) as baseline:
            self.inspect_panel(q, plots.comparison_rows(MATCHED, baseline), expected)

    def test_incomplete_astropy_grid_is_reported(self):
        q = next(q for q in plots.BACKGROUND if q.key == "D")
        with without_astropy_growth_at({Fraction(1, 4)}) as baseline:
            fig, axes = plt.subplots(3, sharex=True)
            try:
                with self.assertRaisesRegex(ValueError, "astropy_growth"):
                    rows = plots.comparison_rows(MATCHED, baseline)
                    plots.curves_panel(axes, q, [r for r in rows if r["quantity"] == q.key])
            finally:
                plt.close(fig)


if __name__ == "__main__":
    unittest.main()
