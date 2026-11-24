"""Check five-panel comparisons against the published models and exact residuals."""

from contextlib import contextmanager
import csv
from fractions import Fraction
from importlib import import_module
import json
from pathlib import Path
from shutil import copyfile
import sys
from tempfile import TemporaryDirectory
from typing import NamedTuple
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
SOURCES = {
    "legacy",
    "matched",
    "sashimi",
    "colossus_spline",
    "colossus_direct",
    "fixed_band",
    "astropy_growth",
}
PANELS = (
    (None, "value"),
    ("legacy", "absolute"),
    ("legacy", "signed"),
    ("matched", "absolute"),
    ("matched", "signed"),
)


class ExpectedInterval(NamedTuple):
    lo: Fraction
    hi: Fraction
    source: str
    status: str


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
        fig, axes = plt.subplots(len(PANELS), sharex=True)
        try:
            plots.curves_panel(axes, q, selected)
            active = set(next(iter(expected.values())))
            xs = [float(x) for x in expected]
            styles = {}
            for (base, mode), ax in zip(PANELS, axes, strict=True):
                with self.subTest(quantity=q.key, base=base, mode=mode):
                    model_lines = [line for line in ax.lines if line.get_gid() in SOURCES]
                    self.assertCountEqual([line.get_gid() for line in model_lines], active)
                    expected_zeros = {}
                    if base is not None:
                        self.assertEqual(ax.get_yscale(), "log" if mode == "absolute" else "linear")
                        self.assertNotRegex(ax.get_ylabel() + ax.get_title(), "[%％]")
                    for line in model_lines:
                        source = line.get_gid()
                        self.assertTrue(line.get_visible())
                        self.assertNotEqual(line.get_alpha(), 0)
                        style = (line.get_color(), line.get_marker(), line.get_linestyle())
                        self.assertEqual(style, styles.setdefault(source, style))
                        np.testing.assert_array_equal(line.get_xdata(), xs)
                        if mode == "value":
                            target = [float(v[source]) * q.factor for v in expected.values()]
                        else:
                            ratios = [(v[source] - v[base]) / v[base] for v in expected.values()]
                            if mode == "absolute":
                                target = [float(abs(r)) if r else np.nan for r in ratios]
                                zero_xs = [x for x, r in zip(xs, ratios, strict=True) if not r]
                                if zero_xs:
                                    expected_zeros[f"zero-{source}"] = zero_xs
                            else:
                                target = [float(r) for r in ratios]
                        data = np.asarray(line.get_ydata(), dtype=float)
                        np.testing.assert_array_equal(np.isnan(data), np.isnan(target))
                        np.testing.assert_allclose(
                            data,
                            target,
                            rtol=4 * np.finfo(float).eps,
                            atol=0,
                            equal_nan=True,
                            err_msg=f"{q.key}: {source} against {base} ({mode})",
                        )
                    markers = [
                        line for line in ax.lines if (line.get_gid() or "").startswith("zero-")
                    ]
                    self.assertCountEqual(
                        [line.get_gid() for line in markers], expected_zeros.keys()
                    )
                    for marker in markers:
                        self.assertTrue(marker.get_visible())
                        self.assertNotIn(marker.get_marker(), (None, "None", ""))
                        self.assertEqual(marker.get_transform(), ax.get_xaxis_transform())
                        np.testing.assert_array_equal(
                            marker.get_xdata(), expected_zeros[marker.get_gid()]
                        )
            np.testing.assert_array_equal(axes[1].get_ylim(), axes[3].get_ylim())
            np.testing.assert_array_equal(axes[2].get_ylim(), axes[4].get_ylim())
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

    def certificate_values(self, q, data=MATCHED):
        coordinate = "z" if q.group == "background" else "mass_msun"
        old_name = "kernel" if q.group == "conditional" else q.group
        old = document(BASELINE, old_name)
        primary = {
            "F": "kernel-bounds",
            "dS_dM": "derivative-full",
            "sigma": "full",
            "S": "full",
            "M200": "host-bounds",
            "dM_dz": "host-bounds",
        }.get(q.key, "growth")
        legacy_file = f"../eps-plots/{old_name}.json"
        reports = [("legacy", old, legacy_file, False)]
        if q.key == "H":
            reports.append(("matched", old, legacy_file, False))
        else:
            reports.append(("matched", document(data, primary), primary + ".json", False))
            fine_path = data / (primary + "-fine.json")
            if fine_path.exists():
                reports.append(("matched", json.loads(fine_path.read_text()), fine_path.name, True))
        result = {"matched": {}, "legacy": {}}
        for model, cert, source, refined in reports:
            table_checked = cert.get("kernel_replayed", cert.get("kernel_replay") == "passed")
            for x, row in by_coordinate(cert["rows"], coordinate).items():
                box = row.get(q.key)
                if not box or row.get("resolved") is False or row.get("status") == "unresolved":
                    continue
                checked = row.get(
                    "end_to_end_kernel_replayed", row.get("kernel_replayed", table_checked)
                )
                if checked or refined:
                    # The published 7680-cell reports evaluated the proved integrator at runtime.
                    # A row-level full replay, when present, supplies a stronger status.
                    result[model][x] = ExpectedInterval(
                        Fraction(box["lo"]),
                        Fraction(box["hi"]),
                        source,
                        "kernel_checked" if checked else "proved_evaluator_runtime",
                    )
        return result

    def assert_bars(self, ax, gid, wanted):
        segments = {}
        for collection in ax.collections:
            if collection.get_gid() != gid:
                continue
            self.assertTrue(collection.get_visible())
            self.assertNotEqual(collection.get_alpha(), 0)
            for segment in collection.get_segments():
                self.assertEqual(segment.shape, (2, 2))
                x = segment[0, 0]
                self.assertEqual(x, segment[1, 0])
                self.assertNotIn(x, segments, "duplicate interval at one mass/redshift")
                segments[x] = segment[:, 1]
        self.assertEqual(segments.keys(), wanted.keys())
        if wanted:
            np.testing.assert_allclose(
                [segments[x] for x in wanted],
                list(wanted.values()),
                rtol=4 * np.finfo(float).eps,
                atol=0,
            )

    def inspect_intervals(self, q, records, sources, certified):
        selected = [r for r in records if r["quantity"] == q.key]
        fig, axes = plt.subplots(len(PANELS), sharex=True)
        try:
            plots.curves_panel(axes, q, selected)
            for model, expected in certified.items():
                physical = {
                    float(x): sorted((float(i.lo) * q.factor, float(i.hi) * q.factor))
                    for x, i in expected.items()
                }
                self.assert_bars(axes[0], f"interval-value-{model}", physical)
                for row in selected:
                    interval = expected.get(Fraction(row["coordinate"]))
                    ends = None if interval is None else (interval.lo, interval.hi)
                    self.assertEqual(plots.plotting_enclosure(row, model), ends)
                    if model == "matched":
                        displayed = tuple(
                            row[f"matched_plot_{side}_exact"] for side in ("lo", "hi")
                        )
                        if interval is None:
                            self.assertEqual(displayed, ("", ""))
                        else:
                            self.assertEqual(tuple(map(Fraction, displayed)), ends)
                            self.assertEqual(row["matched_plot_enclosure_source"], interval.source)
                            self.assertEqual(row["matched_plot_enclosure_status"], interval.status)
            for (base, mode), ax in zip(PANELS[1:], axes[1:], strict=True):
                plot_min, plot_max = ax.get_ylim()
                absolute = mode == "absolute"
                geometry = "abs" if absolute else "signed"
                for model, bounds in certified.items():
                    with self.subTest(quantity=q.key, base=base, mode=mode, interval=model):
                        wanted, zeros = {}, []
                        for x, interval in bounds.items():
                            lo, hi = interval.lo, interval.hi
                            ref = sources[x][base]
                            a, b = sorted(((lo - ref) / ref, (hi - ref) / ref))
                            nominal = (sources[x][model] - ref) / ref
                            self.assertLessEqual(a, nominal)
                            self.assertLessEqual(nominal, b)
                            if absolute:
                                contains_zero = a <= 0 <= b
                                if contains_zero:
                                    zeros.append(float(x))
                                lower = Fraction(0) if contains_zero else min(abs(a), abs(b))
                                upper = max(abs(a), abs(b))
                                if not upper or float(upper) < plot_min or float(lower) > plot_max:
                                    continue
                                ends = (max(plot_min, float(lower)), min(plot_max, float(upper)))
                            else:
                                ends = (float(a), float(b))
                            wanted[float(x)] = ends

                        self.assert_bars(ax, f"interval-{geometry}-{model}", wanted)
                        indicators = [
                            line for line in ax.lines if line.get_gid() == f"interval-zero-{model}"
                        ]
                        shown = []
                        for marker in indicators:
                            self.assertTrue(marker.get_visible())
                            self.assertEqual(marker.get_transform(), ax.get_xaxis_transform())
                            self.assertNotIn(marker.get_marker(), (None, "", "None"))
                            shown.extend(marker.get_xdata())
                        self.assertCountEqual(shown, zeros)
        finally:
            plt.close(fig)

    def test_certified_interval_glyphs_against_both_models(self):
        for q in (*plots.BACKGROUND, *plots.MASS, *plots.HOST):
            with self.subTest(quantity=q.key):
                self.inspect_intervals(
                    q, self.records, self.source_values(q), self.certificate_values(q)
                )

    def test_singleton_and_unchecked_intervals(self):
        # Exercise a point enclosure and an unchecked proposal without changing the audited inputs.
        q = next(q for q in plots.MASS if q.key == "F")
        sources = self.source_values(q)
        with TemporaryDirectory() as tmp:
            data = Path(tmp)
            for path in MATCHED.glob("*.json"):
                if path.name != "kernel-bounds.json" and not path.name.endswith("-fine.json"):
                    copyfile(path, data / path.name)
            kernel = document(MATCHED, "kernel-bounds")
            first, second = kernel["rows"][:2]
            x = Fraction(first["mass_msun"])
            prediction = sources[x]["matched"]
            self.assertNotEqual(prediction, sources[x]["legacy"])
            first["F"] = {"lo": str(prediction), "hi": str(prediction)}
            second["end_to_end_kernel_replayed"] = False
            (data / "kernel-bounds.json").write_text(json.dumps(kernel))
            rows = plots.comparison_rows(data, BASELINE)
            self.inspect_intervals(q, rows, sources, self.certificate_values(q, data))
            # A missing host refinement must also leave its checked primary history drawable.
            host = plots.HOST[0]
            self.inspect_intervals(
                host, rows, self.source_values(host), self.certificate_values(host, data)
            )

    def test_missing_refined_row_uses_primary_only_at_that_mass(self):
        q = next(q for q in plots.MASS if q.key == "F")
        with TemporaryDirectory() as tmp:
            data = Path(tmp)
            for p in MATCHED.glob("*.json"):
                if p.name not in ("kernel-bounds.json", "kernel-bounds-fine.json"):
                    copyfile(p, data / p.name)
            primary = document(MATCHED, "kernel-bounds")
            primary["rows"][1]["end_to_end_kernel_replayed"] = False
            (data / "kernel-bounds.json").write_text(json.dumps(primary))
            refined = document(MATCHED, "kernel-bounds-fine")
            refined["rows"][0].update(F=None, status="unresolved")
            (data / "kernel-bounds-fine.json").write_text(json.dumps(refined))
            rows = plots.comparison_rows(data, BASELINE)
            expected = self.certificate_values(q, data)
            m0, m1 = [Fraction(r["mass_msun"]) for r in refined["rows"][:2]]
            self.assertEqual(expected["matched"][m0].source, "kernel-bounds.json")
            self.assertEqual(expected["matched"][m1].source, "kernel-bounds-fine.json")
            self.inspect_intervals(q, rows, self.source_values(q), expected)

    def test_csv_keeps_both_integral_meshes_and_the_selected_interval(self):
        primary = by_coordinate(document(MATCHED, "kernel-bounds")["rows"], "mass_msun")
        fine = by_coordinate(document(MATCHED, "kernel-bounds-fine")["rows"], "mass_msun")
        with TemporaryDirectory() as tmp:
            plots.save_tables(Path(tmp), self.records)
            with (Path(tmp) / "matched-samples.csv").open(newline="") as f:
                exported = by_coordinate(
                    (r for r in csv.DictReader(f) if r["quantity"] == "F"), "coordinate"
                )
        self.assertEqual(exported.keys(), primary.keys())
        for x, row in exported.items():
            for prefix, box, source, state in (
                ("matched", primary[x]["F"], "kernel-bounds.json", "kernel_checked"),
                (
                    "matched_fine",
                    fine[x]["F"],
                    "kernel-bounds-fine.json",
                    "proved_evaluator_runtime",
                ),
                (
                    "matched_plot",
                    fine[x]["F"],
                    "kernel-bounds-fine.json",
                    "proved_evaluator_runtime",
                ),
            ):
                self.assertEqual(
                    [row[f"{prefix}_{k}_exact"] for k in ("lo", "hi")], [box["lo"], box["hi"]]
                )
                self.assertEqual(row[prefix + "_enclosure_source"], source)
                self.assertEqual(row[prefix + "_enclosure_status"], state)
            # Populated endpoints alone cannot qualify an otherwise unresolved refinement.
            stale = dict(row, matched_fine_enclosure_status="unresolved")
            self.assertEqual(
                plots.plotting_enclosure(stale, "matched"),
                tuple(Fraction(primary[x]["F"][k]) for k in ("lo", "hi")),
            )

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
        fine = by_coordinate(document(MATCHED, "kernel-bounds-fine")["rows"], "mass_msun")[x]["F"]
        self.assertNotEqual(fine, bound)
        self.assertNotEqual(midpoint(fine), values[x]["matched"])
        expected = tuple(Fraction(fine[k]) for k in ("lo", "hi"))
        self.assertEqual(plots.plotting_enclosure(row, "matched"), expected)
        self.assertEqual(row["matched_fine_enclosure_source"], "kernel-bounds-fine.json")
        self.assertEqual(row["matched_fine_enclosure_status"], "proved_evaluator_runtime")
        self.assertEqual(row["matched_enclosure_status"], "kernel_checked")

    def test_coincident_hubble_references_preserve_binary64_residual(self):
        q = next(q for q in plots.BACKGROUND if q.key == "H")
        sources = self.source_values(q)
        zero = sources[0]
        self.assertEqual(zero["legacy"], zero["matched"])
        self.assertEqual(float(zero["matched"]), float(zero["colossus_spline"]))
        self.assertNotEqual(zero["matched"], zero["colossus_spline"])
        # The reference lines coincide, but the binary64 Colossus residual is nonzero.
        self.inspect_panel(q, self.records, sources)

    def test_all_predictions_equal_gives_only_zero_markers_on_log_axes(self):
        q = next(q for q in plots.BACKGROUND if q.key == "D")
        row = next(
            r for r in self.records if r["quantity"] == "D" and r["coordinate"] == "0"
        ).copy()
        # At z=0, replace the one rounded SASHIMI prediction by the exact normalization D(0)=1.
        row.update(sashimi=1.0, rel_sashimi_vs_legacy=0.0, rel_sashimi_vs_matched=0.0)
        values = dict(self.source_values(q)[0], sashimi=Fraction(1))
        self.assertEqual(set(values.values()), {Fraction(1)})
        self.inspect_panel(q, [row], {Fraction(0): values})

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
            fig, axes = plt.subplots(len(PANELS), sharex=True)
            try:
                with self.assertRaisesRegex(ValueError, "astropy_growth"):
                    rows = plots.comparison_rows(MATCHED, baseline)
                    plots.curves_panel(axes, q, [r for r in rows if r["quantity"] == q.key])
            finally:
                plt.close(fig)


if __name__ == "__main__":
    unittest.main()
