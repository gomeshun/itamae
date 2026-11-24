"""Plot the matched spectrum against Colossus on the published SASHIMI-C EPS grid."""

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
from fractions import Fraction as Q
import gzip
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator, MaxNLocator, NullLocator
import numpy as np

if __package__:
    from .plot_eh_bounds import error_budget, intervals, save_and_plot, point_proofs
    from .plot_eps import BACKGROUND, MASS, HOST, INK, MUTED, GRID, style
else:
    from plot_eh_bounds import error_budget, intervals, save_and_plot, point_proofs
    from plot_eps import BACKGROUND, MASS, HOST, INK, MUTED, GRID, style


BASES = ("legacy", "matched")
KERNEL_CHECKED = "kernel_checked"
EVALUATOR_RUNTIME = "proved_evaluator_runtime"
FINE_CELLS = 7680
SERIES_KEYS = (
    *BASES,
    "sashimi",
    "colossus_spline",
    "colossus_direct",
    "fixed_band",
    "astropy_growth",
)


def point_value(row, key):
    """Preserve rational reference values and the actual binary64 bits of implementations."""
    if key in BASES:
        return Q(row[f"{key}_reference_exact"])
    value = row[key]
    return Q(float(value)) if isinstance(value, str) else Q(value)


def relative_difference(row, key, base):
    """Dimensionless signed difference (value − reference) / reference."""
    if base not in BASES:
        raise ValueError(f"Unknown reference {base!r}")
    ref = point_value(row, base)
    if not ref:
        raise ValueError(f"Cannot form a relative difference with zero {base} reference")
    return (point_value(row, key) - ref) / ref


def active_series(rows):
    """Return complete series for a single quantity. Reject partially joined coordinates."""
    if not rows:
        return ()
    result = []
    for key in SERIES_KEYS:
        present = [r[key] is not None and r[key] != "" for r in rows]
        if all(present):
            result.append(key)
        elif any(present):
            raise ValueError(f"Incomplete {key} series for {rows[0]['quantity']}")
    return tuple(result)


def load(p, name):
    return json.loads((p / (name + ".json")).read_text())


def keyed(rows, field):
    rs = {Q(str(r[field])): r for r in rows}
    if len(rs) != len(rows):
        raise ValueError(f"Duplicate coordinate in {field}")
    return rs


def midpoint(box):
    lo, hi = Q(box["lo"]), Q(box["hi"])
    if lo > hi:
        raise ValueError("Inverted rational interval")
    return (lo + hi) / 2


def validate_models(variance, growth, reference):
    c = variance["cosmology"]
    ref = reference["metadata"]["parameters"]
    paired = [
        (100 * Q(str(c["h"])), ref["H0"]),
        (c["omega_m"], ref["Om0"]),
        (c["omega_b"], ref["Ob0"]),
        (c["ns"], ref["ns"]),
        (c["sigma8"], ref["sigma8"]),
        (c["Tcmb_K"], ref["Tcmb0"]),
        (growth["parameters"]["omega_m"], ref["Om0"]),
        (growth["parameters"]["omega_lambda"], 1 - Q(str(ref["Om0"]))),
        (growth["parameters"]["radiation"], 0),
        (
            c["rho_crit0_coefficient_Msun_h2_Mpc3"],
            reference["metadata"]["rho_crit_0_msun_h2_per_mpc3"],
        ),
    ]
    if any(Q(str(a)) != Q(str(b)) for a, b in paired) or not ref["flat"] or ref["relspecies"]:
        raise ValueError("Incompatible background or spectrum parameters")
    k_band = variance["quadrature"]["k_band_h_per_Mpc"]
    k_ref = reference["metadata"]["matched_bands"]["fixed_band"]["anchor_k_h_mpc"]
    if k_band != k_ref or variance["normalization"]["R8_Mpc_h"] != 8:
        raise ValueError("Incompatible integration band or normalization radius")


def join(data, baseline):
    var = load(data, "lean")
    bg = load(data, "growth")
    if bg["kernel_replayed"] is not True:
        raise ValueError("Expected the checked growth table")
    col = load(data, "colossus")
    validate_models(var, bg, col)
    down = load(data, "downstream")
    grid = load(baseline, "grid")
    mass = set(map(Q, grid["masses_msun"]))
    zs = set(map(Q, grid["redshifts"]))
    native = load(baseline, "sashimi")
    astro = keyed(load(baseline, "external")["astropy"]["background"], "z")
    old_files = {"background": "background", "variance": "variance", "conditional": "kernel"}
    old = {group: load(baseline, name) for group, name in old_files.items()}
    # The H(z) certificate is shared only when the two flat matter + Λ backgrounds agree.
    legacy_bg = old["background"]["parameters"]
    if any(Q(str(legacy_bg[k])) != Q(str(var["cosmology"][k])) for k in ("h", "omega_m")):
        raise ValueError("The two models have different H(z) parameters")
    if (
        old["background"]["kernel_replayed"] is not True
        or old["variance"]["kernel_replay"] != "passed"
        or old["conditional"]["kernel_replayed"] is not True
    ):
        raise ValueError("Expected kernel-checked SASHIMI-C background, variance and EPS intervals")
    records = []
    for q in BACKGROUND + MASS + HOST:
        t = q.group
        is_bg = t == "background"
        fld = "z" if is_bg else "mass_msun"
        coords = zs if is_bg else (mass - {10**12} if t == "conditional" else mass)
        orig = keyed(old[t]["rows"], fld)
        sashimi = keyed(native[t], fld)
        comp = keyed(col[t], fld)
        primary_rows = (
            down["matched"][t]
            if q in HOST or t == "conditional"
            else bg["rows"]
            if is_bg
            else var["mass_rows"]
        )
        new = orig if q.key == "H" else keyed(primary_rows, fld)
        if any(set(x) != coords for x in (orig, sashimi, comp, new)):
            raise ValueError(f"Grid mismatch for {q.key}")
        host_refs = {m: keyed(down[m][t], fld) for m in ("spline", "direct")} if q in HOST else {}
        for x in sorted(coords):
            old_box = orig[x][q.key]
            old_value = midpoint(old_box)
            val = new[x][q.key]
            matched = midpoint(val) if isinstance(val, dict) else Q(val)
            verified_lo = val["lo"] if isinstance(val, dict) else ""
            verified_hi = val["hi"] if isinstance(val, dict) else ""
            cr = comp[x]
            if q in HOST:
                rs, rd = (host_refs[m][x][q.key] for m in ("spline", "direct"))
                rf = None
            elif t == "conditional":
                rs, rd, rf = cr["F_spline"], cr["F_direct"], cr.get("F_fixed_band")
                for mode, observed in (("spline", rs), ("direct", rd)):
                    independent = keyed(down[mode][t], fld)[x]["F"]
                    if not np.isclose(observed, independent, rtol=1e-12, atol=0):
                        raise ValueError(f"Inconsistent conditional inputs at {x}")
            else:
                rs, rd = (cr[m][q.key] for m in ("spline", "direct"))
                rf = cr.get("fixed_band", {}).get(q.key)
            row = dict(
                quantity=q.key,
                coordinate=str(x),
                unit_coordinate="z" if is_bg else "Msun",
                legacy=float(old_value),
                legacy_reference_exact=str(old_value),
                legacy_reference_kind="rational_midpoint",
                legacy_lo_exact=old_box["lo"],
                legacy_hi_exact=old_box["hi"],
                legacy_enclosure_source=f"../eps-plots/{old_files[t]}.json",
                legacy_enclosure_status="kernel_checked",
                sashimi=sashimi[x][q.key],
                astropy_growth=astro[x].get(q.key) if is_bg else None,
                matched=float(matched),
                matched_reference_exact=str(matched),
                matched_reference_kind="rational_midpoint" if isinstance(val, dict) else "binary64",
                matched_lo_exact=verified_lo,
                matched_hi_exact=verified_hi,
                colossus_spline=rs,
                colossus_direct=rd,
                fixed_band=rf,
            )
            for base in BASES:
                for key in SERIES_KEYS:
                    row[f"rel_{key}_vs_{base}"] = (
                        float(relative_difference(row, key, base)) if row[key] is not None else None
                    )
            for name, key, base in (
                ("before", "colossus_spline", "legacy"),
                ("after", "colossus_spline", "matched"),
                ("direct", "colossus_direct", "matched"),
                ("fixed", "fixed_band", "matched"),
            ):
                row["rel_" + name] = row[f"rel_{key}_vs_{base}"]
            records.append(row)
    return records


def load_tight_reports(data, baseline_reports):
    tight_reports = {}
    for stem in ("full", "kernel-bounds", "derivative-full", "host-bounds"):
        name = stem + "-tight"
        if (data / (name + ".json")).exists():
            tight_reports[name] = load(data, name)
    if "full-tight" in tight_reports:
        report = tight_reports["full-tight"]
        for filename, opener in (("sweep-tight.json.gz", gzip.open), ("sweep-tight.json", open)):
            finite = data / filename
            if not finite.is_file():
                continue
            with opener(finite, "rb") as f:
                if sha256(f.read()).hexdigest() != report["finite_band_source_sha256"]:
                    raise ValueError(f"full-tight inputs differ from {filename}")
        mesh = report["mesh"]
        reference_name = "full-fine" if mesh["cells"] == FINE_CELLS else "full"
        reference = baseline_reports.get(reference_name)
        identity = report["same_grid"]
        if identity["matched_cells"] != mesh["cells"] or (
            reference is not None
            and identity["reference_sweep_sha256"] != reference["finite_band_source_sha256"]
        ):
            raise ValueError("Fixed-cell report does not match the original integration mesh")
    for stem in ("host-bounds-tight", "derivative-full-tight"):
        if stem not in tight_reports:
            continue
        checks = tight_reports[stem]["inputs_sha256"]
        if checks["sweep"] != tight_reports["full-tight"]["finite_band_source_sha256"]:
            raise ValueError(f"{stem} and full-tight use different finite-band sweeps")
        input_key, filename = (
            ("growth", "growth.json")
            if stem == "host-bounds-tight"
            else ("finite_derivative", "derivative-tight.json")
        )
        if checks[input_key] != sha256((data / filename).read_bytes()).hexdigest():
            raise ValueError(f"{stem} inputs differ from {filename}")
    return tight_reports


def _comparison_data(data, baseline):
    rows = join(data, baseline)
    bounds = {k: load(data, k) for k in ("full", "derivative", "kernel-bounds")}
    for stem in (
        "full-fine",
        "derivative-fine",
        "derivative-full",
        "derivative-full-fine",
        "kernel-bounds-fine",
        "kernel-derivative",
    ):
        if (data / (stem + ".json")).exists():
            bounds[stem] = load(data, stem)
    tight_reports = load_tight_reports(data, bounds)
    for suffix in ("", "-fine", "-tight"):
        name = "kernel-bounds" + suffix
        report = (tight_reports if suffix == "-tight" else bounds).get(name)
        if report is not None:
            for k, expected in report["inputs_sha256"].items():
                stems = {"full": "full" + suffix, "growth": "growth", "paired": "paired-children"}
                file_name = stems[k] + ".json"
                if sha256((data / file_name).read_bytes()).hexdigest() != expected:
                    raise ValueError(f"{name} inputs differ from {file_name}")
    mass_boxes = intervals(rows, bounds)
    mass_lookup = {(r["quantity"], r["mass_msun"]): r for r in mass_boxes}
    host_report = load(data, "host-bounds")
    host_rows = keyed(host_report["rows"], "z")
    fine_reports = dict(bounds)
    if (data / "host-bounds-fine.json").exists():
        fine_reports["host-bounds-fine"] = load(data, "host-bounds-fine")
    refined_tables = {
        "fine": refined_inputs(fine_reports, "fine"),
        "tight": refined_inputs(tight_reports, "tight"),
    }
    host_keys = {q.key for q in HOST}
    for r in rows:
        key = r["quantity"]
        coordinate = Q(r["coordinate"])
        result = mass_lookup.get((key, r["coordinate"]))
        if result is not None:
            r["matched_lo_exact"], r["matched_hi_exact"] = result["lo"], result["hi"]
            source, checked = result["inputs"], result["end_to_end_kernel_replay"]
        elif key in host_keys:
            host = host_rows[coordinate]
            box = host[key]
            r["matched_lo_exact"], r["matched_hi_exact"] = (
                (box["lo"], box["hi"]) if box else ("", "")
            )
            source, checked = (
                "host-bounds.json",
                host.get("kernel_replayed", host_report["kernel_replayed"]),
            )
        else:
            source, checked = (
                ("../eps-plots/background.json" if key == "H" else "growth.json"),
                True,
            )
        r["enclosure_source"] = source
        r["enclosure_status"] = (
            (KERNEL_CHECKED if checked else EVALUATOR_RUNTIME)
            if r["matched_lo_exact"]
            else "unresolved"
        )
        r["matched_enclosure_source"] = source
        r["matched_enclosure_status"] = r["enclosure_status"]
        for variant, tables in refined_tables.items():
            name = f"{Path(source).stem}-{variant}"
            sample = tables.get(name, {}).get(coordinate, {})
            box = sample.get(key)
            prefix = "matched_" + variant
            available = box is not None and sample.get("status") not in (
                "unresolved",
                "unavailable",
            )
            r[prefix + "_lo_exact"] = box["lo"] if available else ""
            r[prefix + "_hi_exact"] = box["hi"] if available else ""
            r[prefix + "_enclosure_source"] = name + ".json" if available else ""
            r[prefix + "_enclosure_status"] = "unresolved"
            if available:
                closed = sample.get("end_to_end_kernel_replayed", sample.get("kernel_replayed"))
                r[prefix + "_enclosure_status"] = KERNEL_CHECKED if closed else EVALUATOR_RUNTIME
        tighter_source = Path(r["matched_tight_enclosure_source"]).stem
        r["matched_tight_cells"] = (
            tight_reports[tighter_source]["mesh"]["cells"] if tighter_source else ""
        )
        # Keep the nominal numerical line even when an integral was unresolved.
        # Every published bound is checked against the exact reference used for this chart.
        for lo, hi in (
            ("matched_lo_exact", "matched_hi_exact"),
            ("matched_fine_lo_exact", "matched_fine_hi_exact"),
            ("matched_tight_lo_exact", "matched_tight_hi_exact"),
        ):
            if bool(r[lo]) != bool(r[hi]):
                raise ValueError(f"Incomplete bounds for {key} at {r['coordinate']}")
            if r[lo] and not Q(r[lo]) <= point_value(r, "matched") <= Q(r[hi]):
                raise ValueError(
                    f"Matched reference outside {lo}/{hi} for {key} at {r['coordinate']}"
                )
        selected = plot_enclosure_prefix(r, "matched")
        for field in ("lo_exact", "hi_exact", "enclosure_source", "enclosure_status"):
            r["matched_plot_" + field] = (
                r[selected + "_" + field]
                if selected
                else ("unresolved" if field == "enclosure_status" else "")
            )
        if result is not None:
            for side in ("lo", "hi"):
                result["tight_" + side] = r[f"matched_tight_{side}_exact"]
            result["plot_source"] = r["matched_plot_enclosure_source"]
    return rows, bounds, mass_boxes


def refined_inputs(reports, variant):
    """Recognize evaluated mesh reports whose transfer to each physical quantity was checked."""
    result = {}
    for stem in ("full", "derivative-full", "kernel-bounds", "host-bounds"):
        name = stem + "-" + variant
        if name not in reports:
            continue
        doc = reports[name]
        mesh = doc.get("mesh", {})
        if variant == "fine" and mesh.get("cells", doc.get("cells")) != FINE_CELLS:
            continue
        if variant == "tight" and (
            mesh.get("enclosure_method") != "fixed_cell_enclosure" or not mesh.get("cells")
        ):
            continue
        rows = doc["rows"]
        coordinate = "mass_msun"
        if stem == "full":
            if not (doc.get("propagation_kernel_replayed") and doc.get("radius_kernel_replayed")):
                continue
        elif stem == "derivative-full":
            rows = [
                r
                for r in rows
                if r.get("conditional_implication_kernel_checked") or r.get("kernel_replayed")
            ]
        elif stem == "kernel-bounds":
            if not doc.get("kernel_replayed"):
                continue
        else:
            if not doc.get("witness"):
                continue
            coordinate = "z"
            rows = [r for r in rows if r.get("resolved")]
        result[name] = keyed(rows, coordinate)
    return result


def plot_enclosure_prefix(row, source):
    """Use fixed-cell tightening only when it improves both-sided containment on this row."""
    selected = source if row[source + "_enclosure_status"] == KERNEL_CHECKED else None
    if source != "matched":
        return selected
    for prefix in ("matched_fine", "matched_tight"):
        if row.get(prefix + "_enclosure_status") not in (KERNEL_CHECKED, EVALUATOR_RUNTIME):
            continue
        if not row.get(prefix + "_lo_exact") or not row.get(prefix + "_hi_exact"):
            continue
        if selected and prefix == "matched_tight":
            lo, hi = Q(row[prefix + "_lo_exact"]), Q(row[prefix + "_hi_exact"])
            old_lo, old_hi = Q(row[selected + "_lo_exact"]), Q(row[selected + "_hi_exact"])
            if not (old_lo <= lo <= hi <= old_hi and (lo > old_lo or hi < old_hi)):
                continue
        selected = prefix
    return selected


def plotting_enclosure(row, source):
    selected = plot_enclosure_prefix(row, source)
    if selected is None:
        return None
    return Q(row[selected + "_lo_exact"]), Q(row[selected + "_hi_exact"])


def comparison_rows(data, baseline):
    """Join all coordinates, their exact reference values and available enclosures."""
    return _comparison_data(data, baseline)[0]


def save_tables(out, rows):
    with (out / "matched-samples.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    by_quantity = {}
    for key in dict.fromkeys(r["quantity"] for r in rows):
        subset = [r for r in rows if r["quantity"] == key]
        enclosed = [r for r in subset if r["matched_lo_exact"]]
        checked = [r for r in enclosed if r["enclosure_status"] == "kernel_checked"]

        def inside(r, field):
            return Q(r["matched_lo_exact"]) <= point_value(r, field) <= Q(r["matched_hi_exact"])

        summary = {
            "points": len(subset),
            "enclosures": {
                "returned": len(enclosed),
                "end_to_end_kernel_checked": len(checked),
                "plotted_value_inside": sum(inside(r, "matched") for r in enclosed),
                "colossus_spline_inside": sum(inside(r, "colossus_spline") for r in enclosed),
            },
            "plotted_interval_sources": dict(
                Counter(r["matched_plot_enclosure_source"] for r in subset)
            ),
        }
        candidates = [
            r for r in subset if r["matched_fine_lo_exact"] and r["matched_tight_lo_exact"]
        ]
        ratios = []
        for r in candidates:
            fine_width = Q(r["matched_fine_hi_exact"]) - Q(r["matched_fine_lo_exact"])
            if fine_width:
                width = Q(r["matched_plot_hi_exact"]) - Q(r["matched_plot_lo_exact"])
                ratios.append(float(width / fine_width))
        if ratios:
            summary["plotted_to_original_fine_width_ratio"] = {
                "candidate_points": len(candidates),
                "strictly_narrower": sum(
                    r["matched_plot_enclosure_source"] == r["matched_tight_enclosure_source"]
                    for r in candidates
                ),
                "nonzero_fine_widths": len(ratios),
                "minimum": min(ratios),
                "median": float(np.median(ratios)),
                "maximum": max(ratios),
            }
        summary["references"] = {}
        for base in BASES:
            lo, hi = base + "_lo_exact", base + "_hi_exact"
            available = [r for r in subset if r[lo] and r[hi]]
            residuals = {}
            for series in active_series(subset):
                field = f"rel_{series}_vs_{base}"
                worst = max(subset, key=lambda r: abs(r[field]))
                residuals[series] = {
                    "coordinate": worst["coordinate"],
                    "signed_relative_difference": worst[field],
                }
            summary["references"][base] = {
                "reference_kind": subset[0][base + "_reference_kind"],
                "enclosures": len(available),
                "reference_inside": sum(
                    Q(r[lo]) <= point_value(r, base) <= Q(r[hi]) for r in available
                ),
                "residuals": residuals,
            }
        for name in ("before", "after", "direct", "fixed"):
            k = "rel_" + name
            available = [r for r in subset if r[k] is not None]
            if available:
                worst = max(available, key=lambda r: abs(r[k]))
                summary[name] = {
                    "coordinate": worst["coordinate"],
                    "signed_relative_difference": worst[k],
                }
        by_quantity[key] = summary
    (out / "matched-summary.json").write_text(json.dumps(by_quantity, indent=2) + "\n")


SERIES_STYLE = {
    "legacy": ("Lean · SASHIMI モデル", "#3d497e", "-", "s"),
    "matched": ("Lean · Colossus モデル", "#166d65", "-", "^"),
    "sashimi": ("SASHIMI-C", "#bd5b1b", (0, (5, 2)), "x"),
    "colossus_spline": ("Colossus · 既定補間", "#823e78", (0, (6, 2)), "o"),
    "colossus_direct": ("Colossus · 直接計算", "#916313", (0, (2, 2)), "D"),
    "fixed_band": ("Colossus · 同一帯域で再積分", "#246fa0", (0, (5, 2, 1, 2)), "+"),
    "astropy_growth": ("Astropy（成長量は積分）", "#555b34", (0, (2, 2, 5, 2)), "v"),
}


def line_style(key):
    label, color, ls, marker = SERIES_STYLE[key]
    return dict(
        label=label,
        color=color,
        ls=ls,
        marker=marker,
        mfc="none",
        mew=1.05,
        ms=4.8,
        lw=1.65 if key in BASES else 1.1,
    )


@dataclass
class ResidualIntervals:
    x: np.ndarray
    lo: np.ndarray
    hi: np.ndarray

    @property
    def magnitude_lo(self):
        return np.maximum(0, np.maximum(self.lo, -self.hi))

    @property
    def magnitude_hi(self):
        return np.maximum(np.abs(self.lo), np.abs(self.hi))


def ratio_for_plot(r):
    """Refuse to display a nonzero exact ratio as a zero due to float conversion."""
    y = float(r)
    if not np.isfinite(y) or (y == 0 and r != 0):
        raise ValueError(f"Relative difference exceeds plotting precision: {r}")
    return y


def residual_intervals(rows, source, base):
    """Map each checked numerator interval to the residual of its mathematical value."""
    xs, lower, upper = [], [], []
    for row in rows:
        enclosure = plotting_enclosure(row, source)
        if enclosure is None:
            continue
        b = point_value(row, base)
        lo, hi = sorted((end - b) / b for end in enclosure)
        xs.append(float(Q(row["coordinate"])))
        lower.append(ratio_for_plot(lo))
        upper.append(ratio_for_plot(hi))
    return ResidualIntervals(np.array(xs), np.array(lower), np.array(upper))


def interval_bars(ax, xs, lo, hi, source, gid, *, log_floor=None):
    """Draw checked endpoints at sampled x only. A log boundary is never an endpoint."""
    if log_floor is not None:
        visible = hi >= log_floor
        xs, lo, hi = xs[visible], lo[visible], hi[visible]
    if not len(xs):
        return
    color = SERIES_STYLE[source][1]
    # A slim stroke over a wider one keeps identical intervals from masking each other.
    width, cap, layer = (0.9, 2.7, 1.35) if source == "legacy" else (1.8, 3.8, 1.20)
    drawn_lo = lo if log_floor is None else np.maximum(lo, log_floor)
    ax.vlines(xs, drawn_lo, hi, color=color, lw=width, alpha=0.29, zorder=layer, gid=gid)
    real_lo = np.ones(len(lo), dtype=bool) if log_floor is None else lo >= log_floor
    ax.plot(
        np.r_[xs[real_lo], xs],
        np.r_[lo[real_lo], hi],
        ls="none",
        marker="_",
        ms=cap,
        mew=0.9,
        color=color,
        alpha=0.65,
        zorder=layer,
    )


def exact_zero_strip(ax, zeros, boxes):
    """Separate exact equality of point values from membership of zero in a certificate."""
    for bottom, height, label, background in (
        (-0.20, 0.17, "一致", "#f1f3f5"),
        (-0.38, 0.15, "0 を含む", "#eaf3f1"),
    ):
        ax.add_patch(
            Rectangle(
                (0, bottom),
                1,
                height,
                transform=ax.transAxes,
                fc=background,
                ec="none",
                clip_on=False,
            )
        )
        ax.text(
            -0.017,
            bottom + height / 2,
            label,
            transform=ax.transAxes,
            ha="right",
            va="center",
            color=MUTED,
            fontsize=8,
        )
    for i, (key, xs) in enumerate(zeros.items()):
        lane = -0.20 + (i + 0.5) * 0.17 / len(zeros)
        _, color, _, marker = SERIES_STYLE[key]
        ax.plot(
            xs,
            np.full(len(xs), lane),
            transform=ax.get_xaxis_transform(),
            clip_on=False,
            scaley=False,
            ls="none",
            marker=marker,
            color=color,
            mfc="none",
            ms=3,
            mew=0.65,
            gid=f"zero-{key}",
        )
    for key, lane in zip(BASES, (-0.269, -0.343), strict=True):
        box = boxes[key]
        xs = box.x[(box.lo <= 0) & (box.hi >= 0)]
        if len(xs):
            ax.plot(
                xs,
                np.full(len(xs), lane),
                transform=ax.get_xaxis_transform(),
                clip_on=False,
                scaley=False,
                ls="none",
                marker="|",
                color=SERIES_STYLE[key][1],
                ms=3,
                mew=1.0,
                gid=f"interval-zero-{key}",
            )


def set_residual_scales(log_axes, linear_axes, positive, linear_peak):
    """Use the same scale for a given quantity across its two choices of reference."""
    if positive:
        log_min, log_max = np.log10([min(positive), max(positive)])
        margin = np.clip((log_max - log_min) * 0.045, 0.12, 0.55)
        log_limits = (10 ** (log_min - margin), 10 ** (log_max + margin))
    else:
        log_limits = (0.1, 10)

    for ax in log_axes:
        ax.set_ylim(*log_limits)
        ax.set_yscale("log", base=10)
        if positive:
            ax.yaxis.set_major_locator(LogLocator(base=10, numticks=5))
        else:
            ax.set_yticks([])
            ax.text(
                0.5,
                0.5,
                "全標本で相対差 0",
                transform=ax.transAxes,
                ha="center",
                va="center",
                color=MUTED,
                fontsize=10,
            )
        ax.yaxis.set_minor_locator(NullLocator())

    for ax in linear_axes:
        limit = 1.10 * linear_peak if linear_peak else 1
        ax.set_ylim(-limit, limit)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4, symmetric=True))
        # Keep the exponent on the actual ticks where the signed residuals are sub-ULP.
        if 0 < linear_peak < 0.001:
            ax.yaxis.set_major_formatter(lambda y, _: f"{y:g}".replace("-", "−"))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.axhline(0, color="#97a7ae", lw=0.85, zorder=1)


def curves_panel(axes, q, ps, *, zoom=False):
    top, abs_sashimi, signed_sashimi, abs_colossus, signed_colossus = axes
    x = np.array([float(Q(p["coordinate"])) for p in ps])
    series = active_series(ps)
    for base in reversed(BASES):
        verified = [
            (
                float(Q(p["coordinate"])),
                *sorted(float(end) * q.factor for end in box),
            )
            for p in ps
            if (box := plotting_enclosure(p, base)) is not None
        ]
        if verified:
            x_verified, lo, hi = zip(*verified, strict=True)
            interval_bars(
                top,
                np.array(x_verified),
                np.array(lo),
                np.array(hi),
                base,
                f"interval-value-{base}",
            )

    styles = {
        key: dict(
            **line_style(key),
            gid=key,
            markevery=(SERIES_KEYS.index(key), max(len(SERIES_KEYS), len(ps) // 8)),
        )
        for key in series
    }
    for key, style in styles.items():
        top.plot(x, [float(point_value(p, key)) * q.factor for p in ps], **style)
    positive = []
    boxes_by_base = {}
    panel_pairs = (
        ("legacy", abs_sashimi, signed_sashimi),
        ("matched", abs_colossus, signed_colossus),
    )
    for base, log_ax, lin_ax in panel_pairs:
        zeros = {}
        for key, style in styles.items():
            signed = np.array([ratio_for_plot(relative_difference(p, key, base)) for p in ps])
            is_zero = signed == 0
            magnitude = np.abs(signed)
            positive.extend(magnitude[~is_zero])
            if is_zero.any():
                zeros[key] = x[is_zero]
                magnitude[is_zero] = np.nan
            log_ax.plot(x, magnitude, **style)
            lin_ax.plot(x, signed, **style)
        boxes = {key: residual_intervals(ps, key, base) for key in BASES}
        boxes_by_base[base] = boxes
        exact_zero_strip(log_ax, zeros, boxes)

    nominal_peak = max(positive, default=0.0)
    for boxes in boxes_by_base.values():
        for box in boxes.values():
            for ends in (box.magnitude_lo, box.magnitude_hi):
                positive.extend(ends[ends > 0])
    linear_peak = nominal_peak if zoom else max(positive, default=0.0)
    set_residual_scales(
        (abs_sashimi, abs_colossus), (signed_sashimi, signed_colossus), positive, linear_peak
    )
    for base, log_ax, lin_ax in panel_pairs:
        for key in reversed(BASES):
            box = boxes_by_base[base][key]
            interval_bars(lin_ax, box.x, box.lo, box.hi, key, f"interval-signed-{key}")
            interval_bars(
                log_ax,
                box.x,
                box.magnitude_lo,
                box.magnitude_hi,
                key,
                f"interval-abs-{key}",
                log_floor=log_ax.get_ylim()[0],
            )

    top.set_title(q.title, loc="left", fontsize=14, weight="bold", pad=8)
    top.set_ylabel(q.unit)
    if q.log_y:
        top.set_yscale("log")
        top.yaxis.set_major_locator(LogLocator(numticks=5))
    else:
        top.yaxis.set_major_locator(MaxNLocator(nbins=4))
    for ax in axes:
        ax.grid(axis="y", color=GRID, lw=0.65)
        ax.set_axisbelow(True)
        ax.margins(x=0.025)
    for ax in axes[:-1]:
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    if q.group == "background":
        signed_colossus.set_xticks(range(8))
        signed_colossus.set_xlabel("赤方偏移 $z$")
    else:
        signed_colossus.set_xscale("log")
        signed_colossus.xaxis.set_major_locator(LogLocator(base=10, numticks=8))
        signed_colossus.xaxis.set_minor_locator(NullLocator())
        signed_colossus.set_xlabel(r"物理質量 $M$  [$M_\odot$]")


def curves(rows, quantities, title, subtitle, foot, *, zoom=False):
    groups = [[r for r in rows if r["quantity"] == q.key] for q in quantities]
    fig = plt.figure(figsize=(14.4, 10.8))
    fig.text(
        0.029,
        0.977,
        title + ("　｜　線形段を拡大" if zoom else ""),
        fontsize=22,
        weight="bold",
        va="top",
    )
    fig.text(0.031, 0.930, subtitle, fontsize=11, color=MUTED, va="top")
    available = {key for group in groups for key in active_series(group)}
    handles = [Line2D([], [], **line_style(k)) for k in SERIES_KEYS if k in available]
    matched_sources = {plot_enclosure_prefix(row, "matched") for group in groups for row in group}
    methods = []
    for method, label in (
        ("matched_tight", "固定セル改良"),
        ("matched_fine", "7,680 分割"),
        ("matched", "カーネル検査"),
    ):
        if method in matched_sources:
            methods.append(label)
    matched_label = "Colossus · " + "／".join(methods)
    for key, model in zip(BASES, ("SASHIMI · カーネル検査", matched_label), strict=True):
        handles.append(
            Line2D(
                [],
                [],
                ls="none",
                marker="|",
                ms=10,
                mew=2.4,
                color=SERIES_STYLE[key][1],
                label=f"Lean 区間 · {model}",
            )
        )
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.026, 0.907),
        frameon=False,
        ncol=4,
        handlelength=3.6,
        columnspacing=2.8,
        fontsize=10,
    )
    # The gaps after the logarithmic axes include the dedicated exact-zero strips.
    grid = fig.add_gridspec(
        9,
        len(groups),
        height_ratios=[1.30, 0.23, 1, 0.56, 1, 0.30, 1, 0.56, 1],
        left=0.169,
        right=0.981,
        top=0.812,
        bottom=0.13,
        wspace=0.26,
        hspace=0,
    )
    for col, (q, ps) in enumerate(zip(quantities, groups, strict=True)):
        top = fig.add_subplot(grid[0, col])
        axes = [top, *(fig.add_subplot(grid[row, col], sharex=top) for row in (2, 4, 6, 8))]
        curves_panel(axes, q, ps, zoom=zoom)
    for row, base, metric in (
        (2, "legacy", "相対差の絶対値 · 対数"),
        (4, "legacy", "符号付き相対差 · 線形"),
        (6, "matched", "相対差の絶対値 · 対数"),
        (8, "matched", "符号付き相対差 · 線形"),
    ):
        box = grid[row, 0].get_position(fig)
        y = (box.y0 + box.y1) / 2
        ref_name = "SASHIMI 準拠 Lean" if base == "legacy" else "Colossus 準拠 Lean"
        fig.text(
            0.023,
            y + 0.007,
            ref_name + " 基準",
            fontsize=10,
            weight="bold",
            color=SERIES_STYLE[base][1],
            va="bottom",
        )
        fig.text(0.023, y - 0.002, metric, fontsize=10, color=INK, va="top")
    fig.text(
        0.031,
        0.055,
        "縦棒は各標本の区間。積分和は証明済み評価器で計算し、上下限への伝播をカーネルで検査。各点の計算法を比較 CSV に記録。"
        if any(s in matched_sources for s in ("matched_fine", "matched_tight"))
        else "縦棒は各標本でカーネル検査した Lean の区間。",
        fontsize=10,
        color=INK,
    )
    fig.text(
        0.031,
        0.033,
        "相対差は r = (x − b) / b。基準値 b を固定して区間を換算。「一致」は x = b、「0 を含む」は区間に b が含まれる標本。"
        + ("　線形段の区間全幅は通常版に掲載。" if zoom else ""),
        fontsize=9.5,
        color=MUTED,
    )
    fig.text(0.031, 0.013, foot, fontsize=9, color=MUTED)
    return fig


def write_figures(out, filename, figures):
    with PdfPages(out / filename) as pdf:
        for name, fig in figures:
            fig.savefig(out / (name + ".png"), dpi=160)
            pdf.savefig(fig)
            plt.close(fig)


def _render(data, baseline, out):
    rows, bounds, mass_boxes = _comparison_data(data, baseline)
    out.mkdir(parents=True, exist_ok=True)
    save_tables(out, rows)
    style()
    pages = [
        (
            "matched-mass",
            MASS[:2],
            "密度揺らぎの標準偏差と分散",
            r"$z=0$ の $\sigma(M)$ と $S(M)=\sigma^2(M)$　　物理質量 $10^{-6}$–$10^{12}\,M_\odot$",
            "SASHIMI は Ludlow フィット、Colossus は Eisenstein–Hu と球形トップハット。同一帯域の系列は Colossus の関数を再積分し、基準半径で再規格化。",
        ),
        (
            "matched-derivative-kernel",
            MASS[2:],
            "分散の質量微分と EPS 条件付き密度",
            r"親 $10^{12}\,M_\odot$・$z=0$、子 $z=1$ の条件付き密度",
            "左上は負の質量微分の大きさ。相対差は元の符号で計算。右の条件付き密度は各処方の分散と成長を同じ Yang 式に入力。",
        ),
        (
            "matched-background",
            BACKGROUND[:2],
            "宇宙膨張と線形成長",
            r"平坦 ΛCDM　　$h=0.674$　　$\Omega_m=0.315$　　$D(0)=1$",
            "膨張率は両 Lean モデルで同じ式。Astropy の成長因子はその背景膨張率から積分して評価。",
        ),
        (
            "matched-growth-barrier",
            BACKGROUND[2:],
            "成長因子の赤方偏移微分と崩壊障壁",
            r"平坦 ΛCDM　　$h=0.674$　　$\Omega_m=0.315$　　$\delta_c(z)=1.686/D(z)$",
            "両処方で同じ崩壊閾値を採用。Astropy の成長量はその背景膨張率から評価。",
        ),
        (
            "matched-growth-host",
            HOST,
            "ホストハローの質量履歴と質量変化率",
            r"$M_{200}(0)=10^{12}\,M_\odot$　　各処方の分散と成長を同じ SASHIMI-C の質量履歴の式に入力",
            "右の質量変化率は SASHIMI-C のビリアル質量フィットを含む。",
        ),
    ]
    write_figures(
        out,
        "matched-comparison.pdf",
        ((name, curves(rows, qs, title, sub, foot)) for name, qs, title, sub, foot in pages),
    )
    with PdfPages(out / "matched-linear-zoom.pdf") as pdf:
        for _, qs, title, subtitle, foot in pages:
            if qs[0] in MASS or qs[0] in HOST:
                fig = curves(rows, qs, title, subtitle, foot, zoom=True)
                pdf.savefig(fig)
                plt.close(fig)
    diagnostics = [
        ("matched-enclosures", save_and_plot(mass_boxes, bounds, out)),
        ("matched-precision", error_budget(mass_boxes)),
        (
            "matched-transfer-growth",
            point_proofs(
                load(data, "transfer"),
                load(data, "growth"),
                load(data, "colossus"),
                load(baseline, "external"),
                load(data, "transfer-samples"),
                load(baseline, "background"),
            ),
        ),
    ]
    write_figures(out, "matched-diagnostics.pdf", diagnostics)


def render(data, baseline, out):
    out.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=out) as tmp:
        stage = Path(tmp)
        _render(data, baseline, stage)
        for path in stage.iterdir():
            path.replace(out / path.name)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    for arg in ("data", "baseline", "output"):
        cli.add_argument(arg, type=Path)
    args = cli.parse_args()
    render(args.data, args.baseline, args.output)
