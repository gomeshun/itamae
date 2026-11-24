"""Plot the matched spectrum against Colossus on the published SASHIMI-C EPS grid."""

import argparse
import csv
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, MaxNLocator, NullLocator
import numpy as np

if __package__:
    from .plot_eh_bounds import error_budget, intervals, save_and_plot, point_proofs
    from .plot_eps import BACKGROUND, MASS, HOST, INK, MUTED, GRID, style
else:
    from plot_eh_bounds import error_budget, intervals, save_and_plot, point_proofs
    from plot_eps import BACKGROUND, MASS, HOST, INK, MUTED, GRID, style


BASES = ("legacy", "matched")
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
    for suffix in ("", "-fine"):
        name = "kernel-bounds" + suffix
        if name in bounds:
            for k, expected in bounds[name]["inputs_sha256"].items():
                stems = {"full": "full" + suffix, "growth": "growth", "paired": "paired-children"}
                file_name = stems[k] + ".json"
                if sha256((data / file_name).read_bytes()).hexdigest() != expected:
                    raise ValueError(f"{name} inputs differ from {file_name}")
    mass_boxes = intervals(rows, bounds)
    mass_lookup = {(r["quantity"], r["mass_msun"]): r for r in mass_boxes}
    host_report = load(data, "host-bounds")
    host_rows = keyed(host_report["rows"], "z")
    host_fine = keyed(load(data, "host-bounds-fine")["rows"], "z")
    for r in rows:
        key = r["quantity"]
        r["matched_fine_lo_exact"], r["matched_fine_hi_exact"] = "", ""
        result = mass_lookup.get((key, r["coordinate"]))
        if result is not None:
            r["matched_lo_exact"], r["matched_hi_exact"] = result["lo"], result["hi"]
            r["matched_fine_lo_exact"], r["matched_fine_hi_exact"] = (
                result["fine_lo"],
                result["fine_hi"],
            )
            source, checked = result["inputs"], result["end_to_end_kernel_replay"]
        elif key in {q.key for q in HOST}:
            host = host_rows[Q(r["coordinate"])]
            box = host[key]
            tight = host_fine[Q(r["coordinate"])][key]
            if tight:
                r["matched_fine_lo_exact"], r["matched_fine_hi_exact"] = tight["lo"], tight["hi"]
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
            ("kernel_checked" if checked else "proved_evaluator_runtime")
            if r["matched_lo_exact"]
            else "unresolved"
        )
        r["matched_enclosure_source"] = source
        r["matched_enclosure_status"] = r["enclosure_status"]
        # Keep the nominal numerical line even when an integral was unresolved.
        # Every published bound is checked against the exact reference used for this chart.
        for lo, hi in (
            ("matched_lo_exact", "matched_hi_exact"),
            ("matched_fine_lo_exact", "matched_fine_hi_exact"),
        ):
            if bool(r[lo]) != bool(r[hi]):
                raise ValueError(f"Incomplete bounds for {key} at {r['coordinate']}")
            if r[lo] and not Q(r[lo]) <= point_value(r, "matched") <= Q(r[hi]):
                raise ValueError(
                    f"Matched reference outside {lo}/{hi} for {key} at {r['coordinate']}"
                )
    return rows, bounds, mass_boxes


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
LINEAR_CORE_PCT = 1e-6


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


def residual_axes(axes, peak):
    """Keep both reference models on one signed percentage scale per observable."""
    for ax in axes:
        if peak <= 10 * LINEAR_CORE_PCT:
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, symmetric=True))
            ax.ticklabel_format(
                axis="y", style="sci", scilimits=(-3, 3), useOffset=False, useMathText=True
            )
            limit = 1.12 * peak if peak else LINEAR_CORE_PCT
            ax.set_ylim(-limit, limit)
        else:
            ax.set_yscale("symlog", base=10, linthresh=LINEAR_CORE_PCT, linscale=0.65)
            outer_exp = int(np.floor(np.log10(peak)))
            inner_exp = (outer_exp + int(np.log10(LINEAR_CORE_PCT))) // 2
            outer, inner = 10.0**outer_exp, 10.0**inner_exp
            a, b = ["1" if e == 0 else rf"10^{{{e}}}" for e in (outer_exp, inner_exp)]
            ax.set_yticks(
                [-outer, -inner, 0, inner, outer],
                ["$-" + a + "$", "$-" + b + "$", "0", "$" + b + "$", "$" + a + "$"],
            )
            ax.set_ylim(-2 * peak, 2 * peak)
        ax.yaxis.set_minor_locator(NullLocator())


def curves_panel(ax_tuple, q, ps):
    top, mid, low = ax_tuple
    x = [float(Q(p["coordinate"])) for p in ps]
    for base in BASES:
        verified = [
            (
                float(Q(p["coordinate"])),
                *sorted(float(Q(p[f"{base}_{side}_exact"])) * q.factor for side in ("lo", "hi")),
            )
            for p in ps
            if p[base + "_enclosure_status"] == "kernel_checked"
        ]
        if verified:
            x_verified, lo, hi = zip(*verified, strict=True)
            top.vlines(
                x_verified, lo, hi, color=SERIES_STYLE[base][1], lw=1.4, alpha=0.27, zorder=1
            )
    peak = 0
    for key in active_series(ps):
        style = line_style(key)
        style["gid"] = key
        # Stagger the markers, never the data, so coincident series remain identifiable.
        style["markevery"] = (SERIES_KEYS.index(key), max(len(SERIES_KEYS), len(ps) // 8))
        top.plot(x, [float(point_value(p, key)) * q.factor for p in ps], **style)
        for ax, base in zip((mid, low), BASES, strict=True):
            pct = [float(100 * relative_difference(p, key, base)) for p in ps]
            ax.plot(x, pct, **style)
            peak = max(peak, max(map(abs, pct)))
    residual_axes((mid, low), peak)
    top.set_title(q.title, loc="left", fontsize=13.5, weight="bold", pad=8)
    top.set_ylabel(q.unit)
    if q.log_y:
        top.set_yscale("log")
        top.yaxis.set_major_locator(LogLocator(numticks=5))
    else:
        top.yaxis.set_major_locator(MaxNLocator(nbins=4))
    mid.set_ylabel("Lean · SASHIMI 基準\n相対差 [%]", color=SERIES_STYLE["legacy"][1], fontsize=8.5)
    low.set_ylabel(
        "Lean · Colossus 基準\n相対差 [%]", color=SERIES_STYLE["matched"][1], fontsize=8.5
    )
    for ax in ax_tuple:
        ax.grid(axis="y", color=GRID, lw=0.65)
        ax.set_axisbelow(True)
        ax.margins(x=0.025)
    for ax in (top, mid):
        ax.tick_params(axis="x", which="both", labelbottom=False)
    if q.group == "background":
        low.set_xticks(range(8))
        low.set_xlabel("赤方偏移 $z$")
    else:
        low.set_xscale("log")
        low.xaxis.set_major_locator(LogLocator(base=10, numticks=8))
        low.xaxis.set_minor_locator(NullLocator())
        low.set_xlabel(r"質量 $M$  [$M_\odot$]")


def curves(rows, quantities, title, subtitle, foot):
    groups = [[r for r in rows if r["quantity"] == q.key] for q in quantities]
    nrows = (len(quantities) + 1) // 2
    fig = plt.figure(figsize=(13.0, 5.85 * nrows + 1.8))
    fig.text(0.096, 0.967, title, fontsize=21, weight="bold", va="top")
    fig.text(0.096, 0.926 if nrows == 2 else 0.900, subtitle, fontsize=10.5, color=MUTED, va="top")
    available = {key for group in groups for key in active_series(group)}
    keys = [key for key in SERIES_KEYS if key in available]
    fig.legend(
        handles=[Line2D([], [], **line_style(k)) for k in keys],
        loc="upper left",
        bbox_to_anchor=(0.091, 0.903 if nrows == 2 else 0.864),
        frameon=False,
        ncol=3,
        handlelength=3.5,
        columnspacing=2.4,
        fontsize=10,
    )
    cells = fig.add_gridspec(
        nrows,
        2,
        left=0.10,
        right=0.981,
        top=0.813 if nrows == 2 else 0.718,
        bottom=0.121 if nrows == 2 else 0.205,
        wspace=0.29,
        hspace=0.36,
    )
    for cell, q, group in zip(cells, quantities, groups, strict=True):
        sub = cell.subgridspec(3, 1, height_ratios=[1.33, 1, 1], hspace=0.18)
        top = fig.add_subplot(sub[0])
        mid = fig.add_subplot(sub[1], sharex=top)
        low = fig.add_subplot(sub[2], sharex=top)
        curves_panel((top, mid, low), q, group)
    fig.text(
        0.096,
        0.066 if nrows == 2 else 0.095,
        "上段の縦線は両 Lean モデルの保証区間。中・下段は 100 × (各系列 / 指定モデルの Lean 代表値 − 1)。",
        fontsize=10,
        color=INK,
    )
    h_note = "膨張率のみ線形軸。" if any(q.key == "H" for q in quantities) else ""
    fig.text(
        0.096,
        0.044 if nrows == 2 else 0.064,
        h_note
        + rf"相対差は符号付き対数軸（ゼロ付近の $\pm 10^{{{int(np.log10(LINEAR_CORE_PCT))}}}$ % は線形）。同じ量の中・下段は共通目盛。",
        fontsize=9.5,
        color=MUTED,
    )
    fig.text(0.096, 0.022 if nrows == 2 else 0.032, foot, fontsize=9.5, color=MUTED)
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
    comparison = [
        (
            "matched-mass",
            curves(
                rows,
                MASS,
                "質量分散と条件付き密度",
                r"$z=0$ の質量分散　　親 $10^{12}M_\odot$・$z=0$、子 $z=1$ の EPS 条件付き密度",
                "同一帯域の参照値は Colossus の伝達関数を再積分し、基準半径でも再規格化。条件付き密度は各方式の分散と成長を同じ Yang 式に入力。",
            ),
        ),
        (
            "matched-background",
            curves(
                rows,
                BACKGROUND,
                "宇宙膨張・線形成長・崩壊障壁",
                r"平坦 ΛCDM　　$h=0.674$　　$\Omega_m=0.315$　　$D(0)=1$",
                "膨張率は両 Lean モデルで同じ式。Astropy の成長量にはその膨張率から求めた成長積分を使用。",
            ),
        ),
        (
            "matched-growth-host",
            curves(
                rows,
                HOST,
                "ホストハローの質量履歴と質量変化率",
                r"$M_{200}(0)=10^{12}M_\odot$　　各方式の分散と成長を同じ SASHIMI-C の質量処方に入力",
                "質量変化率は SASHIMI-C が用いるビリアル質量フィットを含む。重なる線は異なる位置のマーカーで表示。",
            ),
        ),
    ]
    write_figures(out, "matched-comparison.pdf", comparison)
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
