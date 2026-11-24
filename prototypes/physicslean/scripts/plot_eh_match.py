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

MATCH, OLD, SPLINE, DIRECT, SAME = "#163d59", "#8b6976", "#c26311", "#00817e", "#4374ae"


def load(p, name):
    return json.loads((p / (name + ".json")).read_text())


def keyed(rows, field):
    rs = {Q(str(r[field])): r for r in rows}
    if len(rs) != len(rows):
        raise ValueError(f"Duplicate coordinate in {field}")
    return rs


def midpoint(box):
    return float((Q(box["lo"]) + Q(box["hi"])) / 2)


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
    old = {
        group: load(baseline, name)["rows"]
        for group, name in [
            ("background", "background"),
            ("variance", "variance"),
            ("conditional", "kernel"),
        ]
    }
    records = []
    for q in BACKGROUND + MASS + HOST:
        t = q.group
        is_bg = t == "background"
        fld = "z" if is_bg else "mass_msun"
        coords = zs if is_bg else (mass - {10**12} if t == "conditional" else mass)
        orig = keyed(old[t], fld)
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
        if any(set(x) != coords for x in (orig, comp, new)):
            raise ValueError(f"Grid mismatch for {q.key}")
        host_refs = {m: keyed(down[m][t], fld) for m in ("spline", "direct")} if q in HOST else {}
        for x in sorted(coords):
            old_value = midpoint(orig[x][q.key])
            val = new[x][q.key]
            matched = midpoint(val) if isinstance(val, dict) else val
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
                legacy=old_value,
                sashimi=sashimi[x][q.key],
                astropy_growth=astro[x].get(q.key) if is_bg else None,
                matched=matched,
                matched_lo_exact=verified_lo,
                matched_hi_exact=verified_hi,
                colossus_spline=rs,
                colossus_direct=rd,
                fixed_band=rf,
            )
            for name, ref, value in [
                ("before", rs, old_value),
                ("after", rs, matched),
                ("direct", rd, matched),
                ("fixed", rf, matched),
            ]:
                row["rel_" + name] = None if ref is None else float((Q(ref) - Q(value)) / Q(value))
            records.append(row)
    return records


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
            return Q(r["matched_lo_exact"]) <= Q(r[field]) <= Q(r["matched_hi_exact"])

        summary = {
            "points": len(subset),
            "enclosures": {
                "returned": len(enclosed),
                "end_to_end_kernel_checked": len(checked),
                "plotted_value_inside": sum(inside(r, "matched") for r in enclosed),
                "colossus_spline_inside": sum(inside(r, "colossus_spline") for r in enclosed),
            },
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


def header(fig, title, subtitle):
    fig.text(0.09, 0.975, title, fontsize=20, weight="bold", va="top")
    fig.text(0.09, 0.925, subtitle, fontsize=10.2, color=MUTED, va="top")


def curves_panel(ax_tuple, q, ps):
    top, before, after = ax_tuple
    xx = np.array([float(Q(p["coordinate"])) for p in ps])

    def series(k):
        return np.array([p[k] for p in ps], dtype=float)

    top.set_title(q.title, loc="left", fontsize=12.7, weight="bold", pad=5)
    if q in HOST:
        for x, p in zip(xx, ps, strict=True):
            if p["matched_lo_exact"]:
                lo, hi = sorted(
                    float(Q(p[k])) * q.factor for k in ("matched_lo_exact", "matched_hi_exact")
                )
                top.vlines(x, lo, hi, color="#74aac4", lw=1.25, zorder=0)
            if p["matched_fine_lo_exact"]:
                lo, hi = sorted(
                    float(Q(p[k])) * q.factor
                    for k in ("matched_fine_lo_exact", "matched_fine_hi_exact")
                )
                top.vlines(x, lo, hi, color="#176a61", lw=1.8, zorder=1)

    top.plot(xx, series("legacy") * q.factor, color=OLD, lw=1.55, ls=(0, (5, 2)))
    top.plot(
        xx,
        series("sashimi") * q.factor,
        color=OLD,
        marker="x",
        ms=4,
        ls="none",
        markevery=(1, max(1, len(xx) // 9)),
    )
    top.plot(xx, series("matched") * q.factor, color=MATCH, lw=2)
    top.plot(
        xx,
        series("colossus_spline") * q.factor,
        color=SPLINE,
        marker="o",
        mfc="none",
        ms=4,
        mew=1,
        ls="none",
        markevery=max(1, len(xx) // 14),
    )
    if all(p["astropy_growth"] is not None for p in ps):
        top.plot(
            xx,
            series("astropy_growth") * q.factor,
            color="#315c39",
            marker="1",
            ls="none",
            ms=7,
            markevery=(1, 4),
        )
    top.set_ylabel(q.unit)
    if q.log_y:
        top.set_yscale("log")
        top.yaxis.set_major_locator(LogLocator(numticks=5))
    else:
        top.yaxis.set_major_locator(MaxNLocator(nbins=4))
    before.axhline(0, color="#8e9eaa", lw=0.7)
    before.plot(xx, 100 * series("rel_before"), color=OLD, lw=1.55, ls=(0, (4, 2)))
    before.plot(xx, 100 * series("rel_after"), color=SPLINE, lw=1.55)
    before.set_ylabel("既定値との差 [%]", fontsize=8.5)
    after.axhline(0, color="#8e9eaa", lw=0.7)
    after.plot(xx, 100 * series("rel_after"), color=SPLINE, lw=1.55)
    after.plot(xx, 100 * series("rel_direct"), color=DIRECT, lw=1.2, ls=(0, (3, 2)))
    if all(p["rel_fixed"] is not None for p in ps):
        after.plot(
            xx,
            100 * series("rel_fixed"),
            color=SAME,
            ls="none",
            marker="+",
            ms=3.4,
            markevery=max(1, len(xx) // 22),
        )
    after.set_ylabel("整合後 [%]", fontsize=8.5)
    for a in (before, after):
        a.yaxis.set_major_locator(MaxNLocator(nbins=3))
    for a in ax_tuple:
        a.grid(axis="y", color=GRID, lw=0.65)
        a.set_axisbelow(True)
        a.margins(x=0.022)
    for a in (top, before):
        a.tick_params(axis="x", which="both", labelbottom=False)
    if q.group == "background":
        after.set_xticks(range(8))
        after.set_xlabel("赤方偏移 $z$")
    else:
        after.set_xscale("log")
        after.xaxis.set_major_locator(LogLocator(base=10, numticks=8))
        after.xaxis.set_minor_locator(NullLocator())
        after.set_xlabel(r"質量 $M$  [$M_\odot$]")


def curves(rows, quantities, title, subtitle, foot, fixed):
    fig = plt.figure(figsize=(12.0, 10.8))
    header(fig, title, subtitle)
    labels = [
        (OLD, (0, (5, 2)), "前回の SASHIMI-C 処方"),
        (MATCH, "-", "Lean 整合モデル"),
        (SPLINE, "-", "Colossus 既定（補間）"),
        (DIRECT, (0, (3, 2)), "Colossus 直接計算"),
    ]
    handles = [Line2D([], [], color=c, lw=1.8, ls=ls, label=label) for c, ls, label in labels]
    handles.append(Line2D([], [], color=OLD, ls="none", marker="x", label="SASHIMI-C 実装"))
    if not fixed:
        handles.append(
            Line2D([], [], color="#315c39", ls="none", marker="1", ms=7, label="Astropy + 成長積分")
        )
    if fixed:
        handles.append(Line2D([], [], color=SAME, ls="none", marker="+", label="同一帯域・再積分"))
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.085, 0.905),
        frameon=False,
        ncol=3,
        handlelength=2.6,
        columnspacing=1.25,
        fontsize=9.4,
        labelspacing=0.42,
    )
    cells = fig.add_gridspec(
        2, 2, left=0.09, right=0.982, top=0.825, bottom=0.134, hspace=0.45, wspace=0.255
    )
    for cell, q in zip(cells, quantities, strict=True):
        trio = cell.subgridspec(3, 1, height_ratios=[1.65, 1, 1], hspace=0.10)
        a = fig.add_subplot(trio[0])
        b = fig.add_subplot(trio[1], sharex=a)
        c = fig.add_subplot(trio[2], sharex=a)
        curves_panel((a, b, c), q, [p for p in rows if p["quantity"] == q.key])
    fig.text(
        0.09,
        0.055,
        "中段は Colossus 既定値 / 各モデル − 1。下段は基準を Lean 整合モデルにそろえて拡大。",
        fontsize=9.7,
        color=INK,
    )
    fig.text(0.09, 0.027, foot, fontsize=9, color=MUTED)
    return fig


def _render(data, baseline, out):
    rows = join(data, baseline)
    out.mkdir(parents=True, exist_ok=True)
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
    save_tables(out, rows)
    style()
    figs = [
        (
            "matched-mass",
            curves(
                rows,
                MASS,
                "質量分散と EPS 密度のモデルをそろえる",
                r"バリオンを含む Eisenstein–Hu (1998)　 球形トップハット　 $h=0.674$　 $\Omega_m=0.315$　 $\Omega_b=0.049$　 $n_s=0.965$　 $\sigma_8=0.811$",
                "分散は z=0。条件付き密度は親 $10^{12}$ 太陽質量・z=0、子 z=1。下段の直接計算と同一帯域積分の条件は比較表に記載。",
                True,
            ),
        ),
        (
            "matched-growth-host",
            curves(
                rows,
                [BACKGROUND[1], BACKGROUND[2], *HOST],
                "線形成長とホスト質量に伝わる差",
                r"平坦な物質＋Λの成長解を $D(0)=1$ で規格化　 ホスト $M_{200}(0)=10^{12}\,M_\odot$",
                "ホストの縦線は960分割（青）・7,680分割（緑）の区間。比較値は各実装の分散・成長率を同じ質量処方に入力。",
                False,
            ),
        ),
    ]
    figs.append(("matched-enclosures", save_and_plot(mass_boxes, bounds, out)))
    figs.append(("matched-precision", error_budget(mass_boxes)))
    figs.append(
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
        )
    )
    with PdfPages(out / "matched-comparison.pdf") as pdf:
        for name, fig in figs:
            fig.savefig(out / (name + ".png"), dpi=160)
            pdf.savefig(fig)
            plt.close(fig)


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
