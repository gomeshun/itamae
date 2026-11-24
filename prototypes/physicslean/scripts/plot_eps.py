"""Plot kernel-checked SASHIMI-C EPS inputs against the three numerical references.

Run with the directory containing grid, background, variance, kernel, sashimi and
external JSON files. Comparisons are performed on exact rationals before rendering.
"""

import argparse
import csv
from dataclasses import dataclass
from fractions import Fraction
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, MaxNLocator, NullLocator
import numpy as np


LEAN, NATIVE, COLOSSUS, ASTROPY = "#183b59", "#92517c", "#bd6209", "#087d7d"
BAND, INK, MUTED, GRID = "#53a4cf", "#172f41", "#49616f", "#d8e1e5"
MODELS = {
    "sashimi": dict(color=NATIVE, marker="o", mfc="none", ms=3.4, mew=0.85, ls="none"),
    "colossus": dict(color=COLOSSUS, ls=(0, (5, 2)), lw=1.65),
    "astropy": dict(color=ASTROPY, ls=(0, (2, 1.5)), lw=1.5, marker="x", ms=4, markevery=(1, 4)),
}
LABELS = {"sashimi": "SASHIMI-C", "colossus": "Colossus", "astropy": "Astropy + 成長積分"}


@dataclass(frozen=True)
class Quantity:
    key: str
    group: str
    title: str
    unit: str
    short: str
    factor: float = 1.0
    log_y: bool = False


BACKGROUND = [
    Quantity("H", "background", "膨張率", r"$H(z)$  [km s$^{-1}$ Mpc$^{-1}$]", r"$H$"),
    Quantity("D", "background", "線形成長因子", r"$D(z)$", r"$D$"),
    Quantity(
        "dD_dz",
        "background",
        "成長因子の赤方偏移微分",
        r"$\mathrm{d}D/\mathrm{d}z$",
        r"$\mathrm{d}D/\mathrm{d}z$",
    ),
    Quantity(
        "delta_c", "background", "球対称崩壊の障壁", r"$\delta_c(z)=1.686/D(z)$", r"$\delta_c$"
    ),
]
MASS = [
    Quantity("sigma", "variance", "密度揺らぎの標準偏差", r"$\sigma(M,0)$", r"$\sigma$"),
    Quantity("S", "variance", "密度揺らぎの分散", r"$S(M)=\sigma^2(M,0)$", r"$S$"),
    Quantity(
        "dS_dM",
        "variance",
        "分散の質量微分の大きさ",
        r"$-\mathrm{d}S/\mathrm{d}M$  [$M_\odot^{-1}$]",
        r"$\mathrm{d}S/\mathrm{d}M$",
        -1,
        True,
    ),
    Quantity(
        "F", "conditional", "条件付き初通過分布", r"$F$  [単位 $S$ あたり]", r"$F$", log_y=True
    ),
]
HOST = [
    Quantity(
        "M200",
        "background",
        "主ハローの質量履歴",
        r"$M_{200}(z)$  [$10^{11}\, M_\odot$]",
        r"$M_{200}$",
        1e-11,
    ),
    Quantity(
        "dM_dz",
        "background",
        "SASHIMI-C の質量成長率",
        r"$\mathtt{dMdz}$  [$10^{11}\, M_\odot$ / 赤方偏移]",
        r"$\mathtt{dMdz}$",
        1e-11,
    ),
]
QUANTITIES = BACKGROUND + MASS + HOST


@dataclass
class Point:
    x: Fraction
    lo: Fraction
    hi: Fraction
    refs: dict

    @property
    def mid(self):
        return (self.lo + self.hi) / 2

    @property
    def radius(self):
        return (self.hi - self.lo) / (2 * abs(self.mid))

    def residual(self, model):
        # Fraction(float) preserves the exact binary64 value, including its rounding.
        return (Fraction(self.refs[model]) - self.mid) / self.mid


def indexed(rows, field):
    points = {}
    for row in rows:
        key = Fraction(row[field])
        if key in points:
            raise ValueError(f"Repeated {field}={key}")
        points[key] = row
    return points


def load_data(root):
    def read(name):
        return json.loads((root / (name + ".json")).read_text())

    grid, bg, var, cond = (read(name) for name in ("grid", "background", "variance", "kernel"))
    if (
        bg["kernel_replayed"] is not True
        or var["kernel_replay"] != "passed"
        or cond["kernel_replayed"] is not True
    ):
        raise ValueError(
            "The grid requires a successful concrete Lean replay for all three datasets"
        )
    if any(Fraction(r["z"]) != 0 for r in var["rows"]):
        raise ValueError("Expected a z=0 variance grid")
    inputs = {"background": bg, "variance": var, "conditional": cond}
    refs = {
        "sashimi": read("sashimi"),
        **{k: v for k, v in read("external").items() if k in MODELS},
    }
    zs = sorted(map(Fraction, grid["redshifts"]))
    ms = sorted(map(Fraction, grid["masses_msun"]))
    expected = dict(background=zs, variance=ms, conditional=[m for m in ms if m < 10**12])
    result = {}
    for q in QUANTITIES:
        field = "z" if q.group == "background" else "mass_msun"
        cert = indexed(inputs[q.group]["rows"], field)
        if sorted(cert) != expected[q.group]:
            raise ValueError(f"Incomplete certificate grid for {q.key}")
        needed = ["sashimi"]
        if q in BACKGROUND:
            needed.extend(("colossus", "astropy"))
        elif q in MASS:
            needed.append("colossus")
        ref = {label: indexed(refs[label][q.group], field) for label in needed}
        if any(sorted(rows) != expected[q.group] for rows in ref.values()):
            raise ValueError(f"Incomplete reference grid for {q.key}")
        points = []
        for x in expected[q.group]:
            box = cert[x][q.key]
            lo, hi = (Fraction(box[end]) for end in ("lo", "hi"))
            values = {k: rows[x][q.key] for k, rows in ref.items()}
            if lo > hi or lo * hi <= 0 or not all(map(math.isfinite, values.values())):
                raise ValueError(f"Invalid interval or reference for {q.key} at {x}")
            points.append(Point(x, lo, hi, values))
        result[q.key] = points
    return result


def tables(out, data):
    with (out / "samples.csv").open("w", newline="") as stream:
        w = csv.writer(stream, lineterminator="\n")
        w.writerow(
            [
                "quantity",
                "z",
                "mass_Msun",
                "lo_exact",
                "hi_exact",
                "relative_half_width",
                "SASHIMI_C_binary64",
                "Colossus_binary64",
                "Astropy_based_binary64",
                "SASHIMI_C_in_interval",
            ]
        )
        for q in QUANTITIES:
            for p in data[q.key]:
                z = p.x if q.group == "background" else int(q.group == "conditional")
                mass = p.x if q.group != "background" else (10**12 if q in HOST else "")
                w.writerow(
                    [
                        q.key,
                        z,
                        mass,
                        p.lo,
                        p.hi,
                        float(p.radius),
                        *(p.refs.get(k, "") for k in MODELS),
                        abs(p.residual("sashimi")) <= p.radius,
                    ]
                )
    results = {}
    for q in QUANTITIES:
        ps = data[q.key]
        results[q.key] = {
            "n": len(ps),
            "max_relative_full_width": float(2 * max(p.radius for p in ps)),
            "sashimi_inside": sum(abs(p.residual("sashimi")) <= p.radius for p in ps),
            "max_abs_relative_difference": {
                m: float(max(abs(p.residual(m)) for p in ps)) for m in ps[0].refs
            },
        }
    (out / "summary.json").write_text(json.dumps(results, indent=2) + "\n")


def style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Noto Sans CJK JP", "DejaVu Sans"],
            "mathtext.fontset": "stixsans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#9eb0b9",
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.labelsize": 10,
            "axes.unicode_minus": True,
            "pdf.fonttype": 42,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def heading(fig, title, sub, models):
    fig.text(0.09, 0.974, title, va="top", fontsize=22, weight="bold")
    fig.text(0.09, 0.918, sub, va="top", fontsize=10.6, color=MUTED)
    handles = [
        Line2D([], [], color=LEAN, lw=2.2, label="ITAMAE 区間中心"),
        Line2D(
            [], [], color=BAND, marker="|", ms=11, mew=2, ls="none", label="Lean 保証区間（各点）"
        ),
    ]
    handles.extend(Line2D([], [], **MODELS[k], label=LABELS[k]) for k in models)
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.082, 0.872),
        frameon=False,
        ncol=len(handles),
        fontsize=9.5,
        columnspacing=1.5,
        handlelength=2.7,
        borderaxespad=0,
    )


def facet(fig, cell, q, ps):
    pair = cell.subgridspec(2, 1, height_ratios=(2.25, 1), hspace=0.055)
    top = fig.add_subplot(pair[0])
    delta = fig.add_subplot(pair[1], sharex=top)
    x = [float(p.x) for p in ps]
    models = list(ps[0].refs)
    top.set_title(q.title, loc="left", fontsize=12, weight="bold", pad=5)
    top.plot(x, [float(p.mid) * q.factor for p in ps], color=LEAN, lw=2.4)
    top.vlines(
        x,
        [float(p.lo) * q.factor for p in ps],
        [float(p.hi) * q.factor for p in ps],
        color=BAND,
        linewidth=2.1,
        zorder=3,
    )
    for model in models:
        top.plot(
            x,
            [p.refs[model] * q.factor for p in ps],
            **MODELS[model],
            zorder=5 if model == "sashimi" else 3,
        )
    top.set_ylabel(q.unit)
    if q.log_y:
        top.set_yscale("log")
        top.yaxis.set_major_locator(LogLocator(numticks=6))
    else:
        top.yaxis.set_major_locator(MaxNLocator(nbins=4))
    max_rel = max(max([p.radius] + [abs(p.residual(m)) for m in models]) for p in ps)
    if max_rel < Fraction("1e-8"):
        power = math.floor(math.log10(float(max_rel)))
        scale = 10.0 ** (-power)
        label = rf"相対差  [$10^{{{power}}}$]"
    else:
        scale, label = 100.0, "相対差  [%]"
    half = np.array([float(p.radius) * scale for p in ps])
    delta.axhline(0, color=LEAN, lw=0.9)
    delta.vlines(x, -half, half, color=BAND, lw=2.2, zorder=3)
    for model in models:
        delta.plot(
            x,
            [float(p.residual(model)) * scale for p in ps],
            **MODELS[model],
            zorder=5 if model == "sashimi" else 3,
        )
    delta.set_ylabel(label, fontsize=8.5)
    delta.yaxis.set_major_locator(MaxNLocator(nbins=3))
    top.tick_params(axis="x", which="both", labelbottom=False)
    if q.group != "background":
        top.set_xscale("log")
        delta.xaxis.set_major_locator(LogLocator(base=10, numticks=8))
        delta.xaxis.set_minor_locator(NullLocator())
        delta.set_xlabel(r"質量 $M$  [$M_\odot$]")
    else:
        delta.set_xticks(range(math.ceil(min(x)), math.floor(max(x)) + 1))
        delta.set_xlabel(r"赤方偏移 $z$")
    for ax in (top, delta):
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.margins(x=0.025)
        ax.set_axisbelow(True)


def curves(data, title, subtitle, quantities, models, note):
    compact = len(quantities) == 2
    fig = plt.figure(figsize=(11.9, 6.5 if compact else 10))
    heading(fig, title, subtitle, models)
    cells = fig.add_gridspec(
        1 if compact else 2,
        2,
        left=0.09,
        right=0.977,
        top=0.76 if compact else 0.79,
        bottom=0.19 if compact else 0.13,
        wspace=0.27,
        hspace=0.55,
    )
    for cell, q in zip(cells, quantities):
        facet(fig, cell, q, data[q.key])
    fig.text(0.09, 0.045, note, fontsize=9.3, color=MUTED, va="bottom")
    return fig


def accuracy(data):
    fig = plt.figure(figsize=(11.9, 7.3))
    fig.text(0.09, 0.968, "保証区間の幅と SASHIMI-C の丸め値", fontsize=21, va="top", weight="bold")
    fig.text(
        0.09,
        0.893,
        "各記号は一つの標本。区間半幅と実装値の差を、同じ区間中心の絶対値で割って比較。",
        fontsize=10.8,
        color=MUTED,
    )
    ax = fig.add_axes((0.145, 0.16, 0.665, 0.62))
    count = fig.add_axes((0.83, 0.16, 0.14, 0.62), sharey=ax)
    count.set_axis_off()
    n = len(QUANTITIES)
    ax.set_ylim(n - 0.5, -0.5)
    ax.set_xscale("log")
    ax.grid(axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    ax.set_yticks(range(n), [q.short for q in QUANTITIES], fontsize=12)
    ax.tick_params(axis="y", length=0, pad=12)
    ax.set_xlabel("相対的な大きさ", labelpad=7)
    ax.spines["left"].set_visible(False)
    handles = []
    count.set_title("実装値が区間内\n（点数）", fontsize=10.5, weight="bold")
    for i, q in enumerate(QUANTITIES):
        ps = data[q.key]
        for offset, vals, col, marker, title in [
            (-0.15, [p.radius for p in ps], BAND, "|", "保証区間の相対半幅"),
            (
                0.15,
                [abs(p.residual("sashimi")) for p in ps],
                NATIVE,
                "o",
                "SASHIMI-C と中心の相対差（絶対値）",
            ),
        ]:
            x = [float(v) for v in vals if v]
            scatter = ax.scatter(
                x,
                [i + offset] * len(x),
                s=47 if marker == "|" else 17,
                marker=marker,
                lw=1.15,
                alpha=0.9,
                label=title,
                **({"color": col} if marker == "|" else {"facecolor": "none", "edgecolor": col}),
            )
            if i == 0:
                handles.append(scatter)
        inside = sum(abs(p.residual("sashimi")) <= p.radius for p in ps)
        count.text(0.50, i, f"{inside} / {len(ps)}", ha="center", va="center", fontsize=10.5)
        if i < n - 1:
            ax.axhline(i + 0.5, color=GRID, lw=0.6)
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=9))
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.082, 0.862),
        frameon=False,
        ncol=2,
        fontsize=10.5,
    )
    fig.text(
        0.09,
        0.043,
        "厳密な有理数で比較してから描画。ゼロは対数軸から省略し、右欄の点数に含めた。"
        "\nH の区間は倍精度の丸め幅より狭い。右欄は各実装値に関する比較で、Lean 証明の合否は全点で合格。",
        fontsize=9.5,
        color=MUTED,
    )
    return fig


def render(src, out):
    data = load_data(src)
    out.mkdir(exist_ok=True, parents=True)
    tables(out, data)
    style()
    panels = [
        (
            "background",
            curves(
                data,
                "赤方偏移とともに変わる成長と崩壊障壁",
                r"平坦な物質＋Λ　 $H_0=67.4$ km s$^{-1}$ Mpc$^{-1}$　 $\Omega_{m0}=0.315$　 $\Omega_{b0}=0.049$　 放射なし",
                BACKGROUND,
                MODELS,
                "下段は（比較値 / 区間中心 − 1）。H は各ライブラリの膨張率。Astropy の成長は H から積分し、障壁は全実装で 1.686/D に統一。\n線は標本点を結ぶ。保証の対象は計算した各点で、区間の幅は精度図で拡大した。",
            ),
        ),
        (
            "mass",
            curves(
                data,
                "質量分散から条件付き EPS カーネルへ",
                r"質量は物理 $M_\odot$　 分散は $z=0$　 カーネルは親 $(10^{12}\,M_\odot,\,z=0)$ と子 $(M,\,z=1)$ の組",
                MASS,
                ("sashimi", "colossus"),
                "下段は（比較値 / 区間中心 − 1）。ITAMAE と SASHIMI-C は Ludlow フィット。Colossus は Eisenstein–Hu スペクトル＋球形トップハット。\nColossus の入力は $\\sigma_8=0.811$, $n_s=0.965$。F は各実装の S と障壁から求めた Yang 式 (14)。親と同質量の点は ΔS=0。",
            ),
        ),
        (
            "host",
            curves(
                data,
                "EPS に渡すホストハローの質量と成長率",
                r"$M_{200}(0)=10^{12}\,M_\odot$、初期赤方偏移 $z_i=0$。各赤方偏移で SASHIMI-C の質量履歴フィットを評価。",
                HOST,
                ("sashimi",),
                "下段の青い縦線が各点の保証区間。比較値はゼロ付近の丸で示した。\n右図は $\\mathrm{d}\\ln M_{200}/\\mathrm{d}z$ にフィットから求めたビリアル質量を掛ける、SASHIMI-C 固有の dMdz 処方。",
            ),
        ),
        ("precision", accuracy(data)),
    ]
    with PdfPages(out / "eps-comparison.pdf", metadata={"Title": "SASHIMI-C EPS | ITAMAE"}) as pdf:
        for name, fig in panels:
            fig.savefig(out / (name + ".png"), dpi=155)
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    render(args.data_dir, args.output_dir)
