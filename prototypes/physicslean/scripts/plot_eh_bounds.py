"""Plot exact rational enclosures separately from differences between implementations."""

import csv
from fractions import Fraction as Q

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullLocator

if __package__:
    from .plot_eps import MASS, GRID, INK, MUTED
else:
    from plot_eps import MASS, GRID, INK, MUTED


def intervals(rows, data):
    """Join each mass to its enclosure, keeping unresolved EPS pairs in the output."""
    sources = {name: {Q(str(r["mass_msun"])): r for r in doc["rows"]} for name, doc in data.items()}
    records = []
    derivative = "derivative-full" if "derivative-full" in sources else "derivative"
    for q, source in zip(MASS, ["full", "full", derivative, "kernel-bounds"], strict=True):
        for p in (row for row in rows if row["quantity"] == q.key):
            r = sources[source][Q(p["coordinate"])]
            fine = sources.get(source + "-fine", {}).get(Q(p["coordinate"]), {})
            tight = fine.get(q.key)
            box = r.get(q.key)
            if box and Q(box["lo"]) > Q(box["hi"]):
                raise ValueError(f"Inverted {q.key} interval at {p['coordinate']}")
            target = {
                "full": "0 < k < infinity",
                "derivative": "finite band",
                "derivative-full": "0 < k < infinity",
                "kernel-bounds": data["kernel-bounds"]["variance_target"],
            }[source]
            closed = r.get("end_to_end_kernel_replayed", r.get("kernel_replayed", False))
            if source == "derivative" and "kernel-derivative" in sources:
                proof = sources["kernel-derivative"][Q(p["coordinate"])]
                if proof["dS_dM"] != box:
                    raise ValueError("Mismatched derivative proof endpoints")
                closed = proof["kernel_replayed"]
            records.append(
                dict(
                    quantity=q.key,
                    mass_msun=p["coordinate"],
                    target=target,
                    inputs=source + ".json",
                    native=p["matched"],
                    colossus=p["colossus_spline"],
                    fixed_band=p["fixed_band"],
                    lo=box["lo"] if box else "",
                    hi=box["hi"] if box else "",
                    status=r.get("status", "enclosed"),
                    end_to_end_kernel_replay=closed,
                    fine_lo=tight["lo"] if tight else "",
                    fine_hi=tight["hi"] if tight else "",
                    fine_status=fine.get("status", "enclosed" if tight else "unavailable"),
                )
            )
    return records


def save_and_plot(records, data, output):
    with (output / "enclosures.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=records[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9.5))
    fig.subplots_adjust(left=0.105, right=0.98, bottom=0.14, top=0.81, hspace=0.39, wspace=0.24)
    fig.text(0.085, 0.964, "質量ごとの実数積分の包含区間", fontsize=20, weight="bold", va="top")
    fig.text(
        0.085,
        0.91,
        "各標本で波数ゼロから無限大までの積分を囲む。縦線の長さは、その質量での数値保証の幅。"
        if "derivative-full" in data
        else "分散と標準偏差は全波数を積分。各パネルに評価範囲を表示。",
        fontsize=10.5,
        color=MUTED,
    )
    blue, dark, orange = "#74aac4", "#163d59", "#bd6209"
    handles = [
        Line2D([], [], color=blue, lw=2, label="960 分割の区間"),
        Line2D([], [], color="#176a61", lw=3, label="7,680 分割の区間"),
        Line2D([], [], color=dark, lw=1.8, label="Lean 数値積分"),
        Line2D([], [], ls="none", marker="o", color=orange, mfc="none", label="Colossus 既定"),
    ]
    if any(not r["lo"] for r in records):
        handles.append(
            Line2D([], [], ls="none", marker="x", color="#9e4856", label="960 分割で未確定")
        )
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.08, 0.891),
        ncol=len(handles),
        frameon=False,
        fontsize=10.3,
    )
    for ax, q in zip(axes.flat, MASS, strict=True):
        group = [r for r in records if r["quantity"] == q.key]
        xx = [float(Q(r["mass_msun"])) for r in group]
        vx, low, high, missing = [], [], [], []
        for x, r in zip(xx, group, strict=True):
            if not r["lo"]:
                missing.append(x)
                continue
            lo, hi = sorted(float(Q(r[k])) * q.factor for k in ("lo", "hi"))
            vx.append(x)
            low.append(lo)
            high.append(hi)
        ax.set_xscale("log")
        ax.set_yscale("log")
        # A mathematically positive density with a zero lower bound extends
        # below any drawable logarithmic ordinate; keep its exact bound in CSV.
        ax.plot(xx, [r["native"] * q.factor for r in group], color=dark, lw=1.7)
        ax.plot(
            xx[::4],
            [r["colossus"] * q.factor for r in group[::4]],
            ls="none",
            marker="o",
            color=orange,
            mfc="none",
            ms=4,
        )
        positive = [v for v in low if v > 0]
        floor = min(positive + [r["native"] * q.factor for r in group]) * 0.7
        ax.vlines(vx, [max(floor, value) for value in low], high, color=blue, lw=1.35, zorder=0)
        fine = [r for r in group if r["fine_lo"]]
        if fine:
            fboxes = [
                sorted(float(Q(r[k])) * q.factor for k in ("fine_lo", "fine_hi")) for r in fine
            ]
            ax.vlines(
                [float(Q(r["mass_msun"])) for r in fine],
                [max(floor, pair[0]) for pair in fboxes],
                [pair[1] for pair in fboxes],
                color="#176a61",
                lw=2.1,
                zorder=1,
            )
        if missing:
            ax.plot(
                missing,
                [0.055] * len(missing),
                transform=ax.get_xaxis_transform(),
                marker="x",
                ls="none",
                color="#9e4856",
                mew=1.6,
                clip_on=False,
            )
        scope = (
            "全波数"
            if q.key in ("sigma", "S") or (q.key == "dS_dM" and "derivative-full" in data)
            else "有限帯域"
        )
        if (
            q.key == "F"
            and data["kernel-bounds"].get("integration_domain_k_h_per_Mpc") == "(0,infinity)"
        ):
            scope = "全波数"
        ax.set_title(f"{q.title}　·　{scope}", loc="left", fontsize=12.3, weight="bold")
        ax.set_xlabel(r"質量 $M$  [$M_\odot$]")
        ax.set_ylabel(q.unit)
        ax.xaxis.set_major_locator(LogLocator(numticks=7))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.grid(axis="y", color=GRID, lw=0.65)
        ax.set_axisbelow(True)
    fig.text(
        0.085,
        0.061,
        "縦線は実数モデルの評価器が返した上下限。標本間を結ぶ紺線は数値積分値。物理フィットの誤差は含まない。",
        fontsize=10.2,
        color=INK,
    )
    fig.text(
        0.085,
        0.032,
        "個別の証明をカーネルで再検査した範囲と、評価器を実行した範囲を添付の検証記録に掲載。",
        fontsize=9.4,
        color=MUTED,
    )
    return fig


def error_budget(records):
    """Compare observed implementation differences to enclosure-derived upper bounds."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.7))
    fig.subplots_adjust(left=0.105, right=0.98, top=0.802, bottom=0.125, hspace=0.385, wspace=0.24)
    fig.text(
        0.086,
        0.969,
        "保証できた数値誤差と、実装間で観測した差",
        fontsize=20,
        weight="bold",
        va="top",
    )
    fig.text(
        0.086,
        0.909,
        "縦軸は Lean 数値積分値に対する絶対相対差。上限は区間の両端のうち、数値積分値から遠い方で決まる。",
        fontsize=10.4,
        color=MUTED,
    )
    colors = ["#658aa7", "#176a61", "#bd6209", "#826380"]
    labels = [
        "960 分割・誤差上限",
        "7,680 分割・誤差上限",
        "Colossus 既定との差",
        "同一帯域・再規格化の差",
    ]
    handles = [
        Line2D([], [], color=c, lw=1.8, ls=ls, label=lab)
        for c, lab, ls in zip(colors, labels, ["-", "-", "--", ":"], strict=True)
    ]
    fig.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.083, 0.884),
        ncol=4,
        frameon=False,
        fontsize=9.9,
    )
    for ax, q in zip(axes.flat, MASS, strict=True):
        subset = [r for r in records if r["quantity"] == q.key]

        def series(r):
            native = Q(r["native"])

            def rel(v):
                return abs((Q(v) - native) / native)

            def bound(lo, hi):
                return max(rel(r[lo]), rel(r[hi])) if r[lo] else None

            return [
                bound("lo", "hi"),
                bound("fine_lo", "fine_hi"),
                rel(r["colossus"]),
                rel(r["fixed_band"]),
            ]

        observations = [series(r) for r in subset]
        for k, (c, ls) in enumerate(zip(colors, ["-", "-", "--", ":"], strict=True)):
            drawable = [
                (float(Q(row["mass_msun"])), float(vals[k]))
                for row, vals in zip(subset, observations, strict=True)
                if vals[k]
            ]
            if drawable:
                x, y = zip(*drawable)
                ax.loglog(x, y, color=c, ls=ls, lw=1.65)
        ax.set_title(q.title, loc="left", fontsize=12.7, weight="bold")
        ax.set_xlabel(r"質量 $M$  [$M_\odot$]")
        ax.set_ylabel("相対値")
        ax.xaxis.set_major_locator(LogLocator(numticks=7))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_locator(LogLocator(numticks=6))
        ax.grid(axis="y", color=GRID, lw=0.65)
    fig.text(
        0.086,
        0.047,
        "実線は各標本で max(|U−v|, |v−L|)/|v|。破線と点線は比較実装値と v の観測差。線は標本間を結ぶ。",
        fontsize=9.7,
        color=MUTED,
    )
    return fig


def point_proofs(transfer, growth, colossus, baseline, curve, background):
    """Plot the kernel-checked physical values and their exact relative widths."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.3), gridspec_kw={"height_ratios": [1.52, 1]})
    fig.subplots_adjust(left=0.065, right=0.975, bottom=0.125, top=0.80, wspace=0.32, hspace=0.38)
    fig.text(0.085, 0.965, "伝達関数・膨張率・崩壊障壁の精度", fontsize=20, weight="bold", va="top")
    fig.text(
        0.085,
        0.90,
        "伝達関数の10波数と成長の29赤方偏移で Lean カーネル検査。下段は厳密な有理数端点から求めた相対区間幅。",
        fontsize=10.5,
        color=MUTED,
    )
    cu = {Q(str(row["k_h_mpc"])): row["T"] for row in colossus["transfer"]}
    for i, (rs, key, coordinate, title, xlab, ylab) in enumerate(
        [
            (
                transfer["rows"],
                "T",
                "k_h_mpc",
                "バリオンを含む伝達関数",
                r"波数 $k$  [$h$/Mpc]",
                "$T(k)$",
            ),
            (
                background["rows"],
                "H",
                "z",
                "Physlib の膨張率",
                "赤方偏移 $z$",
                r"$H(z)$  [km s$^{-1}$ Mpc$^{-1}$]",
            ),
            (
                growth["rows"],
                "delta_c",
                "z",
                "球対称崩壊の障壁",
                "赤方偏移 $z$",
                r"$\delta_c(z)=1.686/D(z)$",
            ),
        ]
    ):
        top, bottom = axes[:, i]
        xs, center, width, ref, zeros = [], [], [], [], []
        cb = {Q(str(r["z"])): r["spline"].get(key) for r in colossus["background"]}
        for r in rs:
            x = Q(r[coordinate])
            b = r if key == "T" else r[key]
            lo, hi = Q(b["lo"]), Q(b["hi"])
            c = (lo + hi) / 2
            xs.append(float(x))
            center.append(float(c))
            ref.append((cu if key == "T" else cb)[x])
            width.append(float((hi - lo) / abs(c)))
            if lo == hi:
                zeros.append(float(x))
        if key == "T":
            for a in (top, bottom):
                a.set_xscale("log")
            top.set_yscale("log")
            wave = [r["k_h_mpc"] for r in curve["rows"]]
            top.plot(
                wave,
                [r["lean"] for r in curve["rows"]],
                lw=1.45,
                color="#163d59",
                label="Lean 数値線",
            )
        top.plot(
            xs,
            center,
            "o" if key == "T" else "o-",
            color="#163d59",
            ms=3.4,
            lw=1.7,
            label="Lean 区間中心",
        )
        top.plot(xs, ref, "o", mfc="none", color="#bd6209", ms=6, label="Colossus")
        if key != "T":
            top.plot(
                xs,
                [r[key] for r in baseline["astropy"]["background"]],
                "1",
                markevery=(1, 4),
                color="#315c39",
                ms=8,
                label="Astropy" if key == "H" else "Astropy + 成長積分",
            )
        top.set_title(title, loc="left", weight="bold", fontsize=13)
        top.set_ylabel(ylab)
        top.legend(frameon=False, fontsize=8.8)
        bottom.set_yscale("log")
        bottom.plot(
            [x for x, w in zip(xs, width) if w],
            [w for w in width if w],
            "o-",
            color="#00817e",
            ms=3.4,
            lw=1,
        )
        if zeros:
            bottom.text(
                0.025,
                0.10,
                "幅ゼロの標本は対数軸外",
                transform=bottom.transAxes,
                fontsize=9,
                color=MUTED,
            )
        bottom.set_ylabel("相対区間幅 $(U-L)/|c|$")
        for a in (top, bottom):
            a.set_xlabel(xlab)
            a.grid(axis="y", color=GRID, lw=0.65)
            a.set_axisbelow(True)
    fig.text(
        0.085,
        0.045,
        "伝達関数の紺線はネイティブ演算で計算。高波数側では、指定した Eisenstein–Hu 式の外挿を評価。",
        fontsize=9.5,
        color=MUTED,
    )
    return fig
