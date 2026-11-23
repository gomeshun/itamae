"""Variance of a tabulated (optionally FDM-modified) spectrum for an interval of radii."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import time

if __package__:
    from .spectrum_table import PROJECT, _read_table, number
    from .table_replay import check_proof, literal
else:
    from spectrum_table import PROJECT, _read_table, number
    from table_replay import check_proof, literal


def replay_window(
    raw_rows,
    layout,
    h,
    model,
    interpolation,
    m22,
    radius,
    band,
    tol,
    depth,
    bits,
    degree,
    lo,
    hi,
    segments,
    chunk_size,
    memory_mb,
    timeout,
):
    """Check independent cell sums, then certify their sum and the reported endpoints."""

    def box(lo, hi):
        return f"⟨{literal(lo)}, {literal(hi)}⟩"

    ratio = (
        f"(fun k => SashimiSpectrum.fdm? {literal(m22)} k {bits} {degree})"
        if model == "fdm"
        else "(fun _ => some (RationalInterval.point 1))"
    )
    mode = ".logLog" if interpolation == "log-log" else ".linear"
    a, b = map(literal, band)
    lines = [
        "import SashimiTableWindow",
        "import SashimiTableRequest",
        "open FloatLib.Numerics SashimiTable SashimiTableWindow",
        "deriving instance DecidableEq for Cell",
        "set_option Elab.async false",
        "set_option maxRecDepth 100000",
        "set_option maxHeartbeats 8000000",
        f"private def power : Knot → Option RationalInterval := weighted? {ratio}",
        f"private def radius : RationalInterval := {box(*radius)}",
        f"private def part := sumCells? 0 {mode} power radius {a} {b} {depth} {bits} {degree}",
    ]
    nodes = []
    chunks = []
    for idx, start in enumerate(range(0, len(raw_rows) - 1, chunk_size)):
        if idx >= len(segments) or segments[idx][0] != start:
            raise RuntimeError("missing or unordered quadrature certificate")
        _, lower, upper = segments[idx]
        knots = ", ".join(box(k, y) for k, y in raw_rows[start : start + chunk_size + 1])
        name = f"N{idx}"
        chunks.append(name)
        lines.extend(
            [
                f"private def rows{name} : List Knot :=",
                f'  (SashimiTableRequest.canonicalRows "{layout}" {literal(h)} [{knots}]).getD []',
                f"private def cells{name} := cells rows{name}",
                f"private theorem cert{name} : part cells{name} = some {box(lower, upper)} := by decide +kernel",
            ]
        )
        nodes.append((name, lower, upper))
    if len(segments) != len(chunks):
        raise RuntimeError("extra quadrature certificates")
    while len(nodes) > 1:
        nxt = []
        for idx in range(0, len(nodes) - 1, 2):
            (x, xl, xu), (y, yl, yu) = nodes[idx : idx + 2]
            name = f"B{len(lines)}"
            lower, upper = xl + yl, xu + yu
            lines.extend(
                [
                    f"private def cells{name} := cells{x} ++ cells{y}",
                    f"private theorem cert{name} : part cells{name} = some {box(lower, upper)} :=",
                    f"  (sumCells?_append_of cert{x} cert{y}).trans (by decide +kernel)",
                ]
            )
            nxt.append((name, lower, upper))
        if len(nodes) % 2:
            nxt.append(nodes[-1])
        nodes = nxt
    name, _, _ = nodes[0]
    # Each pair of chunks shares precisely its boundary knot. The definitionally
    # equal cell lists in the 'hs' proof below bind every pair to the raw rows.
    lines.append(
        "private def knots := "
        + " ++ ".join([f"rows{chunks[0]}"] + [f"rows{c}.tail" for c in chunks[1:]])
    )
    lines.extend(
        [
            f"theorem concrete : integrate? 0 {mode} power knots radius",
            f"    {a} {b} {literal(tol)} {depth} {bits} {degree} = some {box(lo, hi)} := by",
            "  apply integrate_of_parts (by decide +kernel)",
            f"  · rw [show cells knots = cells{name} by first | rfl | decide +kernel]",
            f"    exact cert{name}",
            "  · decide +kernel",
            "#print axioms concrete",
        ]
    )
    check_proof("\n".join(lines) + "\n", PROJECT, max_memory_mb=memory_mb, timeout_seconds=timeout)


def evaluate_window(
    table,
    *,
    layout,
    h,
    model,
    m22,
    radius_lower,
    radius_upper,
    tolerance,
    interpolation="log-log",
    k_min=None,
    k_max=None,
    depth=0,
    bits=80,
    degree=12,
    replay=False,
    replay_memory_mb=4096,
    replay_chunk_cells=16,
    replay_timeout=1800,
):
    if model not in ("cdm", "fdm") or interpolation not in ("log-log", "linear"):
        raise ValueError("unknown transfer or interpolation")
    if (
        any(
            type(i) is not int
            for i in (bits, degree, depth, replay_memory_mb, replay_chunk_cells, replay_timeout)
        )
        or min(replay_memory_mb, replay_chunk_cells, replay_timeout) <= 0
    ):
        raise ValueError(
            "depth, bits, degree and replay memory must be integers with positive replay memory"
        )
    if type(replay) is not bool:
        raise ValueError("replay must be a boolean")
    m22, radius_lower, radius_upper, tolerance = map(
        number, (m22, radius_lower, radius_upper, tolerance)
    )
    k_min = None if k_min is None else number(k_min)
    k_max = None if k_max is None else number(k_max)
    rows, metadata, raw_rows = _read_table(Path(table), layout, h)
    a = rows[0][0] if k_min is None else k_min
    b = rows[-1][0] if k_max is None else k_max
    if not (rows[0][0] <= a <= b <= rows[-1][0] and 0 <= radius_lower <= radius_upper):
        raise ValueError("invalid integration band or radius interval")
    if min(tolerance, depth, bits, degree) < 0 or (model == "fdm" and m22 <= 0):
        raise ValueError("invalid precision, tolerance or FDM mass")

    def parts(*values):
        return [str(i) for value in values for i in (value.numerator, value.denominator)]

    args = [
        model,
        interpolation,
        *parts(m22, radius_lower, radius_upper, a, b, tolerance),
        str(depth),
        str(bits),
        str(degree),
        str(replay_chunk_cells if replay else 0),
    ]
    started = time.monotonic()
    run = subprocess.run(
        ["lake", "env", "lean", "--run", "scripts/TableWindow.lean", *args],
        cwd=PROJECT,
        input="".join(" ".join(parts(*row)) + "\n" for row in rows),
        text=True,
        capture_output=True,
    )
    if run.returncode:
        raise RuntimeError((run.stdout + run.stderr)[-3000:])
    result = dict(
        table=metadata,
        power_transfer=model,
        interpolation=interpolation,
        transfer_placement="at-knots",
        radius_Mpc=list(map(str, (radius_lower, radius_upper))),
        k_band_Mpc_inverse=list(map(str, (a, b))),
        unit="1",
        absolute_width_budget=str(tolerance),
        bits=bits,
        degree=degree,
        depth=depth,
        runtime_seconds=time.monotonic() - started,
    )
    if model == "fdm":
        result["m22"] = str(m22)
    fields = run.stdout.strip().splitlines()[-1].split()
    if fields == ["unresolved"]:
        return dict(result, status="unresolved")
    if len(fields) != 3 or fields[0] != "enclosed":
        raise RuntimeError(run.stdout)
    lo, hi = map(Fraction, fields[1:])
    if not lo <= hi or hi - lo > tolerance:
        raise RuntimeError("inconsistent Lean output")
    result["kernel_replay"] = "not-requested"
    if replay:
        started = time.monotonic()
        segments = []
        for line in run.stdout.splitlines():
            fields = line.split()
            if fields and fields[0] == "part":
                segments.append((int(fields[1]), Fraction(fields[2]), Fraction(fields[3])))
        replay_window(
            raw_rows,
            layout,
            number(metadata["h"]),
            model,
            interpolation,
            m22,
            (radius_lower, radius_upper),
            (a, b),
            tolerance,
            depth,
            bits,
            degree,
            lo,
            hi,
            segments,
            replay_chunk_cells,
            replay_memory_mb,
            replay_timeout,
        )
        result.update(
            kernel_replay="passed",
            replayed_chunks=len(segments),
            replay_seconds=time.monotonic() - started,
        )
    return dict(result, status="enclosed", lower=str(lo), upper=str(hi), width=str(hi - lo))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("table", type=Path)
    ap.add_argument(
        "--format", dest="layout", choices=["canonical", "camb", "wmap7"], required=True
    )
    ap.add_argument("--h", type=Fraction)
    ap.add_argument("--model", choices=["cdm", "fdm"], required=True)
    ap.add_argument("--m22", type=Fraction, default=Fraction(1))
    ap.add_argument("--k-min", type=Fraction)
    ap.add_argument("--k-max", type=Fraction)
    ap.add_argument("--radius-lower", type=Fraction, required=True)
    ap.add_argument("--radius-upper", type=Fraction, required=True)
    ap.add_argument("--tolerance", type=Fraction, required=True)
    ap.add_argument("--interpolation", choices=["linear", "log-log"], default="log-log")
    ap.add_argument("--depth", type=int, default=0)
    ap.add_argument("--bits", type=int, default=80)
    ap.add_argument("--degree", type=int, default=12)
    ap.add_argument(
        "--replay", action="store_true", help="also reduce the full computation in the kernel"
    )
    ap.add_argument("--replay-memory-mb", type=int, default=4096)
    ap.add_argument("--replay-chunk-cells", type=int, default=16)
    ap.add_argument(
        "--replay-timeout", type=int, default=1800, help="total kernel replay budget in seconds"
    )
    result = evaluate_window(**vars(ap.parse_args()))
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "enclosed" else 2)
