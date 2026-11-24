"""Evaluate independent blocks of prescribed cells and mass radii in Lean."""

from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path
import subprocess

PROJECT = Path(__file__).resolve().parents[1]


def edges(steps):
    left = Fraction(1, 10**8)
    end = Fraction(10**10)
    while left < end:
        right = min(2 * left, end)
        spacing = (right - left) / steps
        for j in range(steps):
            yield (left + j * spacing, left + (j + 1) * spacing)
        left = right


def blocks(seq, jobs):
    size = (len(seq) + jobs - 1) // jobs
    return [seq[i : i + size] for i in range(0, len(seq), size)]


def run(worker, data, bits, degree):
    result = subprocess.run(
        ["lake", "env", "lean", "--run", f"scripts/{worker}.lean", str(bits), str(degree)],
        cwd=PROJECT,
        text=True,
        input=data,
        capture_output=True,
    )
    if result.returncode:
        raise RuntimeError(result.stderr + result.stdout[-1400:])
    return result.stdout


def power_cells(subdivisions, bits, degree, jobs):
    requested = list(edges(subdivisions))

    def evaluate(block):
        data = "".join(f"{a} {b}\n" for a, b in block)
        rows = run("EhFixedWeights", data, bits, degree).splitlines()
        if len(rows) != len(block):
            raise RuntimeError("incomplete cell evaluation")
        for line, pair in zip(rows, block, strict=True):
            tag, a, b, pl, pu, al, au = line.split()
            if (
                tag != "cell"
                or (Fraction(a), Fraction(b)) != pair
                or Fraction(pl) > Fraction(pu)
                or Fraction(al) > Fraction(au)
            ):
                raise RuntimeError("invalid fixed-cell response")
        return "\n".join(rows) + "\n"

    with ThreadPoolExecutor(jobs) as pool:
        return "".join(pool.map(evaluate, blocks(requested, jobs)))


def mass_totals(cells, masses, bits, degree, jobs):
    def evaluate(block):
        data = cells + "".join(f"mass {m}\n" for m in block)
        return run("EhFixedTotals", data, bits, degree).splitlines()

    with ThreadPoolExecutor(jobs) as pool:
        calculated = list(pool.map(evaluate, blocks(masses, jobs)))
    if any(not rows or rows[0] != calculated[0][0] for rows in calculated):
        raise RuntimeError("inconsistent shared normalization anchor")
    return (
        cells
        + calculated[0][0]
        + "\n"
        + "".join(line + "\n" for rows in calculated for line in rows[1:])
    )
