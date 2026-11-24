"""Time the matched mass-grid numerical workloads on the same machine.

After lake build, use a Python environment containing the reference's
colossus==1.4.0, NumPy and SciPy: python scripts/eh_timings.py --output benchmark.json
"""

import argparse
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import tempfile
import time

from colossus.cosmology import cosmology
from colossus.lss import peaks
import numpy as np

from eh_plot_variance import PROJECT, _read_output, _source
from eps_eh_reference import PARAMS, PS


def lean_pass(masses, mesh):
    program = (
        _source(masses, [mesh], entrypoint="runWorkload")
        + """
def main : IO Unit := do
  let start ← IO.monoNanosNow
  runWorkload
  let stop ← IO.monoNanosNow
  IO.eprintln s!"elapsed_ns {stop - start}"
"""
    )
    with tempfile.TemporaryDirectory() as tmp:
        request = Path(tmp) / "Timing.lean"
        request.write_text(program)
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        start = time.perf_counter()
        run = subprocess.run(
            ["lake", "env", "lean", "--run", str(request)],
            cwd=PROJECT,
            capture_output=True,
            text=True,
            check=True,
        )
        wall = time.perf_counter() - start
        end = resource.getrusage(resource.RUSAGE_CHILDREN)
    (result,) = _read_output(run.stdout, [mesh], masses)
    (elapsed,) = [
        line.split()[1] for line in run.stderr.splitlines() if line.startswith("elapsed_ns ")
    ]
    timing = {
        "workload_wall_seconds": int(elapsed) / 1e9,
        "process_wall_seconds": wall,
        "process_cpu_seconds": end.ru_utime + end.ru_stime - usage.ru_utime - usage.ru_stime,
        "scope": "Lean binary64, one Simpson mesh plus R8, radii, S and dS/dM and row output. Process timer also includes Lake, imports and elaboration with prebuilt oleans.",
    }
    return result["rows"], timing


def spline_pass(masses, repeats):
    cosmo = cosmology.Cosmology(name="eps_eh_spline", **PARAMS)
    previous = cosmology.current_cosmo
    cosmology.setCurrent(cosmo)
    try:
        start_cpu, start = time.process_time(), time.perf_counter()
        # Isolate the same first-use log-sigma table builder that public sigma() calls.
        cosmo._sigmaInterpolator(0, "tophat", False, None, None, False, PS)
        build = time.perf_counter() - start
        build_cpu = time.process_time() - start_cpu
        times = []
        start_cpu = time.process_time()
        for _ in range(repeats):
            start = time.perf_counter()
            radii = peaks.lagrangianR(masses * cosmo.h)
            s = cosmo.sigma(radii, filt="tophat", ps_args=PS) ** 2
            slope = cosmo.sigma(radii, filt="tophat", derivative=True, ps_args=PS)
            grad = 2 * s * slope / (3 * masses)
            times.append(time.perf_counter() - start)
        query_cpu = (time.process_time() - start_cpu) / repeats
    finally:
        cosmology.setCurrent(previous)
    return np.column_stack((s, grad)), {
        "cold_table_wall_seconds": build,
        "cold_table_cpu_seconds": build_cpu,
        "cached_batch_wall_median_seconds": statistics.median(times),
        "cached_batch_wall_min_seconds": min(times),
        "cached_batch_cpu_mean_seconds": query_cpu,
        "cached_repetitions": repeats,
        "scope": "Default z=0 top-hat spline. Fresh Cosmology with disk persistence disabled; table timing includes required power-spectrum construction. Cached array query includes mass-to-radius, S and spline dS/dM. Package imports and Cosmology constructor excluded.",
    }


def relative_difference(measured, expected):
    return dict(zip(("S", "dS_dM"), np.max(np.abs((measured - expected) / expected), axis=0)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, default=PROJECT / "results/eps-eh")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--repeats", type=int, default=100)
    args = p.parse_args()
    if args.repeats < 1:
        p.error("--repeats must be positive")
    lean = json.loads((args.results / "lean.json").read_text())
    reference = json.loads((args.results / "colossus.json").read_text())
    versions = {name: version(name) for name in ("colossus", "numpy", "scipy")}
    if versions["colossus"] != "1.4.0" or reference["metadata"]["parameters"] != PARAMS:
        raise ValueError("expected pinned Colossus reference cosmology")
    rows = lean["mass_rows"] + [lean["host_progenitor"]]
    spline = reference["variance"] + [reference["host_progenitor"]]
    masses = [float(r["mass_msun"]) for r in rows]
    if masses != [float(r["mass_msun"]) for r in spline]:
        raise ValueError("reference mass grids differ")
    ppd = lean["quadrature"]["panels_per_decade"]
    band = lean["quadrature"]["k_band_h_per_Mpc"]
    low, high = np.log10(band)
    if low != int(low) or high != int(high):
        raise ValueError("expected integer-decade mesh limits")
    load_start = os.getloadavg()
    computed, lt = lean_pass(masses[:-1], (int(low), int(high), ppd))
    if [r["mass_msun"] for r in computed] != masses:
        raise ValueError("Lean progenitor differs from the requested mass")
    sc, st = spline_pass(np.asarray(masses), args.repeats)

    def values(data):
        return np.array([[r["S"], r["dS_dM"]] for r in data])

    errors = {
        "lean_vs_stored_lean": relative_difference(values(computed), values(rows)),
        "spline_vs_stored_spline": relative_difference(sc, values([r["spline"] for r in spline])),
    }
    if any(v > 1e-12 or not np.isfinite(v) for method in errors.values() for v in method.values()):
        raise ValueError(f"reproduction of stored results failed: {errors}")
    info = Path("/proc/cpuinfo")
    model = (
        next(
            (
                line.partition(":")[2].strip()
                for line in info.read_text().splitlines()
                if line.startswith("model name")
            ),
            platform.machine(),
        )
        if info.exists()
        else platform.machine()
    )
    quota = Path("/sys/fs/cgroup/cpu.max")
    report = {
        "workload": {
            "grid_masses": len(masses) - 1,
            "extra_progenitor_msun": masses[-1],
            "mass_range_msun": [min(masses), max(masses)],
            "redshift": 0,
            "k_band_h_per_Mpc": band,
            "panels_per_decade": ppd,
            "panels": (int(high) - int(low)) * ppd,
        },
        "environment": {
            "cpu": model,
            "available_affinity_cpus": len(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else os.cpu_count(),
            "cgroup_cpu_max": quota.read_text().strip() if quota.exists() else None,
            "load_average_start": load_start,
            "load_average_end": os.getloadavg(),
            "python": platform.python_version(),
            "packages": versions,
            "lean": subprocess.check_output(
                ["lake", "env", "lean", "--version"], cwd=PROJECT, text=True
            ).strip(),
            "measurement": "Single host, unisolated scheduling; concurrent proof workers can affect wall time.",
        },
        "lean": lt,
        "colossus": st,
        "max_relative_reproduction_error": errors,
        "scope": "Floating-point numerical predictions. Kernel proof checking and the initial Lake build are separate workloads and excluded.",
    }
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
