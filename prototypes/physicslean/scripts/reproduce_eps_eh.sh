#!/usr/bin/env bash
# Rebuild the matched-model figures and replay their main analytic grid.
set -euo pipefail
out=$(realpath -m "${1:?Usage: scripts/reproduce_eps_eh.sh OUTPUT_DIRECTORY}")
jobs=${JOBS:-2}
project_dir=$(cd "$(dirname "$0")/.." && pwd)
cd "$project_dir"
mkdir -p "$out"
grid=results/eps-plots/grid.json
host="$out/proofs/host"
masses="$out/proofs/masses"
slopes="$out/proofs/slopes"
pairs="$out/proofs/pairs"
lake build
python scripts/eh_fixed_block_replay.py "$out/fixed-octave.lean" "$out/fixed-octave.json" --jobs "$jobs"
python scripts/eps_eh_reference.py "$grid" "$out/colossus.json"
python scripts/eh_transfer.py "$out/colossus.json" "$out/transfer.json"
python scripts/eh_transfer_samples.py "$out/transfer-samples.json"
python scripts/eps_eh_growth.py "$out/growth.json"
python scripts/eh_plot_variance.py --grid "$grid" --output "$out/lean.json"
python scripts/eh_timings.py --results "$out" --output "$out/benchmark.json"
python scripts/eps_eh_downstream.py "$out/lean.json" "$out/growth.json" "$out/colossus.json" "$out/downstream.json"
for subdivisions in 16 128; do
    suffix=
    if [[ "$subdivisions" == 128 ]]; then suffix=-fine; fi
    python scripts/eh_sweep.py --subdivisions "$subdivisions" --grid "$grid" --numerical "$out/lean.json" \
        --output "$out/sweep$suffix.json" --transport "$out/sweep-bins$suffix.json"
    python scripts/eh_derivative_bounds.py "$out/sweep-bins$suffix.json" "$grid" "$out/derivative$suffix.json"
done
python scripts/eh_sweep.py --subdivisions 128 --grid "$grid" --numerical "$out/lean.json" \
    --output "$out/sweep-tight.json" --transport "$out/sweep-bins-tight.json" \
    --fixed-cell-slopes "$out/derivative-tight.json" --jobs "$jobs"
python scripts/eh_sweep_certificate.py "$out/sweep-bins.json" "$out/sweep.json" --mass 1e12 \
    --output "$out/kernel-host.json" --module-output "$host" --jobs "$jobs"
python scripts/eh_sweep_masses.py --sweep "$out/sweep.json" --modules "$host" --output "$masses" --jobs "$jobs" \
    > "$out/mass-replay.json"
python scripts/eh_full_grid.py "$out/sweep.json" "$out/full.json" --mass-modules "$masses" --anchor-modules "$host"
python scripts/eh_full_grid.py "$out/sweep-fine.json" "$out/full-fine.json"
python scripts/eh_derivative_certificate.py prepare "$slopes" --mesh "$out/sweep-bins.json" --derivative "$out/derivative.json"
python scripts/eh_derivative_certificate.py panels "$slopes" --jobs "$jobs"
python scripts/eh_derivative_certificate.py finish "$slopes" --shared "$host" --jobs "$jobs"
python scripts/eh_derivative_certificate.py report "$slopes" --output "$out/kernel-derivative.json"
python scripts/eh_full_derivative.py "$out/sweep.json" "$out/derivative-full.json" --shared "$host" --slopes "$slopes" --certificates "$out/proofs/full-slopes"
python scripts/eh_full_derivative.py "$out/sweep-fine.json" "$out/derivative-full-fine.json" --finite-evaluation "$out/derivative-fine.json" --certificates "$out/proofs/full-slopes-fine"
python scripts/eh_paired_windows.py "$out/sweep.json" "$host" "$pairs" "$out/paired-children.json" --indices 71 70 69 --jobs "$jobs"
python scripts/eh_paired_windows.py "$out/sweep.json" "$host" "$pairs" "$out/paired-73.json" --indices 73 --jobs "$jobs"
python scripts/eh_host_bounds.py "$out/sweep-fine.json" "$out/growth.json" "$out/host-bounds-fine.json"
python scripts/eh_host_bounds.py "$out/sweep.json" "$out/growth.json" "$out/host-bounds.json" --paired "$out/paired-73.json" --modules "$pairs" "$masses" "$host"
python scripts/eps_eh_kernel_bounds.py "$out/full.json" "$out/growth.json" "$out/kernel-bounds.json" --paired "$out/paired-children.json" --modules "$pairs" "$masses" "$host"
python scripts/eps_eh_kernel_bounds.py "$out/full-fine.json" "$out/growth.json" "$out/kernel-bounds-fine.json"
python scripts/eh_tight_downstream.py "$out/sweep-tight.json" "$out/growth.json" "$out" \
    --reference-sweep "$out/sweep-fine.json" --finite-slope "$out/derivative-tight.json"
python scripts/plot_eh_match.py "$out" results/eps-plots "$out"
python scripts/eh_width_certificate.py --csv "$out/matched-samples.csv" --output "$out/fixed-widths.lean"
