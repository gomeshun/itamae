"""Compare native transfer curves using identical binary64 k inputs."""

import argparse
import json
from pathlib import Path
import struct
import subprocess

from colossus.cosmology import power_spectrum
import numpy as np

if __package__:
    from .eps_eh_reference import PARAMS
else:
    from eps_eh_reference import PARAMS


def calculate():
    # Include a denser sampling through the baryon acoustic scales.
    wave = np.unique(np.concatenate([np.geomspace(1e-6, 1e8, 337), np.geomspace(0.015, 0.5, 161)]))
    encoded = "".join(str(int.from_bytes(struct.pack(">d", k), "big")) + "\n" for k in wave)
    run = subprocess.run(
        ["lake", "env", "lean", "--run", "scripts/EhTransferSamples.lean"],
        cwd=Path(__file__).resolve().parents[1],
        input=encoded,
        capture_output=True,
        text=True,
        check=True,
    )
    computed = [
        struct.unpack(">d", int(bits).to_bytes(8, "big"))[0] for bits in run.stdout.splitlines()
    ]
    h = PARAMS["H0"] / 100
    samples = []
    for k, lean in zip(wave, computed, strict=True):
        col = float(
            power_spectrum.modelEisenstein98(k, h, PARAMS["Om0"], PARAMS["Ob0"], PARAMS["Tcmb0"])
        )
        samples.append(dict(k_h_mpc=float(k), lean=lean, colossus=col))
    return dict(
        target="Native binary64 evaluations at the same input bits; sampled curve", rows=samples
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(), indent=2, allow_nan=False) + "\n")
