"""Generate the independent z=0 linear CDM spectrum used by public CI.

Run with CAMB 1.6.4, numpy and OMP_NUM_THREADS=1. Write the output to the
filename given on the command line. Input k is in h/Mpc, output P is in
(Mpc/h)^3. The finite band is wholly inside the requested CAMB domain.
"""

import argparse
from pathlib import Path


def generate(path: Path) -> None:
    import camb
    import numpy as np

    if camb.__version__ != "1.6.4":
        raise ValueError("this reference uses camb==1.6.4")
    params = camb.CAMBparams()
    params.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.1200, mnu=0.06, omk=0, tau=0.054)
    params.InitPower.set_params(As=2.1e-9, ns=0.965)
    params.set_matter_power(redshifts=[0], kmax=70)
    results = camb.get_results(params)
    power = results.get_matter_power_interpolator(
        var1="delta_cdm", var2="delta_cdm", hubble_units=True, k_hunit=True, nonlinear=False
    )
    k = np.geomspace(1e-4, 1e2, 3000)
    with path.open("w", encoding="ascii", newline="\n") as output:
        output.write("# Independent CAMB 1.6.4 z=0 CDM auto power; k[h/Mpc] P[(Mpc/h)^3]\n")
        for wave, value in zip(k, power.P(0, k)):
            output.write(f"{wave:.9e} {value:.9e}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    generate(parser.parse_args().output)
