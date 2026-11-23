# Spectrum fixtures

`wmap7-cdm.dat` is the unchanged 3,001-row SASHIMI-W input from Shin'ichiro Ando's public repository at `99dfc3632eec0126080c0273ddf78f84fe09216c`. The five scalar header entries precede columns P and k. The SHA-256 is `417a0e8e9feafae8c6afb4b89e63aef80718c9641ae95994d25a6cd97438f8c9`. The same license accompanies this file and its two-knot excerpt. [Upstream data](https://github.com/shinichiroando/sashimi-w/blob/99dfc3632eec0126080c0273ddf78f84fe09216c/WMAP7_camb_matterpower_z0_extrapolated.dat) · [MIT license](wmap7-LICENSE.txt)

`camb_cdm0.dat` is an independent linear CDM auto-power calculation for public regression tests, with 3,000 logarithmically spaced wavenumbers from 0.0001 to 100 h/Mpc. The input is Planck-like with h=0.674. All knots lie inside the requested CAMB domain. The release's family-table reproduction uses the distinct SASHIMI-F hash `91269da31616bf5f8d506ad138f157c8eceb284d4cd328cd99e730b882a8dd61` and requires access to that source table. [CAMB generator and parameters](../../scripts/make_camb_reference.py) · [Family-table results](../../results/sashimi-table-window.json)

Regenerate the independent fixture in an environment with CAMB 1.6.4 and NumPy, from the Lake project directory.

```sh
OMP_NUM_THREADS=1 python3 scripts/make_camb_reference.py /tmp/camb_cdm0.dat
sha256sum /tmp/camb_cdm0.dat tests/data/camb_cdm0.dat
```

`camb_cdm31.csv` contains dimensionless power at three wavenumbers from a separate z=31 calculation. `fdm-camb-three-knots.dat` contains a short integration test input in CAMB units. [Prompt-cusp generation and validation](../../CUSP.md)
