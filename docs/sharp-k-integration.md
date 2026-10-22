# Stable moving-boundary sharp-k integration

Sharp-k variance is the integral of `k^3 P(k)/(2 pi^2)` in log-wavenumber up to
the mass-dependent physical cutoff, truncated at the configured finite power
domain. The previous integrator stretched the entire Simpson grid whenever the
upper cutoff changed. This moved samples across the tabulated spectrum's
features and caused fluctuations even where physical variance is nearly flat.

Specification `integrated-variance:v3` uses `n_k-1` fixed logarithmic cells and
the union of those cells with optional spectrum interpolation breakpoints, and
Gauss-Legendre integration of order `sharp_k_order` (default 8) within each cell.
Complete-cell contributions are accumulated once per call; only the final
partial cell depends on the requested mass cutoff. All quadrature weights are
positive. There is no projection of variance, clipping of a derivative or
modification of statistical weights. Power evaluation coordinates are limited
to the configured domain solely to avoid one-ulp exponentiation excursions.

The sharp-k derivative remains the physical moving-boundary term
`-k_c^3 P(k_c) D(z)^2/(6 pi^2 M)` inside the finite domain and zero outside it.
Smooth windows retain their Simpson integration and logarithmic finite
difference. Both cell resolution and quadrature order are in the new identity;
v2 cache identities must not be reused. Batch size changes allocation only.

## Independent evidence

For `P(k)=exp(-k)`, the finite-domain integral is an incomplete gamma function.
Two protection cases fail before this change and pass afterward: the cutoff
plateau/integral and consistency of the moving-boundary derivative with an
independent finite difference. They cover 2,000 masses, scalar/array and
partition equivalence, and 257/513-cell-node convergence. A third protection case compares a piecewise-linear spectrum with its exact
quartic primitive; without splitting at its knots the discrepancy was 4.37e-5.
The complete core suite has 106 passing cases, with Ruff/mypy passing and branch coverage 80%.

The WMAP7 table was additionally integrated by adaptive QUADPACK on its own
3,001 original knots, independently of the product's fixed log-k cells. For the
explicit q5 and q10 spectra at particle mass 2 keV, the final union-of-knots
integrator gives a maximum relative discrepancy of 1.12e-15 over 13
representative masses. Changing from 513 to 8,193 fixed grid nodes changes
variance by at most 2.0e-15. The former grid produced positive variance steps
of up to 8.66e-4 across a 1,000-mass scan; the final scan has no positive steps.
Independent finite differences at redshift 0.5 and physical masses 1e9–1e13
solar masses agree with the product derivative to 3.67e-9 (q5) and 1.74e-8
(q10), using a logarithmic step of 1e-4.

The initial fixed-cell-only experiment at commit `6939b32` was insufficient:
integral errors were around 1.4e-7 and derivative errors reached 7.83e-5 when
cells crossed a spectrum interpolation knot. Its original evidence is retained
in `validation/sharp-k-wmap7-20260910.json`; the final independent results are
in `validation/sharp-k-wmap7-knots-20260910.json`. This was a refinement of the
numerical integration, with unchanged power, mass mapping and normalization.

`integration_breakpoints` is an optional one-dimensional array of finite,
positive, strictly increasing wavenumbers. Integrators may use it to partition
the supplied power; smooth custom spectra need not expose it. Tabulated power
returns a copy of its knots, and transfer-modified power forwards the base
knots. The exact breakpoint content and partition algorithm enter the variance
identifier, so even the intermediate v3 result cannot be reused accidentally.

These establish the shared integrator's numerical behavior, not full WDM
scientific acceptance. The downstream W change must remove its variance
projection, update the numerical identity, and separately validate the whole
catalog and observables against the independent corrected reference. The
q5/q10 physical choices and their normalization are unchanged here.
