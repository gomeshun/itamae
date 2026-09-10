# Stable moving-boundary sharp-k integration

Sharp-k variance is the integral of `k^3 P(k)/(2 pi^2)` in log-wavenumber up to
the mass-dependent physical cutoff, truncated at the configured finite power
domain. The previous integrator stretched the entire Simpson grid whenever the
upper cutoff changed. This moved samples across the tabulated spectrum's
features and caused fluctuations even where physical variance is nearly flat.

Specification `integrated-variance:v3` uses `n_k-1` fixed logarithmic cells and
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
partition equivalence, and 257/513-cell-node convergence. The complete core
suite has 105 passing cases, with Ruff/mypy passing and branch coverage 80%.

The WMAP7 table was additionally integrated by adaptive QUADPACK on its own
3,001 original knots, independently of the product's fixed log-k cells. For the
explicit q5 and q10 spectra at particle mass 2 keV, 4,097 product grid nodes
give a maximum relative integral discrepancy of 1.41e-7 over 13 representative
masses. Increasing to 8,193 nodes changes the independently normalized variance
by less than 9.52e-8. The former grid produced positive variance steps of up to
8.66e-4 across a 1,000-mass scan; the fixed-cell scan has no positive steps.

These establish the shared integrator's numerical behavior, not full WDM
scientific acceptance. The downstream W change must remove its variance
projection, update the numerical identity, and separately validate the whole
catalog and observables against the independent corrected reference. The
q5/q10 physical choices and their normalization are unchanged here.
