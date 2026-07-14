# Foundation integration decisions

This document records the consolidation of the two initial implementation
branches into the canonical ITAMAE foundation.

## Canonical branch

`agent/initial-roadmap-implementation` was selected as the implementation base
because its package layout and public imports match the adapters developed in
SASHIMI-C, SASHIMI-SI, SASHIMI-W, and SASHIMI-F. It also provides the clearer
separation among backends, state types, numerical mechanisms, and physical
model protocols.

## Components retained from `agent/implement-roadmap`

The following variant-independent mechanisms were reimplemented with
validation, documentation, and tests:

- logarithmic and redshift grid helpers;
- a vectorized Picard perturbative evolution solver with optional Shanks
  acceleration;
- an accretion-batch broadcasting helper;
- weighted histograms, explicit-generator Poisson realizations, and catalog
  concatenation.

These are mechanisms rather than physical defaults and therefore belong in
ITAMAE.

## Components deliberately deferred

- WDM and FDM transfer functions remain owned by SASHIMI-W and SASHIMI-F.
- Tabulated power spectra, smoothing windows, variance integration, derivative
  evaluation, and cache keys will be designed in roadmap Phase 6 using W/F
  golden fixtures.
- Turning points, radial periods, and phase-space kernels remain Phase 8 work.
  They require an explicit spatial representation and normalization contract.

The older branch remains useful as a prototype, but it must not be merged
wholesale because its flat module layout, backend API, and FDM transfer formula
do not match the canonical contracts.
