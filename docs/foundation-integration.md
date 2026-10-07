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
- log-log interpolation of explicitly identified tabulated power spectra;
- model-supplied power-ratio composition without assuming whether a transfer
  amplitude should be squared;
- spherical top-hat and sharp-k windows;
- chunked variance integration with explicit mass-filter calibration;
- content-addressed variance cache keys and validated atomic cache files;
- spherical-orbit turning points, radial periods, radial PDFs, and normalized
  shell measures in a model-supplied potential.

These are mechanisms rather than physical defaults and therefore belong in
ITAMAE.

## Components deliberately deferred

- WDM and FDM transfer functions remain owned by SASHIMI-W and SASHIMI-F.
- Spectrum normalization, physical/comoving units, factors of ``h``, filter
  defaults, and cache locations remain explicit variant configuration.
- Orbit-infall distributions, dynamical friction, disruption prescriptions,
  and phase-space transport remain model-owned Phase 8 work. The integrated
  orbit kernel records only normalized radial shell probabilities.

The older branch remains useful as a prototype, but it must not be merged
wholesale because its flat module layout, backend API, and FDM transfer formula
do not match the canonical contracts.
