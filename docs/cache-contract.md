# Variance cache integrity

Spectrum identifiers always contain a digest of the canonical little-endian
float64 wavenumber and power arrays, their shapes, interpolation and extrapolation
policies, and any human-readable source label. Reusing a label for changed data
does not reuse its numerical identity. The integrated-variance identity also
includes the smooth-window derivative step. Batch size is deliberately excluded:
it changes allocation, not the calculation.

Variance cache schema 2.0 binds the request key, mass array and variance array
with a SHA-256 payload digest. Both a matching request key and matching payload
digest are required when loading. Nonfinite, negative, misaligned or undersized
tables are rejected. Schema 1.0 files do not provide this check and must be
recomputed; they are never silently promoted. Atomic writes are retained.

This is accidental-corruption detection, not authentication of untrusted files.
The caller still supplies the cosmology/backend, model and settings to the
request-key builder and is responsible for correct identifiers for custom
callable spectra and transfer functions.

## Validation, 2026-09-10

Before this fix all eight contract cases failed: repeated labels aliased changed
tables; extrapolation and derivative steps did not affect identity; changed
finite mass/variance values were accepted; empty and one-node saved grids were
accepted; and the old unverified schema was accepted. After the fix all eight
pass. The full core suite has 88 passing tests. This change does not modify any
physical formula, quadrature, or spectrum values. It invalidates old identities
and cache files so that numerical inputs cannot silently reuse stale results.
