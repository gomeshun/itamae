# Standard calculation provenance

`build_calculation_metadata` creates provenance schema
`itamae:calculation:v2`. The variant supplies a versioned
`calculation_specification` plus the explicit physics choices, power/variance,
solver, backend, canonical units and source identifiers. New products do not
accept `physics_mode`; that key is rejected even when injected through `extra`.
The source-resolution and reserved-field checks are shared with the historical
builder.

`build_migration_metadata` retains its previous output vocabulary. Existing
catalog archives with legacy/consistent mode labels round-trip unchanged and
are never relabeled as a current calculation. Catalog column schema and unit
schema remain 1.0: this change versions provenance, not column meanings.

Eight tests protect coexistence and serialization of old/new records, missing
specifications, and conflicting identity fields. The full suite passes 114
cases; Ruff, formatting and mypy also pass.
