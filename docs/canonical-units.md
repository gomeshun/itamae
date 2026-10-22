# Canonical units and reproducibility metadata

ITAMAE performs large internal calculations with plain NumPy arrays. Public
unit backends validate inputs, convert them to this canonical representation,
and optionally convert outputs to user-selected units.

| Physical type | Canonical unit |
| --- | --- |
| dimensionless | dimensionless |
| mass | `Msun` |
| length | `Mpc` |
| velocity | `km / s` |
| time | `Gyr` |
| density | `Msun / Mpc3` |
| cross section per mass | `cm2 / g` |

The schema version is exposed as
`itamae.units.CANONICAL_UNIT_SCHEMA_VERSION`. Plain values passed through
`NativeUnits` are interpreted in these units. `AstropyUnits` requires
dimensionally compatible `Quantity` inputs and can convert canonical outputs to
equivalent display units.

Every persisted catalog must record:

- `schema_version`;
- `model_identifier`, owned by the SASHIMI variant;
- `backend_identifier`, supplied by `BackendConfig`;
- an optional `source_identifier` such as a git commit or golden-fixture ID.

Backend and model identifiers must be stable enough for cache keys and
regression provenance. Changing cosmology parameters, unit schema, a physical
prescription, or a numerical convention requires a different identifier.

Power-spectrum and variance caches additionally record or hash:

- power-table or source identifier;
- model-supplied transfer/power-ratio identifier;
- smoothing-window identifier and mass-assignment coefficient;
- integration bounds and resolution;
- cosmology/backend identifier;
- exact mass grid;
- cache-schema version.

ITAMAE never infers whether a variant transfer function is an amplitude or a
power ratio. The `TransferModifiedPowerSpectrum` API accepts the latter
explicitly, preventing an unrecorded factor-of-two change in its exponent.
