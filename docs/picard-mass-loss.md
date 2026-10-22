# Shared Picard log-mass integration

`PicardMassLossSolver` solves

`d ln(m)/dz = rate(z) * exp(exponent(z) * ln(m/reference_mass(z)))`.

The caller supplies all three coefficient functions, compatible mass units,
a state identity, and any scientifically validated domain. ITAMAE contains no
SASHIMI coefficients, cosmological calibration, WDM transfer function or
SIDM profile evolution. C and W use endpoint tables; SI requests the complete
mass history needed by its gravothermal model. Both paths use the same
nonlinear Picard update. The existing `PerturbativeEvolutionSolver` and
`solve_evolution` APIs remain available.

```python
from itamae.evolution import PicardDomain, PicardMassLossSolver

solver = PicardMassLossSolver(
    rate=model.rate,
    exponent=model.exponent,
    reference_mass=model.reference_mass,
    z_min=0.0,
    z_max=7.0,
    state_key=model.numerical_state_key,
    domain=PicardDomain((0.0, 7.0), (-24.0, 3.0)),
)
final_mass = solver.endpoint(mass_at_accretion, z_acc, 0.0)
history = solver.history(mass_at_accretion, z_acc, redshifts)
reference = solver.direct(mass_at_accretion, z_acc, redshifts)
diagnostics = solver.diagnostics()
```

The bounds in this example are supplied by the model, not imposed by ITAMAE.
Omitting `domain` disables the domain-based fallback. Output histories have
shape `(n_redshifts, *mass_shape)` and preserve repeated, non-increasing output
redshifts. Integration always starts at the declared accretion redshift even
if it is absent from the output grid. A scalar target omits the time axis.

Endpoint defaults are 48 accretion-redshift nodes, 32 log-mass-ratio nodes,
129 integration points, cubic interpolation and three nonlinear updates.
History defaults are 100 integration points and three updates; large batches
use at least 32 mass nodes, with no more than 0.25 dex between mass nodes.
Endpoint interpolation changes to linear for a two- or three-node axis.
Explicit `PicardMassLossTable` queries reject coordinates outside their table;
`endpoint()` expands a requested table to contain in-domain queries.

A fourth update provides a convergence estimate. Points exceeding the default
`8e-4` iteration tolerance, or outside the caller's domain, are recomputed by
DOP853 integration of log mass through ITAMAE's existing ODE controller
(`rtol=1e-11`, `atol=1e-12`). Every fallback emits `PicardFallbackWarning` and
is counted with its reason and mass count in `diagnostics()`. The estimate
covers iterations only: **it is not a total interpolation or quadrature error
bound**. Validate those errors separately using grid refinement and a direct
reference for the intended physical model and input range.

`state_key()` must return an immutable equality-comparable value covering
all mutable inputs used by the coefficient callbacks. Changing that value
invalidates coefficient/host/table caches and rejects retained stale tables.
Diagnostics describe the current state; counters reset when the state changes.
No process-local callback identity is used as a portable scientific provenance
identifier: save model inputs, solver settings and source revisions as well.

The numerical implementation is derived from the standalone SASHIMI-C helper
at `2ffefc8f02ff1db2f634aae73c1b2a394914a98e` (MIT notice retained), including
its endpoint/history algorithms. The model adapter, domain policy, common ODE
controller, explicit callback contract and diagnostics are separated here so
all variants call one implementation.
