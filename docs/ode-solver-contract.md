# Explicit ODE solver contract

The common evolution entry point now accepts `method="odeint"` so a variant
can preserve its ODEPACK LSODA calculation, rather than silently changing to
RK45 during a structural migration. No SASHIMI law or constants enter ITAMAE.
`rhs(t, y, *args)` and the time-first result shape are shared by both paths.

Unspecified tolerances preserve the previous ITAMAE RK45 values (1e-8 and
1e-10). The explicitly selected odeint path passes None through to SciPy,
which preserves the historical SciPy defaults (approximately 1.49012e-8 for
both). These differ from solve_ivp's LSODA wrapper defaults, so selecting
`method="LSODA"` is not claimed to reproduce odeint.
See [SciPy's odeint contract](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html).

Finite real initial states, monotonic finite grids, positive tolerances and
matching derivative shapes are required. Solver failure, incomplete output
and non-finite callback/output values are explicit errors. An ODEintWarning
is promoted to RuntimeError with its cause, not ignored.

Validation on 2026-09-10: all 18 new contract cases failed before this change.
Exponential evolution in both time directions matches the analytical result
within 1e-7 and the direct odeint call bit for bit under identical dependencies.
The existing four evolution tests still pass. This checks the numerical
adapter; it does not validate a variant's tidal-stripping prescription.
