# Explicit odeint controls in shared evolution

`solve_evolution(..., method="odeint", odeint_options={...})` preserves the
selected LSODA step limits, orders, critical times and dense/banded Jacobian.
The Jacobian uses the same time-first argument convention as the RHS, including
`args`. Returned values and error handling remain controlled by ITAMAE.

The allowed keys are Dfun, col_deriv, ml, mu, tcrit, h0, hmax, hmin, ixpr,
mxstep, mxhnil, mxordn, mxords and printmessg. Pass rtol/atol and args through
the explicit solve_evolution arguments. tfirst and full_output are internal
controller choices. Unknown, duplicate, non-finite or malformed options fail
before evaluating the physical RHS. Options for odeint cannot be passed to a
different solver. printmessg reports a successful diagnostic as a warning;
failed integration or other ODEPACK warnings raise RuntimeError.

Twenty-one initial protection cases failed before the API extension. Direct
SciPy comparisons in both time directions, dense and banded stiff Jacobians,
invalid Jacobians and explicit failure limits now pass. The direct odeint
comparisons are bitwise equal in the same environment. Default solver choices,
ODE tolerances and physical equations are unchanged.
