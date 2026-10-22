# NFW numerical boundaries

The enclosed-mass function is `f(x)=ln(1+x)-x/(1+x)`. Direct subtraction loses
relative precision at small positive radii. For `x<0.01` the implementation now
uses `sum[(-1)^n (n-1) x^n/n, n=2..10]`; the first omitted term is below floating
precision relative to the leading `x^2/2` at this switch. Larger radii keep the
original log1p expression.

Inversion is bracketed in log-radius. An absolute root tolerance in radius must
not map a positive, representable small radius to zero. The bracket follows
`f(x)<=x^2/2` and `f(exp(y+1))>=y`. Zero is exact; finite nonnegative inputs are
required; a requested mass-function value above the maximum representable
radius's value is an explicit error. Subnormal underflow remains a limit of the
floating-point representation, not an alternative physical prescription.

Profile scale radius/density must be finite and positive. Radii are finite and
nonnegative; density requires a strictly positive radius. Central potential uses
the analytic limit without evaluating zero over zero in a discarded branch.

## Evidence, 2026-09-10

The protection suite initially had 13 failures and two passes. The updated
implementation passes all 15 cases, including comparison to 340-digit standard
library Decimal arithmetic at radii down to `1e-150`, independent inverse
recovery, invalid parameters and warning-free central potential. The complete
core suite has 103 passes, with 80% branch coverage; Ruff and mypy pass.

This fixes evaluation of the existing NFW formula; it introduces no profile or
survival prescription. C's independent 50-digit Lambert-W reference separately
tests the physical inverse. Cross-package golden tests remain unchanged.

Profile parameters, radii and enclosed-mass coordinates require real numeric
inputs. Complex values (including zero imaginary parts), strings and booleans
are rejected before float conversion; no imaginary component is discarded.
Smoothing-window arguments follow the same contract. All 34 new non-real input
cases failed before the guard and pass afterward. Physical float arithmetic is
unchanged; existing high-precision inverse and small-radius checks remain.
