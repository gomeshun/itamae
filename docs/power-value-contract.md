# Real power values and aligned transfers

Tabulated power inputs and queries now reject complex, string, boolean and
other nonnumeric arrays before a float conversion can discard information.
Transfer callbacks must return a finite nonnegative scalar or the query shape;
the base spectrum must preserve that shape. Product overflow fails explicitly
without returning infinity or leaking a floating-point warning.

This completes the same real-valued boundary contract already used by the
variance adapters and windows. No supported real-valued spectrum, interpolation,
identifier, unit convention or physical transfer changes. A scalar power ratio
still broadcasts over any query shape.

Fourteen new tests failed before the fix; all 309 core tests now pass, including
existing analytic spectrum/variance tests and optional backends. The package
retains the variant-owned distinction between transfer amplitude and power ratio.
