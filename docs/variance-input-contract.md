# Public variance input and output contract

The callable adapter previously accepted zero/negative masses, invalid redshifts,
negative sigma, non-finite results and unrelated output shapes. Some NumPy complex
arrays were converted to real values with their imaginary part discarded. The
integrated adapter checked intermediate power values, but not all callback types,
window shapes or the final growth-factor multiplication.

Both adapters now require real finite numeric coordinates and positive mass.
Mass and redshift follow NumPy broadcasting. A redshift domain beyond finiteness
belongs to the selected cosmology/growth prescription; the generic adapter does
not impose one. Inputs are checked before invoking supplied callbacks. The
callable adapter preserves the original argument shapes passed to those callbacks.

Sigma is finite and nonnegative; the variance derivative is finite and may have
either sign. No clipping or sign correction is performed. Results must match the
broadcast coordinate shape, with a scalar callback result explicitly permitted
and expanded to that shape. Empty coordinates remain supported. Identifiers must
be nonempty strings and both callbacks must be callable.

Power and window callbacks must return exactly the requested shape. Power is real,
finite and nonnegative; a window may have either sign. A growth callback may return
a scalar or the requested shape and must be real, finite and nonnegative. Complex
integration knots and nonscalar/nonreal configuration values are rejected before
conversion. Overflow, division by zero and invalid arithmetic within public
variance evaluations raise `ValueError` with the underlying exception attached.
Underflow retains NumPy's normal behavior; this does not promise representability
for arbitrary extreme inputs.

The formulas, integration grids, operation grouping, growth normalization,
derivative sign and numerical identifier are unchanged for valid calculations.
This is an input contract change, not evidence that a particular physical
variance prescription is calibrated.

Validation: `tests/test_variance_input_contract.py` initially had 44 failures and
13 passes (including explicit failures for complex scalars that previously raised
an inconsistent exception type). The complete 57-case file passes after the fix.
It covers callback corruption, extra output dimensions, overflow, pre-callback
validation, scalar broadcasting, empty inputs and unchanged analytic results.
The existing independent sharp-k integration/derivative tests retain their
original tolerances.
