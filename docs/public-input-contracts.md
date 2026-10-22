# Direct catalog, unit and native-background input contracts

Population execution validated its stage arrays, but standalone catalog
construction and NPZ loading could bypass those checks. Complex weights and
unit inputs could lose their imaginary part on conversion; non-finite physical
columns and overflowing products/reductions could be exposed as results.
NativeFlatLCDM also accepted invalid parameters/redshifts and rounded its
identifier enough to collide for distinct floating cosmologies.

## Public boundary

- Physical catalog columns accept finite real numeric or boolean arrays, with
  the existing aligned scalar/array shape contract. Object, string and complex
  columns are rejected. Signs and physical meaning still belong to the model.
- Weight factors remain finite and nonnegative. Their product is evaluated in
  insertion order and must be representable without overflow. Underflow retains
  NumPy's existing convention; no clipping or renormalization is added.
- Weighted sums require finite real aligned values and finite products/sums.
  Histograms require finite counts. An empty count histogram is zero; a density
  normalized by zero total weight raises instead of returning NaN.
- NPZ loading validates array values and rejects non-mapping manifests and
  duplicate column/weight names. Saving rechecks arrays and weights. The file
  schema and original metadata identities are unchanged.
- Native/Astropy unit conversion rejects complex or non-finite inputs before
  casting. Unit conversion does not impose model-specific sign constraints.
- NativeFlatLCDM requires finite scalar 0 < omega_m0 <= 1 and h > 0. Its
  matter-plus-nonnegative-Lambda domain excludes empty matter and negative
  Lambda. Real finite z > -1 is supported, including future redshifts needed
  by centered differentiation at zero. Nonrepresentable expansion rates and
  densities raise explicitly. The identifier keeps float round-trip precision.

## Evidence

The initial new independent boundary file produced 29 failures, three passes
and four warnings before the fix. Its repaired cases include invalid saved
columns, duplicate manifests, factor/reduction overflow, complex units and
cosmology input rejection. Einstein-de Sitter growth, density-fraction and
analytic age identities include z=-0.5,0,1,3. A further explicit density-overflow
case covers otherwise finite large redshift.

All 166 core tests pass with the change, with Ruff and mypy checked separately.
C (53), SI (16 standard + 20 equation tests) and F (35) passed the cross-variant run.
W exposed an existing concentration-grid branch which evaluated NaN formation
redshifts before dropping those candidates. The W consumer must exclude those
same non-real candidates before calling the strict cosmology backend; its full
q5/q10 reference agreement is required before this core change merges.

This boundary change does not change physical defaults, insert a survival mask,
repair an invalid model silently, or establish scientific calibration.

The W consumer fix passed all 37 tests with this strict core, including both
independent q5/q10 full-catalog references at their unchanged tolerances. It
removes the same non-real trial candidates before growth evaluation. No
concentration/solver/count default was changed.

A Hypothesis draw also exposed overflow in the NFW inverse upper-bracket
residual for a subnormal target. The residual sign is now represented by the
largest finite float only when the original division would overflow; normal
residuals and the root definition are unchanged. An independent 1e-310 target
uses the exact small-radius asymptote and treats RuntimeWarnings as errors.
