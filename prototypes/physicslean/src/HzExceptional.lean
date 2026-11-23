import HzRoundoff

/-! # Extreme valid inputs of the binary64 Hubble program -/

namespace HzExceptional

open HzRoundoff
open FloatLib.Floats FloatLib.Floats.Formats.BinaryInterchange

/-- In all probes below the exact Physlib result is simply the positive input H₀. -/
theorem lambda_only (h₀ z : ℝ) : HzCertificate.modelH h₀ 0 z = h₀ := by
  simp [HzCertificate.modelH_eq_unconditional]

/-- Input overflow produces the IEEE positive-infinity encoding. -/
theorem input_overflow :
    (evaluate (2 ^ 1024) 0 0).toNatBits = 0x7ff0000000000000 := by
  float_decide [evaluate, hubble64]

/-- RN-even sends half the least subnormal to +0 and preserves the least subnormal. -/
theorem input_underflow :
    (evaluate (1 / 2 ^ 1075) 0 0).toNatBits = 0 ∧
    (evaluate (1 / 2 ^ 1074) 0 0).toNatBits = 1 := by
  constructor <;> float_decide [evaluate, hubble64]

/-- All three inputs are representable, but cubing the scale factor overflows.
The subsequent zero times infinity propagates a NaN through the square root. -/
theorem intermediate_nan :
    (ExecFloat.Binary.ofRat (2 ^ 400) : Binary64).toRat? = some (2 ^ 400) ∧
    (evaluate 1 0 (2 ^ 400)).isNaN = true := by
  constructor <;> float_decide [evaluate, hubble64]

/-- With z = 2^1024 the redshift conversion already produces infinity. -/
theorem nonfinite_redshift_nan :
    (evaluate 1 0 (2 ^ 1024)).isNaN = true := by
  float_decide [evaluate, hubble64]

#print axioms lambda_only
#print axioms input_overflow
#print axioms input_underflow
#print axioms intermediate_nan
#print axioms nonfinite_redshift_nan

end HzExceptional
