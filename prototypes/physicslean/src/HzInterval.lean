import HzCertificate
import FloatLib.Tactic.Interval
import Mathlib.Tactic.NormNum.RealSqrt

/-! Two uses of FloatLib: verify proposed decimal endpoints directly on the real
expression, or compute dyadic endpoints at a specified precision using its proved
rational square-root enclosure. -/

namespace HzInterval

open FloatLib.Numerics

/-- Generate a dyadic enclosure on a grid with `bits` fractional bits.
Invalid physical inputs do not produce an interval. -/
def enclose (h₀ Ω z : ℚ) (bits : Nat) : Option (Interval ℚ) :=
  if 0 < h₀ ∧ 0 ≤ Ω ∧ Ω ≤ 1 ∧ 0 ≤ z then
    some (Interval.sqrtPointBounds (h₀ ^ 2 * ((1 - Ω) + Ω * (1 + z) ^ 3)) bits)
  else none

/-- The returned endpoints enclose Physlib's real prediction, and their
difference is at most one dyadic grid step. -/
theorem enclose_sound {h₀ Ω z : ℚ} {bits : Nat} {I : Interval ℚ}
    (h : enclose h₀ Ω z bits = some I) :
    ((I.lo : ℝ) ≤ HzCertificate.modelH h₀ Ω z ∧ HzCertificate.modelH h₀ Ω z ≤ I.hi) ∧
      I.hi - I.lo ≤ 1 / (2 : ℚ) ^ bits := by
  unfold enclose at h
  split at h
  · rename_i valid
    obtain ⟨hH, hΩ, hΛ, hz⟩ := valid
    cases h
    have hr : 0 ≤ (1 - Ω) + Ω * (1 + z) ^ 3 :=
      add_nonneg (sub_nonneg.mpr hΛ) (mul_nonneg hΩ (pow_nonneg (add_nonneg zero_le_one hz) _))
    have hr' : 0 ≤ (1 - (Ω : ℝ)) + Ω * (1 + (z : ℝ)) ^ 3 := by exact_mod_cast hr
    have hz' : 0 ≤ (z : ℝ) := by exact_mod_cast hz
    have hH' : 0 ≤ (h₀ : ℝ) := by exact_mod_cast hH.le
    have hs : ((h₀ ^ 2 * ((1 - Ω) + Ω * (1 + z) ^ 3) : ℚ) : ℝ) =
        HzCertificate.modelH h₀ Ω z ^ 2 := by
      push_cast
      exact (HzCertificate.modelH_sq _ _ _ hz' hr').symm
    have hroot := Interval.containsReal_sqrtPointBounds _ bits (mul_nonneg (sq_nonneg h₀) hr)
    rw [hs, Real.sqrt_sq (HzCertificate.modelH_nonneg _ _ _ hH' hz')] at hroot
    exact ⟨(Interval.containsReal_some_iff _ _).mp hroot, Interval.sqrtPointBounds_width_le _ _⟩
  · contradiction

/-- Check rational endpoints using the real formula for the Physlib definition.
The RealSqrt normalization extension lets the tactic keep exact endpoints
when the square root and the requested bound are rational. -/
macro "hz_interval" ppSpace "(" "precision" " := " bits:num ")" : tactic =>
  `(tactic|
    (rw [HzCertificate.modelH_eq _ _ _ (by norm_num)]
     constructor <;> interval (precision := $bits)))

macro "hz_interval" : tactic =>
  `(tactic| hz_interval (precision := 64))

/-- The decimal bounds proposed for the fiducial example follow directly from
the Physlib formula by interval arithmetic at 64 fractional bits. -/
theorem certified_interval :
    12324771803 / 100000000 ≤ HzCertificate.modelH 70 (3 / 10) 1 ∧
      HzCertificate.modelH 70 (3 / 10) 1 ≤ 12324771804 / 100000000 := by
  unfold HzCertificate.modelH Cosmology.FLRW.FriedmannEquation.reducedHubble
  constructor <;> interval

#print axioms enclose_sound
#print axioms certified_interval

end HzInterval
