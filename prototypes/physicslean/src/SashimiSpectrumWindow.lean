import SashimiSpectrum
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-! Stable bounds and the smooth mass derivative of the spherical window. -/
namespace SashimiSpectrum
open Real Set Filter MeasureTheory
open scoped Topology
noncomputable section

private theorem windowNumerator (x : ℝ) :
    (∫ t in 0..x, t * sin t) = sin x - x * cos x := by
  have hd (t : ℝ) : HasDerivAt (fun y => sin y - y*cos y) (t*sin t) t := by
    convert (Real.hasDerivAt_sin t).sub ((hasDerivAt_id' t).mul (Real.hasDerivAt_cos t)) using 1
    ring
  simpa using intervalIntegral.integral_eq_sub_of_hasDerivAt (a := 0) (b := x)
    (fun t _ => hd t) ((continuous_id.mul Real.continuous_sin).intervalIntegrable _ _)

/-- Unit bounds remain stable across many decades in kR. -/
theorem topHat_bounds (x : ℝ) (hx : 0 ≤ x) : -1 ≤ topHat x ∧ topHat x ≤ 1 := by
  rcases hx.eq_or_lt with h | hx
  · simp [← h, topHat]
  have hInt : IntervalIntegrable (fun t => t * sin t) volume 0 x :=
    (continuous_id.mul Real.continuous_sin).intervalIntegrable _ _
  have hp := (continuous_id.pow 2).intervalIntegrable (μ := volume) 0 x
  have hlo := intervalIntegral.integral_mono_on hx.le hp.neg hInt (fun t ht => by
    have h := (abs_le.mp ((Real.abs_sin_le_abs).trans_eq (abs_of_nonneg ht.1))).1
    dsimp; nlinarith [mul_le_mul_of_nonneg_left h ht.1])
  have hhi := intervalIntegral.integral_mono_on hx.le hInt hp (fun t ht => by
    dsimp; nlinarith [mul_le_mul_of_nonneg_left (Real.sin_le ht.1) ht.1])
  simp only [windowNumerator, Pi.neg_apply, intervalIntegral.integral_neg] at hlo hhi
  unfold topHat
  rw [ite_eq_right hx.ne']
  constructor
  · apply (le_div_iff₀ (by positivity)).2
    norm_num at hlo; linarith
  · apply (div_le_iff₀ (by positivity)).2
    norm_num at hhi; linarith

/-- One-sided second-order bound for the small-argument branch. -/
theorem topHat_taylor_lower (x : ℝ) (hx : 0 ≤ x) : 1 - x^2/10 ≤ topHat x := by
  rcases hx.eq_or_lt with h | hx
  · simp [← h, topHat]
  have h := intervalIntegral.integral_mono_on (μ := volume) (f := fun t : ℝ => t^2-t^4/6)
    (g := fun t => t*sin t) hx.le ((by fun_prop : Continuous (fun t : ℝ => t^2-t^4/6)).intervalIntegrable _ _)
    ((continuous_id.mul Real.continuous_sin).intervalIntegrable _ _) (fun t ht => by
      nlinarith [mul_le_mul_of_nonneg_left (Real.sin_ge_sub_cube ht.1) ht.1])
  rw [intervalIntegral.integral_sub (f := fun t : ℝ => t^2) (g := fun t : ℝ => t^4/6) ((continuous_id.pow 2).intervalIntegrable _ _)
      (((continuous_id.pow 4).div_const 6).intervalIntegrable _ _),
    intervalIntegral.integral_div, integral_pow, integral_pow, windowNumerator] at h
  rw [topHat, ite_eq_right hx.ne']
  apply (le_div_iff₀ (by positivity)).2
  norm_num at h; nlinarith

theorem topHat_decay (x : ℝ) (hx : 0 < x) : |topHat x| ≤ 3 * (1+x)/x^3 := by
  have h := abs_sub (sin x) (x * cos x)
  rw [abs_mul, abs_of_pos hx] at h
  have h' := mul_le_mul_of_nonneg_left (Real.abs_cos_le_one x) hx.le
  simp [topHat, hx.ne', abs_div, abs_mul, abs_of_pos hx]
  gcongr
  linarith [Real.abs_sin_le_one x]



end
end SashimiSpectrum
