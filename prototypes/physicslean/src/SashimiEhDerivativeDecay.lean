import SashimiEhDerivativeBounds

/-! Global and high-frequency envelopes for the logarithmic-mass window
derivative. Both bounds are independent of the power spectrum. -/
namespace SashimiEhDerivativeBounds
open Real SashimiSpectrum

private theorem abs_shape {x : ℝ} :
    |shape x| ≤ 2 * |topHat x| * (|sin x / x| + |topHat x|) := by
  simp only [shape, abs_mul, abs_of_nonneg (by norm_num : (0:ℝ) ≤ 2), ← mul_assoc]
  exact mul_le_mul_of_nonneg_left (abs_sub _ _) (by positivity)

theorem abs_shape_le_four {x : ℝ} (hx : 0 < x) : |shape x| ≤ 4 := by
  have hw : |topHat x| ≤ 1 := abs_le.mpr (topHat_bounds x hx.le)
  have hs : |sin x / x| ≤ 1 := by
    rw [abs_div,abs_of_pos hx,div_le_one hx]
    simpa [abs_of_pos hx] using (Real.abs_sin_le_abs : |sin x| ≤ |x|)
  calc
    |shape x| ≤ _ := abs_shape
    _ ≤ 2*1*(1+1) := by gcongr
    _ = 4 := by norm_num

theorem abs_shape_decay {x : ℝ} (hx : 1 ≤ x) : |shape x| ≤ 84 / x^3 := by
  have hp : 0 < x := zero_lt_one.trans_le hx
  have hw : |topHat x| ≤ 6 / x^2 := (topHat_decay x hp).trans (by
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith only [mul_nonneg (sq_nonneg x) (sub_nonneg.mpr hx)])
  have hww : |topHat x| ≤ 6 / x := hw.trans (by
    apply div_le_div_of_nonneg_left (by norm_num) hp
    nlinarith only [mul_nonneg hp.le (sub_nonneg.mpr hx)])
  have hs : |sin x / x| ≤ 1 / x := by
    rw [abs_div, abs_of_pos hp]
    exact div_le_div_of_nonneg_right (Real.abs_sin_le_one x) hp.le
  calc
    |shape x| ≤ _ := abs_shape
    _ ≤ 2*(6/x^2)*(1/x + 6/x) :=
      mul_le_mul (mul_le_mul_of_nonneg_left hw (by norm_num)) (add_le_add hs hww) (by positivity) (by positivity)
    _ = 84 / x^3 := by field_simp; ring

#print axioms abs_shape_decay
end SashimiEhDerivativeBounds
