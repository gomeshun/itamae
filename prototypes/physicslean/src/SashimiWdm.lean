import SashimiSpectrumVariance

/-! Published fit conventions and their physical interpretation on a specified
calibration band. A simulation, observation or Boltzmann table enters as data
or as an explicit hypothesis. -/

namespace SashimiSpectrum
open Real Set MeasureTheory

theorem amplitude_mem_unit {α k : ℝ} (hα : 0 ≤ α) (hk : 0 ≤ k) :
    0 < wdmAmplitude α k ∧ wdmAmplitude α k ≤ 1 := by
  unfold wdmAmplitude
  have h : 0 ≤ (α*k) ^ (2*(1.12 : ℝ)) := rpow_nonneg (mul_nonneg hα hk) _
  exact ⟨rpow_pos_of_pos (by linarith) _, rpow_le_one_of_one_le_of_nonpos (by linarith) (by norm_num)⟩

/-- Holding the reference power and alpha fixed, q10 suppresses at least as
much power as q5. Independent sigma8 recalibrations need independent anchors. -/
theorem q10_le_q5 {α k : ℝ} (hα : 0 ≤ α) (hk : 0 ≤ k) :
    wdmPowerRatio .standardT2Q10 α k ≤ wdmPowerRatio .publishedQ5 α k := by
  rw [standard_ratio_eq_amplitude_sq α k hα hk, published_ratio_eq_amplitude]
  obtain ⟨hp, h⟩ := amplitude_mem_unit hα hk
  nlinarith [mul_nonneg hp.le (sub_nonneg.mpr h)]

/-- Empirical error in physical **power**, with a supplied reference spectrum,
convention, domain in Mpc⁻¹ and absolute error in Mpc³. -/
def WdmFitAgreement (mode : WdmConvention) (α : ℝ) (domain : Set ℝ)
    (reference physical error : ℝ → ℝ) : Prop :=
  ∀ k ∈ domain, |physical k - modifiedPower reference (wdmPowerRatio mode α) k| ≤ error k

theorem physical_power_bounds {mode α domain ref physical error}
    (h : WdmFitAgreement mode α domain ref physical error) {k L U : ℝ} (hk : k ∈ domain)
    (hfit : L ≤ modifiedPower ref (wdmPowerRatio mode α) k ∧
      modifiedPower ref (wdmPowerRatio mode α) k ≤ U) :
    L - error k ≤ physical k ∧ physical k ≤ U + error k := by
  obtain ⟨hl, hu⟩ := abs_le.mp (h k hk)
  constructor <;> linarith [hfit.1, hfit.2]

/-- The order survives any common nonnegative integration weight on a finite band. -/
theorem integrated_q10_le_q5 {α a b : ℝ} (hα : 0 ≤ α) (ha : 0 ≤ a) (hab : a ≤ b)
    (weight : ℝ → ℝ) (hw : ∀ k ∈ Icc a b, 0 ≤ weight k)
    (h10 : IntervalIntegrable (fun k => weight k * wdmPowerRatio .standardT2Q10 α k) volume a b)
    (h5 : IntervalIntegrable (fun k => weight k * wdmPowerRatio .publishedQ5 α k) volume a b) :
    (∫ k in a..b, weight k * wdmPowerRatio .standardT2Q10 α k) ≤
      ∫ k in a..b, weight k * wdmPowerRatio .publishedQ5 α k := by
  exact intervalIntegral.integral_mono_on hab h10 h5 fun k hk =>
    mul_le_mul_of_nonneg_left (q10_le_q5 hα (ha.trans hk.1)) (hw k hk)

#print axioms integrated_q10_le_q5
#print axioms physical_power_bounds
end SashimiSpectrum
