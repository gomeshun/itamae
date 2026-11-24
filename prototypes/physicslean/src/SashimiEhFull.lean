import SashimiEisensteinIntegral
import SashimiVarianceTails

/-! Turn a certified continuous EH band integral into an improper-integral
enclosure using transfer envelopes outside the band. -/
namespace SashimiEhFull
open MeasureTheory Set SashimiEisensteinIntegral SashimiEisensteinHu
open SashimiVarianceTails FloatLib.Numerics

noncomputable def spectrumWeight (c : Cosmology) (R : ℝ) :=
  weight (transferReal (coefficients c)) (c.ns:ℝ) R

noncomputable def fullRaw (c : Cosmology) (R : ℝ) : ℝ :=
  ∫ k in Ioi (0:ℝ), spectrumWeight c R k

theorem measurableWeight (c : Cosmology) (R : ℝ) : Measurable (spectrumWeight c R) := by
  have hT : Measurable (transferReal (coefficients c)) := by
    change Measurable (fun k => variableTransfer.eval (env k (coefficients c)))
    apply SashimiAnalyticWindow.measurableFormula
    intro i
    cases i with
    | zero => exact measurable_id
    | succ i => exact measurable_const
  have hw : Measurable SashimiSpectrum.topHat := by
    unfold SashimiSpectrum.topHat
    exact Measurable.ite (measurableSet_singleton 0) measurable_const (by fun_prop)
  unfold spectrumWeight weight
  fun_prop

theorem weight_eq (c : Cosmology) (R : ℝ) {k : ℝ} (hk : 0 < k) :
    spectrumWeight c R k = Real.exp (((2+c.ns : ℚ):ℝ) * Real.log k) *
      (transferReal (coefficients c) k)^2 * SashimiSpectrum.topHat (k*R)^2 := by
  rw [tilt_eq c hk, Real.rpow_add hk]
  norm_num [spectrumWeight, weight]

theorem raw_eq_band (c : Cosmology) {L K : ℚ} (hL : 0 < L) (hLK : L ≤ K) (R : ℝ) :
    raw c L K R = ∫ k in (L:ℝ)..(K:ℝ), spectrumWeight c R k := by
  apply intervalIntegral.integral_congr
  intro k hk
  rw [Set.uIcc_of_le (by exact_mod_cast hLK)] at hk
  exact (weight_eq c R ((by exact_mod_cast hL : (0:ℝ) < L).trans_le hk.1)).symm

/-- Conservative additive tail allowance using a rational lower bound for R. -/
def budget (L K Rmin C₀ C₁ : ℚ) : ℚ :=
  C₀^2*L^3 + 18*C₁^2/(Rmin^4*K^2)

def pad (J : RationalInterval) (b : ℚ) : RationalInterval := ⟨J.lo,J.hi+b⟩

theorem pad_sound {c : Cosmology} {L K Rmin C₀ C₁ : ℚ} {R : ℝ} {J : RationalInterval}
    (hn : 0 ≤ c.ns ∧ c.ns ≤ 1) (hL : 0 < L) (hL1 : L ≤ 1) (hK : 1 ≤ K)
    (hrm : 0 < Rmin) (hr : (Rmin:ℝ) ≤ R) (hKR : 1 ≤ K*Rmin)
    (lo : ∀ k ∈ Ioc (0:ℝ) L, |transferReal (coefficients c) k| ≤ C₀)
    (hi : ∀ k ∈ Ioi (K:ℝ), |transferReal (coefficients c) k| ≤ C₁/k)
    (hint : IntervalIntegrable (spectrumWeight c R) volume L K)
    (hJ : J.Contains (raw c L K R)) :
    (pad J (budget L K Rmin C₀ C₁)).Contains (fullRaw c R) := by
  rw [raw_eq_band c hL (hL1.trans hK)] at hJ
  have hrmp : (0:ℝ) < Rmin := by exact_mod_cast hrm
  have hKp : (0:ℝ) < K := by exact_mod_cast zero_lt_one.trans_le hK
  have hkr : 1 ≤ (K:ℝ)*R := (by exact_mod_cast hKR : (1:ℝ) ≤ (K:ℝ)*(Rmin:ℝ)).trans
    (mul_le_mul_of_nonneg_left hr hKp.le)
  have hh := extend (measurableWeight c R)
    (show 0 ≤ (c.ns:ℝ) ∧ (c.ns:ℝ) ≤ 1 by exact_mod_cast hn)
    (hrmp.trans_le hr) (by exact_mod_cast hL) (by exact_mod_cast hL1)
    (by exact_mod_cast hK) hkr lo hi hint
  refine ⟨hJ.1.trans hh.1, hh.2.trans ?_⟩
  change _ ≤ ((J.hi + budget L K Rmin C₀ C₁ : ℚ):ℝ)
  unfold budget
  push_cast
  have hnum : 0 ≤ (18:ℝ)*(C₁:ℝ)^2 := by positivity
  have hden : 0 < (Rmin:ℝ)^4*(K:ℝ)^2 := by positivity
  have hfrac := div_le_div_of_nonneg_left hnum hden
    (show (Rmin:ℝ)^4*(K:ℝ)^2 ≤ R^4*(K:ℝ)^2 by gcongr)
  dsimp only [spectrumWeight] at hJ
  linarith [hJ.2]

end SashimiEhFull
