import FloatLib.Tactic.Interval
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-! Linear spectra and finite-band mass variance for the SASHIMI clients.
Wavenumbers are in Mpc⁻¹ and P in Mpc³. The supplied CDM spectrum is a
mathematical input (for example, a specified table interpolant). -/

namespace SashimiSpectrum

open Real MeasureTheory Set Filter
open scoped Topology

noncomputable section

/-- Dimensionless power per unit natural logarithm of wavenumber. -/
def delta2 (P : ℝ → ℝ) (k : ℝ) := k^3 * P k / (2 * π^2)

/-- Apply a power ratio; an amplitude is squared before this step. -/
def modifiedPower (P ratio : ℝ → ℝ) (k : ℝ) := P k * ratio k

/-- Gaussian free-streaming cutoff used by the C prompt-cusp path. -/
def cdmCutoff (kfs k : ℝ) := exp (-(k / kfs)^2)

inductive WdmConvention where | publishedQ5 | standardT2Q10 deriving DecidableEq

def WdmConvention.q : WdmConvention → ℝ
  | .publishedQ5 => 5
  | .standardT2Q10 => 10

/-- SASHIMI-W uses its nonbaryonic Ωc for the warm component in the all-warm convention.
The calibration parameter is the warm-species density, with particle mass in keV
and h dimensionless. The result is in Mpc/h. -/
def wdmAlpha (mKeV Ωc h : ℝ) : ℝ :=
  0.049 * (1 / mKeV) ^ (1.11 : ℝ) * (Ωc / 0.25) ^ (0.11 : ℝ) * (h / 0.7) ^ (1.22 : ℝ)

/-- Convert alpha from Mpc/h to Mpc before using canonical k in Mpc⁻¹. -/
def wdmPowerRatio (convention : WdmConvention) (alpha k : ℝ) : ℝ :=
  (1 + (alpha * k) ^ (2 * (1.12 : ℝ))) ^ (-convention.q / (1.12 : ℝ))

/-- Viel transfer AMPLITUDE, before squaring to obtain a power ratio.
It is numerically the factor used directly as a power ratio by legacy SASHIMI-W. -/
def wdmAmplitude (alpha k : ℝ) : ℝ :=
  (1 + (alpha * k) ^ (2 * (1.12 : ℝ))) ^ (-(5 : ℝ) / 1.12)

theorem published_ratio_eq_amplitude (alpha k : ℝ) :
    wdmPowerRatio .publishedQ5 alpha k = wdmAmplitude alpha k := rfl

/-- W's shared half-mode k; P/P₀ at this scale is 1/2 for q5 and 1/4 for q10. -/
def halfMode (alpha : ℝ) : ℝ :=
  ((2 : ℝ) ^ ((1.12 : ℝ) / 5) - 1) ^ (1 / (2 * (1.12 : ℝ))) / alpha

/-- Squaring the published q5 power ratio gives the alternate q10 convention. -/
theorem q10_eq_q5_sq (alpha k : ℝ) (ha : 0 ≤ alpha) (hk : 0 ≤ k) :
    wdmPowerRatio .standardT2Q10 alpha k = wdmPowerRatio .publishedQ5 alpha k ^ 2 := by
  have hb : 0 ≤ 1 + (alpha * k) ^ (2 * (1.12 : ℝ)) := by positivity
  unfold wdmPowerRatio WdmConvention.q
  rw [← Real.rpow_mul_natCast hb]
  norm_num

theorem halfMode_q5 (alpha : ℝ) (ha : 0 < alpha) :
    wdmPowerRatio .publishedQ5 alpha (halfMode alpha) = 1/2 := by
  have hpow : 0 < (2 : ℝ) ^ ((1.12 : ℝ) / 5) - 1 := by
    exact sub_pos.mpr (Real.one_lt_rpow (by norm_num) (by norm_num))
  unfold wdmPowerRatio halfMode WdmConvention.q
  rw [mul_div_cancel₀ _ ha.ne', ← Real.rpow_mul hpow.le]
  norm_num [← Real.rpow_mul (by norm_num : (0 : ℝ) ≤ 2)]

theorem halfMode_q10 (alpha : ℝ) (ha : 0 < alpha) :
    wdmPowerRatio .standardT2Q10 alpha (halfMode alpha) = 1/4 := by
  rw [q10_eq_q5_sq _ _ ha.le (by
    unfold halfMode
    apply div_nonneg _ ha.le
    apply Real.rpow_nonneg
    exact sub_nonneg.mpr (Real.one_le_rpow (by norm_num) (by norm_num))), halfMode_q5 alpha ha]
  norm_num

theorem standard_ratio_eq_amplitude_sq (alpha k : ℝ) (ha : 0 ≤ alpha) (hk : 0 ≤ k) :
    wdmPowerRatio .standardT2Q10 alpha k = wdmAmplitude alpha k ^ 2 :=
  q10_eq_q5_sq alpha k ha hk

theorem halfMode_amplitude (alpha : ℝ) (ha : 0 < alpha) :
    wdmAmplitude alpha (halfMode alpha) = 1/2 := halfMode_q5 alpha ha

/-- The half-POWER wavenumber for the standard amplitude-squared prescription. -/
def standardPowerHalfMode (alpha : ℝ) : ℝ :=
  ((2 : ℝ) ^ ((1.12 : ℝ) / 10) - 1) ^ (1 / (2 * (1.12 : ℝ))) / alpha

theorem standard_power_half (alpha : ℝ) (ha : 0 < alpha) :
    wdmPowerRatio .standardT2Q10 alpha (standardPowerHalfMode alpha) = 1/2 := by
  have hp : 0 < (2 : ℝ) ^ ((1.12 : ℝ)/10) - 1 :=
    sub_pos.mpr (Real.one_lt_rpow (by norm_num) (by norm_num))
  unfold wdmPowerRatio standardPowerHalfMode WdmConvention.q
  rw [mul_div_cancel₀ _ ha.ne', ← Real.rpow_mul hp.le]
  norm_num [← Real.rpow_mul (by norm_num : (0 : ℝ) ≤ 2)]

/-- Nadler et al. FDM power prescription used by F. m22 is in 10⁻²²eV
and k in Mpc⁻¹. Use sinc(0)=1 as in the NumPy implementation. -/
def fdmPowerRatio (m22 k : ℝ) : ℝ :=
  let A := 2.22 * m22 ^ (0.04 - log m22 / 1000)
  let B := 0.16 * m22 ^ (-0.05 : ℝ)
  let x := A * k / (9 * sqrt m22)
  let t := x ^ (2.5 : ℝ)
  ((if t = 0 then 1 else sin t / t) / (1 + B * x ^ (3.5 : ℝ))) ^ 2

/-- Fourier-space spherical top-hat, continued at zero. -/
def topHat (x : ℝ) := if x = 0 then 1 else 3 * (sin x - x * cos x) / x^3

/-- Sharp-k window includes the cutoff point. -/
def sharpK (x : ℝ) : ℝ := if x ≤ 1 then 1 else 0

def massRadius (ρ M : ℝ) := (3 * M / (4 * π * ρ)) ^ (1/3 : ℝ)

def windowIntegrand (P W : ℝ → ℝ) (radius scale s : ℝ) :=
  delta2 P (exp s) * W (exp s * radius / scale)^2

def windowVariance (P W : ℝ → ℝ) (radius scale logLo logHi : ℝ) : ℝ :=
  ∫ s in logLo..logHi, windowIntegrand P W radius scale s

/-- k_c = c/R. The log implementation makes the moving-boundary derivative
explicit. rho and M must be in consistent mass/volume and mass units. -/
def logCutoff (ρ c M : ℝ) := log c - (log (3 / (4 * π * ρ)) + log M) / 3

theorem exp_logCutoff (ρ c M : ℝ) (hr : 0 < ρ) (hc : 0 < c) (hm : 0 < M) :
    exp (logCutoff ρ c M) = c / massRadius ρ M := by
  have hs : (0 : ℝ) < 3 / (4 * π * ρ) := by positivity
  unfold logCutoff massRadius
  have harg : 3 * M / (4 * π * ρ) = (3 / (4 * π * ρ)) * M := by ring
  rw [harg, Real.rpow_def_of_pos (mul_pos hs hm), exp_sub, exp_log hc, log_mul hs.ne' hm.ne']
  ring_nf

/-- Finite tabulation band; a cutoff below the lower edge has zero variance.
The two strict interior inequalities in the derivative exclude boundary kinks. -/
def sharpVariance (P : ℝ → ℝ) (ρ c logLo logHi M : ℝ) : ℝ :=
  ∫ s in logLo..(max logLo (min logHi (logCutoff ρ c M))), delta2 P (exp s)

theorem sharpVariance_hasDerivAt (P : ℝ → ℝ) (ρ c a b M : ℝ) (hM : 0 < M)
    (hc : Continuous (fun s => delta2 P (exp s)))
    (hlo : a < logCutoff ρ c M) (hhi : logCutoff ρ c M < b) :
    HasDerivAt (sharpVariance P ρ c a b)
      (-delta2 P (exp (logCutoff ρ c M)) / (3*M)) M := by
  have hd : HasDerivAt (logCutoff ρ c) (-(1/M)/3) M := by
    unfold logCutoff
    convert (((Real.hasDerivAt_log hM.ne').const_add (log (3 / (4 * π * ρ)))).div_const 3).const_sub (log c) using 1; ring
  have hnear : (fun m => ∫ s in a..logCutoff ρ c m, delta2 P (exp s)) =ᶠ[𝓝 M]
      sharpVariance P ρ c a b := by
    filter_upwards [hd.continuousAt.eventually_lt_const hhi, hd.continuousAt.eventually_const_lt hlo] with m hm hm'
    simp [sharpVariance, min_eq_right hm.le, max_eq_right hm'.le]
  have h := (hc.integral_hasStrictDerivAt a _).hasDerivAt.comp M hd
  convert h.congr_of_eventuallyEq hnear.symm using 1; ring

/-- Propagate growth into the finite-band variance, before taking its square root. -/
def variance (growth S0 : ℝ) := growth^2 * S0
def sigma (growth S0 : ℝ) := growth * sqrt S0

theorem sigma_sq (growth S0 : ℝ) (h : 0 ≤ S0) :
    sigma growth S0 ^ 2 = variance growth S0 := by
  simp [sigma, variance, mul_pow, Real.sq_sqrt h]

/-- d ln sigma / d ln M or d ln sigma / d ln(1+z) given the
corresponding variance derivative. x is M or (1+z), respectively. -/
theorem logSigma_hasDerivAt {S : ℝ → ℝ} {x S' : ℝ}
    (h : HasDerivAt S S' x) (hx : 0 < x) (hs : 0 < S x) :
    HasDerivAt (fun u => log (sqrt (S (exp u)))) (x * S' / (2 * S x)) (log x) := by
  have hx' : exp (log x) = x := exp_log hx
  have hlog := (h.sqrt hs.ne').log (Real.sqrt_pos.mpr hs).ne'
  have hc := HasDerivAt.comp (log x)
    (show HasDerivAt (fun y => log (sqrt (S y))) _ (exp (log x)) from hx'.symm ▸ hlog)
    (Real.hasDerivAt_exp (log x))
  simp only [hx'] at hc
  convert hc using 1
  · rfl
  · field_simp [hs.ne', (Real.sqrt_pos.mpr hs).ne']
    nlinarith [Real.sq_sqrt hs.le, congrArg (fun t : ℝ => x * S' * t) (Real.sq_sqrt hs.le)]

/-- W uses σ8 to normalize the entire variance, using the same finite-band model at M8. -/
def normalize8 (target rawAt8 raw : ℝ) := target^2 * raw / rawAt8

theorem normalize8_anchor (target raw : ℝ) (h : raw ≠ 0) :
    normalize8 target raw raw = target^2 := by unfold normalize8; field_simp

end
end SashimiSpectrum
