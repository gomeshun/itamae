import SashimiSpectrumExpr
import SashimiIntegrable
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals

/-! Whole-panel bounds and the two omitted tails of a finite spectrum table. -/

namespace SashimiSpectrum
open Real MeasureTheory Set FloatLib.Numerics FloatLib.Numerics.Interval
noncomputable section

/-- π represented by a closed expression with a verified arctangent enclosure. -/
def piExpr : Expr := .binary .mul (.const 4) (.unary .atan (.const 1))

@[simp] theorem eval_piExpr (env : Nat → ℝ) : piExpr.eval env = π := by
  change (4 : ℝ) * arctan (1 : ℚ) = π
  norm_num [Real.arctan_one]
  ring

/-- Window and power here are expressions already applied to k.
For sharp-k integration, restrict the endpoints before setting W=1. -/
def varianceExpr (k power window : Expr) : Expr :=
  .binary .mul (.binary .div (.binary .mul (.unary (.pow 3) k) power)
    (.binary .mul (.const 2) (.unary (.pow 2) piExpr))) (.unary (.pow 2) window)

theorem eval_varianceExpr (k power window : Expr) (env : Nat → ℝ) :
    (varianceExpr k power window).eval env =
      (k.eval env)^3 * power.eval env / (2 * π^2) * (window.eval env)^2 := by
  simp only [varianceExpr, Expr.eval, BinaryOp.eval, UnaryOp.eval, eval_piExpr]
  norm_num

/-- Reuse a complete grid enclosure of the interpolation, transfer and window
on each subinterval. All equalities only need to hold on the stated table band. -/
theorem windowVariance_integrate_sound {P W : ℝ → ℝ} {r scale : ℝ}
    {e : Expr} {a b tol : ℚ} {depth bits degree : ℕ} {I : RationalInterval}
    (hab : a ≤ b)
    (he : ∀ s ∈ Icc (a : ℝ) (b : ℝ),
      (e.eval fun _ => s) = windowIntegrand P W r scale s)
    (h : SashimiNumerics.integrate? e a b depth bits tol degree = some I) :
    I.Contains (windowVariance P W r scale a b) ∧ I.hi - I.lo ≤ tol := by
  have hi := SashimiNumerics.integrate_sound_auto h
  have hh : (∫ s in (a : ℝ)..(b : ℝ), e.eval fun _ => s) =
      windowVariance P W r scale a b := by
    apply intervalIntegral.integral_congr
    rw [uIcc_of_le (by exact_mod_cast hab)]
    exact he
  rwa [hh] at hi

/-- If only pointwise ranges of an external spectrum and the window are known,
their squared-window product still certifies each entire panel. -/
theorem panel_bounds {P W : ℝ → ℝ} {r scale s kl ku pl pu wl wu : ℝ}
    (hk : kl ≤ exp s ∧ exp s ≤ ku) (hkl : 0 ≤ kl)
    (hp : pl ≤ P (exp s) ∧ P (exp s) ≤ pu) (hpl : 0 ≤ pl)
    (hw : wl ≤ (W (exp s * r / scale))^2 ∧ (W (exp s * r / scale))^2 ≤ wu) (hwl : 0 ≤ wl) :
    kl^3 * pl / (2*π^2) * wl ≤ windowIntegrand P W r scale s ∧
      windowIntegrand P W r scale s ≤ ku^3 * pu / (2*π^2) * wu := by
  unfold windowIntegrand delta2
  have hP := hpl.trans hp.1
  have hpu := hP.trans hp.2
  have hku := (Real.exp_pos s).le.trans hk.2
  constructor
  · gcongr <;> first | exact hk.1 | exact hp.1 | exact hw.1
  · gcongr <;> first | exact hk.2 | exact hp.2 | exact hw.2

/-- Two integrable exponential envelopes in log k bound the omitted positive
tail contributions. A finite input table by itself supplies neither envelope. -/
theorem full_band_error {f : ℝ → ℝ} {a b cin cout : ℝ} (hf : Integrable f)
    (hf0 : ∀ t, 0 ≤ f t) (hlo : ∀ t ≤ a, f t ≤ cin * exp t)
    (hhi : ∀ t > b, f t ≤ cout * exp (-t)) :
    0 ≤ (∫ t, f t) - (∫ t in a..b, f t) ∧
      (∫ t, f t) - (∫ t in a..b, f t) ≤ cin * exp a + cout * exp (-b) := by
  have hall := intervalIntegral.integral_Iic_add_Ioi
    (b := b) hf.integrableOn hf.integrableOn
  have hband := intervalIntegral.integral_Iic_sub_Iic
    (a := a) (b := b) hf.integrableOn hf.integrableOn
  have heq : (∫ t, f t) - (∫ t in a..b, f t) =
      (∫ t in Iic a, f t) + ∫ t in Ioi b, f t := by linarith
  rw [heq]
  constructor
  · exact add_nonneg (integral_nonneg hf0) (integral_nonneg hf0)
  · apply add_le_add
    · calc (∫ t in Iic a, f t) ≤ ∫ t in Iic a, cin * exp t :=
            setIntegral_mono_on hf.integrableOn ((integrableOn_exp_Iic a).const_mul cin)
              measurableSet_Iic hlo
          _ = _ := by rw [integral_const_mul, integral_exp_Iic]
    · calc (∫ t in Ioi b, f t) ≤ ∫ t in Ioi b, cout * exp (-t) :=
            setIntegral_mono_on hf.integrableOn ((integrableOn_exp_neg_Ioi b).const_mul cout)
              measurableSet_Ioi hhi
          _ = _ := by rw [integral_const_mul, integral_exp_neg_Ioi]


/-- A positive-k cell written in dk and the same cell in d ln k give exactly
the same integral, with the k² and k³ factors fixed by the Jacobian. -/
theorem windowVariance_dk (P W : ℝ → ℝ) (r c a b : ℝ) :
    windowVariance P W r c a b =
      ∫ k in exp a..exp b, k^2 * P k * W (k*r/c)^2 / (2*π^2) := by
  let g := fun k : ℝ => k^2 * P k * W (k*r/c)^2 / (2*π^2)
  have h := intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
    (a := a) (b := b) (g := g) continuous_exp.continuousOn
    (fun s _ => (hasDerivAt_exp s)) (fun s _ => (exp_pos s).le)
  apply Eq.trans _ h
  apply intervalIntegral.integral_congr
  intro s _
  dsimp [windowVariance, windowIntegrand, delta2, g]
  ring

/-- At fixed redshift the consistent variance derivative scales as D².
The historical F implementation instead returns D times the redshift-zero
derivative; those outputs coincide at D=1. -/
theorem variance_mass_derivative {S : ℝ → ℝ} {M dS D : ℝ} (h : HasDerivAt S dS M) :
    HasDerivAt (fun m => variance D (S m)) (D^2 * dS) M :=
  h.const_mul _

/-- Logarithmic slope of the standard deviation with respect to ln(1+z).
This specializes to the analytic derivative of whichever growth fit is supplied. -/
theorem redshift_log_sigma {D : ℝ → ℝ} {z D' S0 : ℝ} (h : HasDerivAt D D' z)
    (hz : 0 < 1+z) (hd : 0 < D z) (hs : 0 < S0) :
    HasDerivAt (fun u => log (sigma (D (exp u - 1)) S0))
      ((1+z)*D' / D z) (log (1+z)) := by
  have hm : exp (log (1+z)) - 1 = z := by rw [exp_log hz]; ring
  have hg := (h.mul_const (sqrt S0)).log (mul_pos hd (sqrt_pos.mpr hs)).ne'
  have hc := HasDerivAt.comp (log (1+z))
    (h := fun u => exp u - 1)
    (show HasDerivAt (fun t => log (D t * sqrt S0)) _ (exp (log (1+z))-1) from by rw [hm]; exact hg)
    ((hasDerivAt_exp _).sub_const 1)
  rw [exp_log hz] at hc
  convert hc using 1
  · rfl
  · field_simp [hd.ne', (sqrt_pos.mpr hs).ne']

/-- The sharp-k derivative in the k_c convention printed by ITAMAE/F.
The finite-band clamp has kinks at the edges: use only for a strict interior cutoff. -/
theorem sharpVariance_derivative (P : ℝ → ℝ) (ρ c a b M : ℝ) (hr : 0 < ρ)
    (hc : 0 < c) (hm : 0 < M) (hP : Continuous (fun s => delta2 P (exp s)))
    (hlo : a < logCutoff ρ c M) (hhi : logCutoff ρ c M < b) :
    deriv (sharpVariance P ρ c a b) M =
      -(c/massRadius ρ M)^3 * P (c/massRadius ρ M) / (6 * π^2 * M) := by
  rw [(sharpVariance_hasDerivAt P ρ c a b M hm hP hlo hhi).deriv,
    exp_logCutoff ρ c M hr hc hm]
  unfold delta2
  ring


end
end SashimiSpectrum
