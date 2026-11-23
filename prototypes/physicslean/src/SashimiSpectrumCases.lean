import SashimiSpectrumVariance
import SashimiSpectrumFit
import SashimiIntegrable

/-! Regression inputs from SASHIMI-W's WMAP7 table and the analytic C/F/SI
prescriptions. σ₈ normalization is separate; this W cell is a contribution to
the unnormalized sharp-k variance inside its passband. -/

namespace SashimiSpectrum.Cases
open Real FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory Set
noncomputable section

/-- Rows 1206–1207 of WMAP7_camb_matterpower_z0_extrapolated.dat,
k in h/Mpc and P in (Mpc/h)³. -/
def wCellPower (k : ℝ) :=
  57.120377 + (k - 0.99999938) * (54.122997 - 57.120377) / (1.0232924 - 0.99999938)

def wCellPowerExpr (k : Expr) : Expr :=
  .binary .add (.const 57.120377)
    (.binary .mul (.binary .sub k (.const 0.99999938))
      (.const ((54.122997 - 57.120377) / (1.0232924 - 0.99999938))))

theorem eval_wCellPower (k : Expr) (env : Nat → ℝ) :
    (wCellPowerExpr k).eval env = wCellPower (k.eval env) := by
  simp [wCellPowerExpr, wCellPower, Expr.eval, BinaryOp.eval]
  norm_num
  ring

/-- Continuously apply published q5 to the linear-interpolated table.
The particle mass is 1.5 keV. h, Ωb, Ωm come from the file headers. -/
def wAlpha : ℝ := wdmAlpha 1.5 (0.27 - 0.0469) 0.7

def wCellIntegrand : Expr :=
  let k := Expr.unary .exp (.var 0)
  let a := wdmAlphaExpr (.const 1.5) (.const (0.27 - 0.0469)) (.const 0.7)
  varianceExpr k (.binary .mul (wCellPowerExpr k) (wdmPowerExpr .publishedQ5 a k)) (.const 1)

theorem wCellIntegrand_eq (s : ℝ) : (wCellIntegrand.eval fun _ => s) =
    delta2 (modifiedPower wCellPower (wdmPowerRatio .publishedQ5 wAlpha)) (exp s) := by
  unfold wCellIntegrand
  rw [eval_varianceExpr]
  simp only [Expr.eval, BinaryOp.eval, UnaryOp.eval, eval_wCellPower]
  have ha : ∀ (env : ℕ → ℝ), (wdmAlphaExpr (.const 1.5) (.const (0.27 - 0.0469)) (.const 0.7)).eval env = wAlpha := by
    intro env
    rw [eval_wdmAlphaExpr (by norm_num [Expr.eval]) (by norm_num [Expr.eval]) (by norm_num [Expr.eval])]
    norm_num [wAlpha, Expr.eval]
  have hpos : 0 < wAlpha := by unfold wAlpha wdmAlpha; positivity
  rw [eval_wdmPowerExpr .publishedQ5 (env := fun _ => s) (k := .unary .exp (.var 0))
    (by rw [ha]; exact hpos) (Real.exp_pos s), ha]
  simp [Expr.eval, UnaryOp.eval, delta2, modifiedPower]

/-- The chosen natural-log panel stays within the data cell. -/
theorem wCell_domain (s : ℝ) (h : s ∈ Icc 0 (1/50)) :
    0.99999938 ≤ exp s ∧ exp s ≤ 1.0232924 := by
  have h0 := h.1
  have h1 := h.2
  constructor <;> interval

set_option maxHeartbeats 4000000 in
/-- Both the finite integration error and rounding of every expression evaluation
are included in these rational bounds; P itself is the specified interpolant. -/
theorem wmap7_q5_cell :
    0.056 ≤ (∫ s in (0 : ℝ)..(1/50), delta2
      (modifiedPower wCellPower (wdmPowerRatio .publishedQ5 wAlpha)) (exp s)) ∧
    (∫ s in (0 : ℝ)..(1/50), delta2
      (modifiedPower wCellPower (wdmPowerRatio .publishedQ5 wAlpha)) (exp s)) ≤ 0.060 := by
  let I : RationalInterval := ⟨380077139/6710886400, 160359181/2684354560⟩
  have hh : SashimiNumerics.integrate? wCellIntegrand 0 (1/50) 1 28 (1/100) 12 = some I := by
    decide +kernel
  obtain ⟨h, _⟩ := SashimiNumerics.integrate_sound_auto hh
  simp only [wCellIntegrand_eq, Rat.cast_zero, Rat.cast_div, Rat.cast_one, Rat.cast_ofNat,
    RationalInterval.Contains] at h
  exact ⟨le_trans (by norm_num [I]) h.1, le_trans h.2 (by norm_num [I])⟩

/-- At αk=1, squaring the amplitude is a large and observable difference. -/
theorem wdm_conventions :
    (0.04530 ≤ wdmAmplitude 1 1 ∧ wdmAmplitude 1 1 ≤ 0.04531) ∧
    (0.00205 ≤ wdmPowerRatio .standardT2Q10 1 1 ∧ wdmPowerRatio .standardT2Q10 1 1 ≤ 0.00206) := by
  simp only [wdmAmplitude, wdmPowerRatio, WdmConvention.q, one_mul, Real.one_rpow]
  norm_num
  rw [Real.rpow_def_of_pos (by norm_num : (0:ℝ) < 2),
      Real.rpow_def_of_pos (by norm_num : (0:ℝ) < 2)]
  constructor <;> constructor <;> interval

theorem fdm_sample :
    0.5262 ≤ fdmPowerRatio 1 (9 / 2.22) ∧ fdmPowerRatio 1 (9 / 2.22) ≤ 0.5263 := by
  norm_num [fdmPowerRatio]
  constructor <;> interval

theorem ludlow_pivot :
    3.77928 ≤ ludlowSigma0 1 1 ∧ ludlowSigma0 1 1 ≤ 3.77929 := by
  norm_num [ludlowSigma0, ludlowLogShape, ludlowDen]

theorem ludlow_pivot_slope :
    -0.10762 ≤ ludlowLogSlope 0 ∧ ludlowLogSlope 0 ≤ -0.10761 := by
  norm_num [ludlowLogSlope, ludlowDen]

end
end SashimiSpectrum.Cases
