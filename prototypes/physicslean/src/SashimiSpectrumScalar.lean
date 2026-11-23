import SashimiSpectrumFit
import SashimiSpectrumVariance

/-! Scalar frontend for table-independent spectrum quantities. String and arity
validation select a closed expression, whose evaluation and final width are checked
by FloatLib. Power ratios and windows include their finite limits at k=0. -/

namespace SashimiSpectrum
open FloatLib.Numerics FloatLib.Numerics.Interval

def halfModeExpr (alpha : Expr) : Expr :=
  let x := positivePowExpr (.const 2) (.const (1.12/5))
  .binary .div (positivePowExpr (.binary .sub x (.const 1)) (.const (1/(2*1.12)))) alpha

theorem eval_halfModeExpr (alpha : Expr) (env : ℕ → ℝ) :
    (halfModeExpr alpha).eval env = halfMode (alpha.eval env) := by
  unfold halfModeExpr halfMode
  simp only [Expr.eval, BinaryOp.eval, positivePowExpr, UnaryOp.eval]
  push_cast
  rw [← Real.rpow_def_of_pos (by norm_num : (0:ℝ) < 2),
    ← Real.rpow_def_of_pos (sub_pos.mpr (Real.one_lt_rpow (by norm_num) (by norm_num)))]

def standardPowerHalfModeExpr (alpha : Expr) : Expr :=
  let x := positivePowExpr (.const 2) (.const (1.12/10))
  .binary .div (positivePowExpr (.binary .sub x (.const 1)) (.const (1/(2*1.12)))) alpha

theorem eval_standardPowerHalfModeExpr (alpha : Expr) (env : ℕ → ℝ) :
    (standardPowerHalfModeExpr alpha).eval env = standardPowerHalfMode (alpha.eval env) := by
  unfold standardPowerHalfModeExpr standardPowerHalfMode
  simp only [Expr.eval, BinaryOp.eval, positivePowExpr, UnaryOp.eval]
  push_cast
  rw [← Real.rpow_def_of_pos (by norm_num : (0:ℝ) < 2),
    ← Real.rpow_def_of_pos (sub_pos.mpr (Real.one_lt_rpow (by norm_num) (by norm_num)))]

def ludlowSlopeExpr (pivot M : Expr) : Expr :=
  let u := Expr.binary .sub (.unary .log pivot) (.unary .log M)
  let t := fun a => Expr.unary .exp (.binary .mul (.const a) u)
  let p := Expr.binary .mul (.const 1.53) (t 0.275)
  let q := Expr.binary .mul (.const 3.36) (t 0.198)
  .binary .add (.const (-0.292)) (.binary .div
    (.binary .add (.binary .mul (.const 0.275) p) (.binary .mul (.const 0.198) q))
    (.binary .add (.binary .add (.const 1) p) q))

theorem eval_ludlowSlopeExpr (pivot M : Expr) (env : ℕ → ℝ) :
    (ludlowSlopeExpr pivot M).eval env =
      ludlowLogSlope (Real.log (pivot.eval env) - Real.log (M.eval env)) := by
  simp [ludlowSlopeExpr, Expr.eval, BinaryOp.eval, UnaryOp.eval, ludlowLogSlope, ludlowDen]
  ring

theorem wdm_zero (convention : WdmConvention) (alpha : ℝ) :
    wdmPowerRatio convention alpha 0 = 1 := by
  norm_num [wdmPowerRatio, Real.zero_rpow]

theorem fdm_zero (m : ℝ) : fdmPowerRatio m 0 = 1 := by
  norm_num [fdmPowerRatio, Real.zero_rpow]

def scalarExpression? (op : String) (args : List ℚ) : Option Expr := do
  match op, args with
  | "cdm-cutoff", [kfs,k] =>
    guard (0 < kfs ∧ 0 ≤ k)
    return cutoffExpr (.const kfs) (.const k)
  | "wdm-alpha", [m,Ω,h] =>
    guard (0 < m ∧ 0 < Ω ∧ 0 < h)
    return wdmAlphaExpr (.const m) (.const Ω) (.const h)
  | "wdm-amplitude", [alpha,k] | "wdm-q5", [alpha,k] | "wdm-q10", [alpha,k] =>
    guard (0 < alpha ∧ 0 ≤ k)
    let q := if op = "wdm-q10" then WdmConvention.standardT2Q10 else .publishedQ5
    return if k = 0 then .const 1 else wdmPowerExpr q (.const alpha) (.const k)
  | "wdm-half-mode", [alpha] =>
    guard (0 < alpha)
    return halfModeExpr (.const alpha)
  | "wdm-standard-power-half", [alpha] =>
    guard (0 < alpha)
    return standardPowerHalfModeExpr (.const alpha)
  | "fdm-ratio", [m,k] =>
    guard (0 < m ∧ 0 ≤ k)
    return if k = 0 then .const 1 else fdmPowerExpr (.const m) (.const k)
  | "top-hat", [x] =>
    guard (0 ≤ x)
    return if x = 0 then .const 1 else topHatExpr (.const x)
  | "sharp-k", [x] =>
    guard (0 ≤ x)
    return .const (if x ≤ 1 then 1 else 0)
  | "delta2", [k,p] =>
    guard (0 ≤ k ∧ 0 ≤ p)
    return .binary .div (.const (k^3*p)) (.binary .mul (.const 2) (.unary (.pow 2) piExpr))
  | "variance", [d,s0] =>
    guard (0 ≤ d ∧ 0 ≤ s0)
    return .const (d^2*s0)
  | "sigma", [d,s0] =>
    guard (0 ≤ d ∧ 0 ≤ s0)
    return .binary .mul (.const d) (.unary .sqrt (.const s0))
  | "ludlow-sigma", [pivot,m] | "ludlow-slope", [pivot,m] =>
    guard (0 < pivot ∧ 0 < m)
    return if op = "ludlow-sigma" then ludlowExpr (.const pivot) (.const m)
      else ludlowSlopeExpr (.const pivot) (.const m)
  | _, _ => none

def scalar? (op : String) (args : List ℚ) (bits : Nat) (tol : ℚ) (degree : Nat) :
    Option RationalInterval := do
  let e ← scalarExpression? op args
  SashimiNumerics.closed? e bits tol degree

theorem scalar_sound {op : String} {args : List ℚ} {bits degree : ℕ} {tol : ℚ}
    {I : RationalInterval} (h : scalar? op args bits tol degree = some I) :
    ∃ e, scalarExpression? op args = some e ∧ I.Contains (e.eval fun _ => 0) ∧ I.hi - I.lo ≤ tol := by
  obtain ⟨e, he, hi⟩ := Option.bind_eq_some_iff.mp h
  exact ⟨e, he, SashimiNumerics.closed_sound hi⟩

end SashimiSpectrum
