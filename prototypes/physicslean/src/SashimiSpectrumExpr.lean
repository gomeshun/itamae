import SashimiSpectrum
import FloatLib.Numerics.Enclosure.Expression.ClosedProof

/-! FloatLib expressions for the windows and transfers on strictly positive
wavenumbers. Their real semantics are proved below. The integration API can
apply these to arbitrary expressions, including a unit conversion or exp(log k). -/

namespace SashimiSpectrum

open FloatLib.Numerics.Interval Real
open Expr

/-- Real exponent at a positive base, reified via exp and log. -/
def positivePowExpr (x y : Expr) := unary .exp (binary .mul (unary .log x) y)

theorem eval_positivePowExpr {x y : Expr} {env : Nat → ℝ} (hx : 0 < x.eval env) :
    (positivePowExpr x y).eval env = x.eval env ^ y.eval env := by
  exact (Real.rpow_def_of_pos hx _).symm

def cutoffExpr (kfs k : Expr) : Expr :=
  unary .exp (unary .neg (unary (.pow 2) (binary .div k kfs)))

theorem eval_cutoffExpr (kfs k : Expr) (env : Nat → ℝ) :
    (cutoffExpr kfs k).eval env = cdmCutoff (kfs.eval env) (k.eval env) := rfl

def topHatExpr (x : Expr) : Expr :=
  binary .div (binary .mul (.const 3) (binary .sub (unary .sin x)
    (binary .mul x (unary .cos x)))) (unary (.pow 3) x)

theorem eval_topHatExpr {x : Expr} {env : Nat → ℝ} (h : x.eval env ≠ 0) :
    (topHatExpr x).eval env = topHat (x.eval env) := by
  simp [topHat, h, topHatExpr, Expr.eval, BinaryOp.eval, UnaryOp.eval]

/-- W's alpha in Mpc/h, evaluated from the three source parameters. -/
def wdmAlphaExpr (mKeV Ωc h : Expr) : Expr :=
  let factor := fun x y p => positivePowExpr (.binary .div x y) (.const p)
  .binary .mul (.binary .mul (.binary .mul (.const 0.049)
    (factor (.const 1) mKeV 1.11)) (factor Ωc (.const 0.25) 0.11)) (factor h (.const 0.7) 1.22)

theorem eval_wdmAlphaExpr {m Ω h : Expr} {env : Nat → ℝ} (hm : 0 < m.eval env)
    (hΩ : 0 < Ω.eval env) (hh : 0 < h.eval env) :
    (wdmAlphaExpr m Ω h).eval env = wdmAlpha (m.eval env) (Ω.eval env) (h.eval env) := by
  simp only [wdmAlphaExpr, wdmAlpha, positivePowExpr, Expr.eval, BinaryOp.eval, UnaryOp.eval]
  push_cast
  norm_num
  repeat' rw [← Real.rpow_def_of_pos (by positivity)]
  rw [Real.rpow_def_of_pos (inv_pos.mpr hm), Real.log_inv, neg_mul]

def wdmPowerExpr (convention : WdmConvention) (alpha k : Expr) : Expr :=
  let q : ℚ := match convention with | .publishedQ5 => 5 | .standardT2Q10 => 10
  positivePowExpr (binary .add (.const 1)
    (positivePowExpr (binary .mul alpha k) (.const (2 * 1.12)))) (.const (-q / 1.12))

theorem eval_wdmPowerExpr (convention : WdmConvention) {alpha k : Expr} {env : Nat → ℝ}
    (ha : 0 < alpha.eval env) (hk : 0 < k.eval env) :
    (wdmPowerExpr convention alpha k).eval env =
      wdmPowerRatio convention (alpha.eval env) (k.eval env) := by
  have hab : 0 < (binary .mul alpha k).eval env := mul_pos ha hk
  have hp : 0 < (binary .add (.const 1)
      (positivePowExpr (binary .mul alpha k) (.const (2 * 1.12)))).eval env := by
    change 0 < 1 + (positivePowExpr _ _).eval env
    rw [eval_positivePowExpr hab]
    exact add_pos zero_lt_one (Real.rpow_pos_of_pos hab _)
  cases convention <;> unfold wdmPowerExpr wdmPowerRatio <;>
    rw [eval_positivePowExpr hp] <;>
    simp [Expr.eval, BinaryOp.eval, WdmConvention.q, eval_positivePowExpr hab]

def fdmPowerExpr (m22 k : Expr) : Expr :=
  let A := binary .mul (.const 2.22) (positivePowExpr m22
    (binary .sub (.const 0.04) (binary .div (unary .log m22) (.const 1000))))
  let B := binary .mul (.const 0.16) (positivePowExpr m22 (.const (-0.05)))
  let x := binary .div (binary .mul A k) (binary .mul (.const 9) (unary .sqrt m22))
  let t := positivePowExpr x (.const 2.5)
  unary (.pow 2) (binary .div (binary .div (unary .sin t) t)
    (binary .add (.const 1) (binary .mul B (positivePowExpr x (.const 3.5)))))

theorem eval_fdmPowerExpr {m22 k : Expr} {env : Nat → ℝ}
    (hm : 0 < m22.eval env) (hk : 0 < k.eval env) :
    (fdmPowerExpr m22 k).eval env = fdmPowerRatio (m22.eval env) (k.eval env) := by
  unfold fdmPowerExpr fdmPowerRatio
  simp only [Expr.eval, BinaryOp.eval, UnaryOp.eval, positivePowExpr]
  simp only [← Real.rpow_def_of_pos hm]
  push_cast
  norm_num
  repeat' rw [← Real.rpow_def_of_pos (by positivity)]
  rw [ite_eq_right (ne_of_gt (Real.rpow_pos_of_pos (by positivity) _))]

end SashimiSpectrum
