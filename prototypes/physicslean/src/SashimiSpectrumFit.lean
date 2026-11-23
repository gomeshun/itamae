import SashimiSpectrumExpr

/-! Cold default variance used in C/SI and the F base class. The derivative
is taken with respect to natural log mass; all decimal fit coefficients are
exact rationals for this specified empirical prescription. -/

namespace SashimiSpectrum
open Real FloatLib.Numerics.Interval
noncomputable section

def ludlowDen (u : ℝ) := 1 + 1.53 * exp (0.275 * u) + 3.36 * exp (0.198 * u)

def ludlowLogShape (u : ℝ) := 22.26 * exp (0.292 * u) / ludlowDen u

/-- The physical pivot mass is 10¹⁰ M_sun/h in the chosen mass units. -/
def ludlowSigma0 (pivot M : ℝ) := ludlowLogShape (log pivot - log M)

theorem ludlowSigma0_eq (pivot M : ℝ) (hp : 0 < pivot) (hm : 0 < M) :
    ludlowSigma0 pivot M =
      22.26 * (pivot / M) ^ (0.292 : ℝ) /
        (1 + 1.53 * (pivot/M) ^ (0.275 : ℝ) + 3.36 * (pivot/M) ^ (0.198 : ℝ)) := by
  simp [ludlowSigma0, ludlowLogShape, ludlowDen, Real.rpow_def_of_pos (div_pos hp hm),
    Real.log_div hp.ne' hm.ne', mul_comm]

def ludlowLogSlope (u : ℝ) :=
  -0.292 + (1.53 * 0.275 * exp (0.275*u) + 3.36 * 0.198 * exp (0.198*u)) / ludlowDen u

/-- d ln σ / d ln M; the factor of two for d ln S / d ln M follows from S=σ². -/
theorem ludlow_log_slope (pivot u : ℝ) :
    HasDerivAt (fun m => log (ludlowSigma0 pivot (exp m)))
      (ludlowLogSlope (log pivot - u)) u := by
  have de (a : ℝ) : HasDerivAt (fun m => exp (a * (log pivot - m)))
      (exp (a * (log pivot - u)) * (-a)) u := by
    convert ((((hasDerivAt_id' u).const_sub (log pivot)).const_mul a).exp) using 1; ring
  let v := log pivot - u
  have hden : 0 < ludlowDen v := by unfold ludlowDen; positivity
  have hn : 0 < 22.26 * exp (0.292*v) := by positivity
  have hd := (((de 0.275).const_mul (1.53 : ℝ)).const_add 1).add ((de 0.198).const_mul (3.36 : ℝ))
  have h := (((de 0.292).const_mul (22.26 : ℝ)).div hd hden.ne').log (div_pos hn hden).ne'
  simp only [ludlowSigma0, log_exp, ludlowLogShape, ludlowDen]
  convert h using 1
  dsimp [ludlowLogSlope, ludlowDen, v] at *
  field_simp
  ring

/-- Numerically evaluate the C/SI prescription through FloatLib. Units cancel
in the ratio pivot/M. Zero mass is excluded by the real semantics theorem. -/
def ludlowExpr (pivot M : Expr) : Expr :=
  let u := Expr.binary .sub (.unary .log pivot) (.unary .log M)
  let t := fun a => Expr.unary .exp (.binary .mul (.const a) u)
  .binary .div (.binary .mul (.const 22.26) (t 0.292))
    (.binary .add (.binary .add (.const 1) (.binary .mul (.const 1.53) (t 0.275)))
      (.binary .mul (.const 3.36) (t 0.198)))

theorem eval_ludlowExpr (pivot M : Expr) (env : Nat → ℝ) :
    (ludlowExpr pivot M).eval env = ludlowSigma0 (pivot.eval env) (M.eval env) := by
  simp [ludlowExpr, Expr.eval, BinaryOp.eval, UnaryOp.eval, ludlowSigma0, ludlowLogShape, ludlowDen]

end
end SashimiSpectrum
