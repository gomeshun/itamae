import SashimiNumerics
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.Tactic.Ring

/-! The Carroll–Press–Turner fitting function used in the four SASHIMI
clients. Numerical certification here is relative to this fit. Its accuracy as
a solution of the perturbation equation is a separate modelling question. -/
namespace SashimiGrowth

open FloatLib.Numerics FloatLib.Numerics.Interval

def omega {K : Type*} [Field K] (m l z : K) : K :=
  1 - l / (l + m * (1 + z)^3)

def omegaPrime {K : Type*} [Field K] (m l z : K) : K :=
  l * (3*m*(1+z)^2) / (l + m*(1+z)^3)^2

def phiRest {K : Type*} [Field K] (x y : K) : K :=
  -y + (1+x/2)*(1+y/70)

noncomputable def phi (x y : ℝ) : ℝ := x ^ (4/7 : ℝ) + phiRest x y

def slopeRest {K : Type*} [Field K] (x : K) : K :=
  1 + (1+(1-x)/70)/2 - (1+x/2)/70

noncomputable def phiSlope (x : ℝ) : ℝ := (4/7) * x ^ (-3/7 : ℝ) + slopeRest x

/-- The legacy W fit retains its own ΩΛ at the reference time, which can
cause D(0) to differ from one. Set l=1-m for normalized native ITAMAE. -/
noncomputable def growth (m l z : ℝ) : ℝ :=
  let x := omega m l z
  (phi m l / m) * (x / (phi x (1-x) * (1+z)))

noncomputable def rateAt (m l z x : ℝ) : ℝ :=
  let dx := omegaPrime m l z
  let p := phi x (1-x)
  let v := p * (1+z)
  (phi m l / m) * ((dx * v - x * ((phiSlope x * dx) * (1+z) + p)) / v^2)

/-- Exact analytic derivative; the conditions exclude singularities of the fit. -/
theorem growth_hasDerivAt (m l z : ℝ)
    (he : l + m * (1+z)^3 ≠ 0) (hx : omega m l z ≠ 0)
    (hp : phi (omega m l z) (1-omega m l z) * (1+z) ≠ 0) :
    HasDerivAt (growth m l) (rateAt m l z (omega m l z)) z := by
  let x := omega m l z
  have ha : HasDerivAt (fun t : ℝ => 1+t) 1 z := by simpa using (hasDerivAt_id z).const_add 1
  have de := ((ha.pow 3).const_mul m).const_add l
  have dx : HasDerivAt (omega m l) (omegaPrime m l z) z := by
    convert! (hasDerivAt_const z (1 : ℝ)).sub ((hasDerivAt_const z l).div de he) using 1
    unfold omegaPrime
    simp only [Pi.pow_apply]
    ring
  have dy := (hasDerivAt_const z (1 : ℝ)).sub dx
  have dr := dx.rpow_const (p := (4/7 : ℝ)) (Or.inl hx)
  have dp : HasDerivAt (fun t => phi (omega m l t) (1-omega m l t)) (phiSlope x * omegaPrime m l z) z := by
    convert! dr.add (dy.neg.add (((dx.div_const 2).const_add 1).mul ((dy.div_const 70).const_add 1))) using 1
    unfold phiSlope slopeRest x
    norm_num
    ring
  convert! ((dx.div (dp.mul ha) hp).const_mul (phi m l / m)) using 1
  simp [rateAt, x]

inductive Quantity where
  | growth | growthDerivative | legacyFDerivative | collapseThreshold
deriving Repr

noncomputable def value (q : Quantity) (m l h z : ℝ) : ℝ :=
  match q with
  | .growth => growth m l z
  | .growthDerivative => rateAt m l z (omega m l z)
  | .legacyFDerivative => rateAt m l z (omega (m/h^2) l z)
  | .collapseThreshold => (1686/1000) / growth m l z

/-- The complementary-density convention in SASHIMI equals the ratio used in ITAMAE. -/
theorem omega_eq_matter_fraction {m l z : ℝ} (he : l + m*(1+z)^3 ≠ 0) :
    omega m l z = m*(1+z)^3 / (l + m*(1+z)^3) := by
  unfold omega
  field_simp
  ring

/-- ITAMAE's native formula evaluates the CPT suppression g and divides it by
g₀(1+z). This identity connects the SASHIMI evaluation order to that definition. -/
theorem growth_eq_cpt (m l z : ℝ) :
    let x := omega m l z
    growth m l z = ((5/2)*x/phi x (1-x)) / (((5/2)*m/phi m l)*(1+z)) := by
  unfold growth
  simp only [div_eq_mul_inv, mul_inv_rev, inv_inv]
  ring

/-- With flat reference fractions and nonzero matter, the CPT normalization gives D(0)=1. -/
theorem growth_zero {m : ℝ} (hm : 0 < m) (hm1 : m ≤ 1) : growth m (1-m) 0 = 1 := by
  have ho : omega m (1-m) 0 = m := by simp [omega]
  have hpos : 0 < phi m (1-m) := by
    have hr := Real.rpow_pos_of_pos hm (4/7 : ℝ)
    have hs := mul_nonneg hm.le (sub_nonneg.mpr hm1)
    unfold phi phiRest
    nlinarith
  simp [growth, ho]
  field_simp

/-- The only transcendental in the fit at rational arguments is a real
power of a positive rational. exp(log(x)*p) is enclosed by FloatLib. -/
private def powerExpr (x p : ℚ) : Expr :=
  .unary .exp (.unary .log (.const x) * .const p)

private def phiExpr (x y : ℚ) : Expr :=
  powerExpr x (4/7) + .const (phiRest x y)

private def slopeExpr (x : ℚ) : Expr :=
  .const (4/7) * powerExpr x (-3/7) + .const (slopeRest x)

def expression (q : Quantity) (m l h z : ℚ) : Expr :=
  let x := omega (if q matches .legacyFDerivative then m/h^2 else m) l z
  let p := phiExpr x (1-x)
  let v := p * .const (1+z)
  let c := .binary .div (phiExpr m l) (.const m)
  let D := c * .binary .div (.const x) v
  match q with
  | .growth => D
  | .collapseThreshold => .binary .div (.const (1686/1000)) D
  | .growthDerivative | .legacyFDerivative =>
    let dx := omegaPrime m l z
    let dv := slopeExpr x * .const dx * .const (1+z) + p
    c * .binary .div (.binary .sub (.const dx * v) (.const x * dv)) (.unary (.pow 2) v)

private theorem eval_mul (a b : Expr) (v : Nat → ℝ) : (a*b).eval v = a.eval v * b.eval v := rfl
private theorem eval_add (a b : Expr) (v : Nat → ℝ) : (a+b).eval v = a.eval v + b.eval v := rfl

private theorem eval_power {x : ℚ} (p : ℚ) (hx : 0 < x) (v : Nat → ℝ) :
    (powerExpr x p).eval v = (x : ℝ) ^ (p : ℝ) := by
  rw [Real.rpow_def_of_pos (by exact_mod_cast hx)]
  rfl

private theorem eval_phi {x : ℚ} (y : ℚ) (hx : 0 < x) (v : Nat → ℝ) :
    (phiExpr x y).eval v = phi x y := by
  simp only [phiExpr, eval_add, eval_power _ hx, Expr.eval]
  unfold phi phiRest
  push_cast
  rfl

private theorem eval_slope {x : ℚ} (hx : 0 < x) (v : Nat → ℝ) :
    (slopeExpr x).eval v = phiSlope x := by
  simp only [slopeExpr, eval_add, eval_power _ hx, Expr.eval]
  unfold phiSlope slopeRest
  push_cast
  rfl

def matterArgument (q : Quantity) (m l h z : ℚ) : ℚ :=
  omega (if q matches .legacyFDerivative then m/h^2 else m) l z

theorem expression_value (q : Quantity) (m l h z : ℚ) (hm : 0 < m)
    (hx : 0 < matterArgument q m l h z) :
    (expression q m l h z).eval (fun _ => 0) = value q m l h z := by
  cases q <;> simp only [matterArgument,↓reduceIte, Bool.false_eq_true] at hx
  all_goals simp only [expression, ↓reduceIte, Bool.false_eq_true,
    eval_phi _ hx, eval_phi _ hm, eval_slope hx, Expr.eval, BinaryOp.eval, UnaryOp.eval]
  all_goals simp [value, growth, rateAt, omega, omegaPrime]

def enclose (q : Quantity) (m l h z tol : ℚ) (bits : Nat := 64)
    (degree : Nat := 24) : Option RationalInterval :=
  if 0 < m ∧ m ≤ 1 ∧ 0 ≤ l ∧ 0 < h ∧ 0 ≤ z ∧ 0 < matterArgument q m l h z then
    SashimiNumerics.closed? (expression q m l h z) bits tol degree
  else none

theorem enclose_sound {q m l h z tol bits degree I}
    (hc : enclose q m l h z tol bits degree = some I) :
    I.Contains (value q m l h z) ∧ I.hi - I.lo ≤ tol := by
  unfold enclose at hc
  split at hc
  · rename_i valid
    obtain ⟨hm, _, _, _, _, hx⟩ := valid
    exact expression_value _ _ _ _ _ hm hx ▸ SashimiNumerics.closed_sound hc
  · contradiction

/-- Two-point growth ratio for W, or the secant approximation selected by the
F consistent adapter. The caller supplies the exact stored redshifts; the
secant's approximation error versus the analytic derivative is separate. -/
inductive TwoPoint where
  | ratio | secant
deriving Repr

def twoPoint? (q : TwoPoint) (m l lo hi tol : ℚ) (bits : Nat := 64)
    (degree : Nat := 24) : Option RationalInterval :=
  if 0 < m ∧ m ≤ 1 ∧ 0 ≤ l ∧ -1 < lo ∧ -1 < hi ∧ 0 < omega m l lo ∧ 0 < omega m l hi ∧
      ((q matches .secant) → lo ≠ hi) then
    let a := expression .growth m l 1 lo
    let b := expression .growth m l 1 hi
    let e := match q with
      | .ratio => Expr.binary .div b a
      | .secant => .binary .div (.binary .sub b a) (.const (hi-lo))
    SashimiNumerics.closed? e bits tol degree
  else none

theorem twoPoint_sound {q m l lo hi tol bits degree I}
    (hc : twoPoint? q m l lo hi tol bits degree = some I) :
    I.Contains (match q with
      | .ratio => growth m l hi / growth m l lo
      | .secant => (growth m l hi - growth m l lo) / ((hi : ℝ)-lo)) ∧ I.hi - I.lo ≤ tol := by
  cases q <;> simp only [twoPoint?] at hc
  all_goals
    split at hc
    · rename_i valid
      obtain ⟨hm, _, _, _, _, hlo, hhi, _⟩ := valid
      have ha := expression_value .growth m l 1 lo hm hlo
      have hb := expression_value .growth m l 1 hi hm hhi
      obtain ⟨h, hw⟩ := SashimiNumerics.closed_sound hc
      exact ⟨by simpa only [Expr.eval, BinaryOp.eval, ha, hb, value, Rat.cast_sub] using h, hw⟩
    · contradiction

#print axioms twoPoint_sound
#print axioms growth_hasDerivAt
#print axioms enclose_sound

end SashimiGrowth
