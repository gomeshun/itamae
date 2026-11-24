import SashimiSpectrumVariance

/-! Fixed two-epoch EPS conditional mass-fraction density per dS in Yang et al.,
Eq. (14). Both variances are at the common reference epoch; redshift occurs in
the two extrapolated spherical-collapse barriers. -/
namespace EpsKernel
open FloatLib.Numerics FloatLib.Numerics.Interval

noncomputable def density (d s : ℝ) : ℝ :=
  d * Real.exp (-(d^2) / (2*s)) / (Real.sqrt (2*Real.pi) * s ^ (3/2 : ℝ))

private def expression : Expr :=
  let d := Expr.var 0
  let s := Expr.var 1
  let exponent := Expr.binary .div (.unary .neg (.unary (.pow 2) d)) (.const 2 * s)
  .binary .div (d * .unary .exp exponent)
    (.unary .sqrt (.const 2 * SashimiSpectrum.piExpr) * ((.unary .sqrt s : Expr) * s))

private theorem eval_expression (v : Nat → ℝ) (hs : 0 < v 1) :
    expression.eval v = density (v 0) (v 1) := by
  change _ / (Real.sqrt (((2:ℚ):ℝ) * SashimiSpectrum.piExpr.eval v) * (Real.sqrt (v 1) * v 1)) = _
  rw [SashimiSpectrum.eval_piExpr]
  unfold density
  rw [show (3/2 : ℝ) = 1/2 + 1 by norm_num, Real.rpow_add_one hs.ne', ← Real.sqrt_eq_rpow]
  norm_num [Expr.eval, UnaryOp.eval, BinaryOp.eval]

/-- The barrier at the child epoch (higher z) comes first, followed by the barrier at the parent epoch and the two reference variances.
Unordered inputs or a nonpositive lower difference leave the result unresolved. -/
def enclose? (childBarrier parentBarrier child parent : RationalInterval) (bits degree : Nat) :
    Option RationalInterval := do
  let d := childBarrier.sub parentBarrier
  let s := child.sub parent
  guard (childBarrier.lo ≤ childBarrier.hi ∧ parentBarrier.lo ≤ parentBarrier.hi ∧ child.lo ≤ child.hi ∧ parent.lo ≤ parent.hi
    ∧ 0 < d.lo ∧ 0 < s.lo)
  expression.gridEnclosureWith? [d,s] bits degree

theorem sound {childBarrier parentBarrier child parent bits degree I} {b₁ b₀ s₁ s₀ : ℝ}
    (hl : childBarrier.Contains b₁) (he : parentBarrier.Contains b₀) (hc : child.Contains s₁)
    (hp : parent.Contains s₀) (h : enclose? childBarrier parentBarrier child parent bits degree = some I) :
    0 < b₁-b₀ ∧ 0 < s₁-s₀ ∧ I.Contains (density (b₁-b₀) (s₁-s₀)) := by
  unfold enclose? at h
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  have bounds : childBarrier.lo ≤ childBarrier.hi ∧ parentBarrier.lo ≤ parentBarrier.hi ∧ child.lo ≤ child.hi ∧
      parent.lo ≤ parent.hi ∧ 0 < (childBarrier.sub parentBarrier).lo ∧ 0 < (child.sub parent).lo := by
    by_contra hn
    simp [guard, hn] at hg
  have hd := RationalInterval.contains_sub hl he
  have hs := RationalInterval.contains_sub hc hp
  have dpos : (0 : ℝ) < (childBarrier.sub parentBarrier).lo := by exact_mod_cast bounds.2.2.2.2.1
  have spos : (0 : ℝ) < (child.sub parent).lo := by exact_mod_cast bounds.2.2.2.2.2
  let v : Nat → ℝ := fun | 0 => b₁-b₀ | _ => s₁-s₀
  have hr := Expr.contains_gridEnclosureWith? (values := v) (h := h) (by
    intro i box hi
    cases i with
    | zero => simp at hi; subst box; exact hd
    | succ i => cases i with
      | zero => simp at hi; subst box; exact hs
      | succ i => simp at hi)
  rw [eval_expression v (spos.trans_le hs.1)] at hr
  exact ⟨dpos.trans_le hd.1, spos.trans_le hs.1, hr⟩

#print axioms sound
end EpsKernel
