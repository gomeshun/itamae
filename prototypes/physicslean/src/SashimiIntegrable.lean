import SashimiNumerics
import Mathlib.MeasureTheory.Function.SpecialFunctions.Basic

/-! Successful enclosures provide a bound on each closed cell. Combined with
measurability of the FloatLib expression this implies Lebesgue integrability,
including for expressions whose operations are discontinuous elsewhere. -/

namespace SashimiNumerics

open FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory Set

/-- Use FloatLib's default interpretation of registered functions and the integration variable. -/
theorem measurable_eval : ∀ e : Expr, Measurable (fun t : ℝ => e.eval fun _ => t)
  | .const _ => by simpa only [Expr.eval] using (measurable_const : Measurable (fun _ : ℝ => _))
  | .var _ => by simp only [Expr.eval]; exact measurable_id
  | .call _ _ => by simpa only [Expr.eval] using (measurable_const : Measurable (fun _ : ℝ => (0 : ℝ)))
  | .unary op a => by
    have ha := measurable_eval a
    cases op <;> simp only [Expr.eval, UnaryOp.eval, Real.tan_eq_sin_div_cos,
      Real.tanh_eq_sinh_div_cosh] <;> fun_prop
  | .binary op a b => by
    have ha := measurable_eval a
    have hb := measurable_eval b
    cases op <;> simp only [Expr.eval, BinaryOp.eval] <;> fun_prop
  | .ternary op a b c => by
    have ha := measurable_eval a
    have hb := measurable_eval b
    have hc := measurable_eval c
    cases op
    exact (ha.mul hb).add hc
termination_by structural e => e

/-- A measurable function bounded by every accepted leaf is integrable across the whole mesh. -/
theorem rectangles_integrable {f : ℝ → ℝ} (hm : Measurable f)
    {range : ℚ → ℚ → Option RationalInterval}
    (hbox : ∀ a b I, range a b = some I → ∀ t ∈ Icc (a : ℝ) (b : ℝ), I.Contains (f t))
    {n : Nat} {a b : ℚ} {I : RationalInterval} (hab : a ≤ b)
    (h : rectangles? range a b n = some I) : IntervalIntegrable f volume a b := by
  induction n generalizing a b I with
  | zero =>
    obtain ⟨J, hj, _⟩ := Option.map_eq_some_iff.mp h
    have hab' : (a : ℝ) ≤ b := by exact_mod_cast hab
    rw [intervalIntegrable_iff_integrableOn_Icc_of_le hab']
    have hC : IntegrableOn (fun _ : ℝ => max (-(J.lo : ℝ)) (J.hi : ℝ))
        (Icc (a : ℝ) (b : ℝ)) volume := integrableOn_const (isCompact_Icc.measure_lt_top.ne)
    apply hC.mono' hm.aestronglyMeasurable
    filter_upwards [ae_restrict_mem (μ := volume) (s := Icc (a : ℝ) (b : ℝ)) measurableSet_Icc] with t ht
    obtain ⟨hl, hu⟩ := hbox _ _ _ hj t ht
    rw [Real.norm_eq_abs, abs_le]
    exact ⟨le_trans (neg_le_neg (le_max_left _ _)) (by simpa using hl), hu.trans (le_max_right _ _)⟩
  | succ n ih =>
    obtain ⟨J, hj, rest⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨K, hk, _⟩ := Option.bind_eq_some_iff.mp rest
    exact (ih (by linarith) hj).trans (ih (by linarith) hk)

theorem integrable_of_integrate {e : Expr} {a b tol : ℚ} {depth bits degree : Nat}
    {I : RationalInterval} (h : integrate? e a b depth bits tol degree = some I) :
    IntervalIntegrable (fun t => e.eval fun _ => t) volume a b := by
  unfold integrate? at h
  split at h
  · obtain ⟨J, hj, _⟩ := Option.bind_eq_some_iff.mp h
    exact rectangles_integrable (measurable_eval e) (exprRange_sound e bits degree) ‹a ≤ b› hj
  · contradiction

/-- A successful execution certifies a finite integral and its requested width. -/
theorem integrate_sound_auto {e : Expr} {a b tol : ℚ} {depth bits degree : Nat}
    {I : RationalInterval} (h : integrate? e a b depth bits tol degree = some I) :
    I.Contains (∫ t in (a : ℝ)..(b : ℝ), e.eval fun _ => t) ∧ I.hi - I.lo ≤ tol :=
  integrate_sound (integrable_of_integrate h) h

#print axioms integrate_sound_auto

end SashimiNumerics
