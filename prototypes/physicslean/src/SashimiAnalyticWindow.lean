import SashimiTableWindow
import SashimiNormalization

/-! Integral enclosures for a continuous-spectrum formula on a specified finite k band.
Coefficient enclosures are shared by every cell. No interpolation between sampled powers
is used. The numerator and σ₈ denominator use the same weighted spectrum. -/
namespace SashimiAnalyticWindow

open FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory Set SashimiSpectrum

/-- Index 0 is the wavenumber, higher indices stand for constant model coefficients. -/
def environment (values : Nat → ℝ) (k : ℝ) : Nat → ℝ
  | 0 => k
  | i+1 => values i

def Covers (boxes : List RationalInterval) (values : Nat → ℝ) : Prop :=
  ∀ (i : Nat) (I : RationalInterval), boxes[i]? = some I → I.Contains (values i)

theorem measurableFormula (env : Nat → ℝ → ℝ) (hv : ∀ i, Measurable (env i)) :
    ∀ e : Expr, Measurable (fun k : ℝ => e.eval fun i => env i k)
  | .const _ => by simpa only [Expr.eval] using (measurable_const : Measurable (fun _ : ℝ => _))
  | .var i => hv i
  | .call _ _ => by simpa only [Expr.eval] using (measurable_const : Measurable (fun _ : ℝ => (0 : ℝ)))
  | .unary op a => by
    have ha := measurableFormula env hv a
    cases op <;> simp only [Expr.eval, UnaryOp.eval, Real.tan_eq_sin_div_cos,
      Real.tanh_eq_sinh_div_cosh] <;> fun_prop
  | .binary op a b => by
    have ha := measurableFormula env hv a
    have hb := measurableFormula env hv b
    cases op <;> simp only [Expr.eval, BinaryOp.eval] <;> fun_prop
  | .ternary op e f g => by
    cases op
    exact ((measurableFormula env hv e).mul (measurableFormula env hv f)).add (measurableFormula env hv g)
termination_by structural e => e

noncomputable def integrand (weighted : Expr) (coeff : Nat → ℝ) (radius k : ℝ) : ℝ :=
  weighted.eval (environment coeff k) * topHat (k*radius)^2

theorem measurableIntegrand (e : Expr) (c : Nat → ℝ) (r : ℝ) :
    Measurable (integrand e c r) := by
  have hf : Measurable (fun k : ℝ => e.eval (environment c k)) := by
    apply measurableFormula (fun i k => environment c k i)
    intro i
    cases i with
    | zero => exact measurable_id
    | succ n => exact measurable_const
  have hw : Measurable topHat := by
    unfold topHat
    exact Measurable.ite (measurableSet_singleton 0) measurable_const (by fun_prop)
  exact hf.mul ((hw.comp (measurable_id.mul measurable_const)).pow_const 2)

/-- Enclose weighted power and the top-hat window throughout one positive k cell. -/
def range? (weighted : Expr) (boxes : List RationalInterval) (R : RationalInterval)
    (bits degree : Nat) (l u : ℚ) : Option RationalInterval :=
  if 0 ≤ l ∧ 0 ≤ R.lo then
    (weighted.gridEnclosureWith? (⟨l,u⟩ :: boxes) bits degree).map fun P =>
      P.mul (SashimiTableWindow.windowSq (l*R.lo) (u*R.hi) bits degree)
  else none

theorem range_sound {e boxes R bits degree c r} (hc : Covers boxes c) (hr : R.Contains r)
    (l u : ℚ) (I : RationalInterval) (h : range? e boxes R bits degree l u = some I)
    (k : ℝ) (hk : k ∈ Icc (l:ℝ) u) : I.Contains (integrand e c r k) := by
  unfold range? at h
  by_cases hd : 0 ≤ l ∧ 0 ≤ R.lo
  · rw [ite_eq_left hd] at h
    obtain ⟨P, hP, hh⟩ := Option.map_eq_some_iff.mp h
    cases hh
    apply RationalInterval.contains_mul
    · apply Expr.contains_gridEnclosureWith? (h := hP)
      intro i J hJ
      cases i with
      | zero => simp at hJ; subst J; exact hk
      | succ n => exact hc n J (by simpa using hJ)
    · have hl : 0 ≤ (l:ℝ) := by exact_mod_cast hd.1
      have hR : 0 ≤ (R.lo:ℝ) := by exact_mod_cast hd.2
      apply SashimiTableWindow.windowSq_sound (mul_nonneg hd.1 hd.2)
      · push_cast
        exact mul_le_mul hk.1 hr.1 hR (hl.trans hk.1)
      · push_cast
        exact mul_le_mul hk.2 hr.2 (hR.trans hr.1) ((hl.trans hk.1).trans hk.2)
  · simp [hd] at h

def segment? (weighted : Expr) (boxes : List RationalInterval) (R : RationalInterval)
    (depth bits degree : Nat) (a b : ℚ) : Option RationalInterval :=
  if a ≤ b then SashimiNumerics.rectangles? (range? weighted boxes R bits degree) a b depth
  else none

theorem segment_sound {e boxes R depth bits degree a b I c r}
    (hc : Covers boxes c) (hr : R.Contains r) (h : segment? e boxes R depth bits degree a b = some I) :
    IntervalIntegrable (integrand e c r) volume a b ∧
      I.Contains (∫ k in (a:ℝ)..b, integrand e c r k) := by
  unfold segment? at h
  split at h
  · rename_i hle
    have hf := SashimiNumerics.rectangles_integrable (measurableIntegrand e c r)
      (range_sound hc hr) hle h
    exact ⟨hf, SashimiNumerics.rectangles_sound (range_sound hc hr) hle hf h⟩
  · contradiction

/-- Nonuniform mesh outside recursively bisected cells. Each endpoint names the next
closed interval; no tails or intervals between successive endpoints are omitted. -/
def band? (weighted : Expr) (boxes : List RationalInterval) (R : RationalInterval)
    (depth bits degree : Nat) (a : ℚ) : List ℚ → Option RationalInterval
  | [] => some (.point 0)
  | b :: rest => do
    let x ← segment? weighted boxes R depth bits degree a b
    let y ← band? weighted boxes R depth bits degree b rest
    return x.add y

theorem band_sound {e boxes R depth bits degree a c r I} (hc : Covers boxes c) (hr : R.Contains r)
    (ends : List ℚ) (h : band? e boxes R depth bits degree a ends = some I) :
    IntervalIntegrable (integrand e c r) volume a (ends.getLastD a) ∧
    I.Contains (∫ k in (a:ℝ)..(ends.getLastD a), integrand e c r k) := by
  induction ends generalizing a I with
  | nil =>
    cases h
    simp [List.getLastD, RationalInterval.point, RationalInterval.Contains]
  | cons b rest ih =>
    simp only [band?] at h
    obtain ⟨U, hU, hrem⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨V, hV, hh⟩ := Option.bind_eq_some_iff.mp hrem
    cases hh
    have h1 := segment_sound hc hr hU
    have h2 := ih hV
    rw [List.getLastD_cons]
    exact ⟨h1.1.trans h2.1, (intervalIntegral.integral_add_adjacent_intervals h1.1 h2.1) ▸
      RationalInterval.contains_add h1.2 h2.2⟩

#print axioms band_sound

end SashimiAnalyticWindow
