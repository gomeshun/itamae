import SashimiSpectrumExpr
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds

/-! Stable FDM knot enclosure. Taylor bounds avoid dividing a tiny sine by
an interval containing zero. At large arguments the unit sine bound supplies
a decaying envelope without resolving rapidly oscillating phase. -/
namespace SashimiSpectrum
open Real FloatLib.Numerics FloatLib.Numerics.Interval

def fdmX (m k : Expr) : Expr :=
  let a := Expr.const 2.22 * positivePowExpr m (.binary .sub (.const 0.04)
    (.binary .div (.unary .log m) (.const 1000)))
  .binary .div (a * k) (.const 9 * .unary .sqrt m)

def fdmT (m k : Expr) : Expr :=
  let x := fdmX m k
  .unary (.pow 2) x * .unary .sqrt x

def fdmD (m k : Expr) : Expr :=
  let b := Expr.const 0.16 * positivePowExpr m (.const (-0.05))
  let x := fdmX m k
  .const 1 + b * (.unary (.pow 3) x * .unary .sqrt x)

private noncomputable def sinc (t : ℝ) := if t = 0 then 1 else sin t / t

private theorem geometry {m k : Expr} {v : ℕ → ℝ} (hm : 0 < m.eval v) (hk : 0 ≤ k.eval v) :
    0 ≤ (fdmT m k).eval v ∧ 1 ≤ (fdmD m k).eval v ∧
    fdmPowerRatio (m.eval v) (k.eval v) = (sinc ((fdmT m k).eval v) / (fdmD m k).eval v)^2 := by
  let x := (fdmX m k).eval v
  have hx : 0 ≤ x := by
    dsimp [x, fdmX, positivePowExpr, Expr.eval, UnaryOp.eval, BinaryOp.eval]
    positivity
  have hpow (a : ℕ) : x ^ ((a : ℝ)+1/2) = x^a * sqrt x := by
    rw [Real.rpow_add_of_nonneg hx (by positivity) (by norm_num), Real.rpow_natCast, Real.sqrt_eq_rpow]
  have h5 := hpow 2
  have h7 := hpow 3
  norm_num1 at h5 h7
  have ht : (fdmT m k).eval v = x^2 * sqrt x := rfl
  have hd : (fdmD m k).eval v = 1 + (0.16 * m.eval v ^ (-0.05 : ℝ)) * (x^3 * sqrt x) := by
    simp [fdmD, Expr.eval, BinaryOp.eval, UnaryOp.eval, eval_positivePowExpr hm, x]
  have heqx : x = 2.22 * (m.eval v) ^ (0.04 - log (m.eval v) / 1000) * k.eval v / (9 * sqrt (m.eval v)) := by
    simp [x, fdmX, Expr.eval, BinaryOp.eval, UnaryOp.eval, eval_positivePowExpr hm]
  rw [ht, hd]
  refine ⟨by positivity, le_add_of_nonneg_right (by positivity), ?_⟩
  dsimp only [fdmPowerRatio]
  rw [← heqx]
  rw [show (2.5 : ℝ) = 5/2 by norm_num, show (3.5 : ℝ) = 7/2 by norm_num, h5, h7]
  rfl

private def lowExpr (t d : Expr) : Expr :=
  .unary (.pow 2) (.binary .div (.binary .max (.const 0)
    (.binary .sub (.const 1) (.binary .div (.unary (.pow 2) t) (.const 6)))) d)

private def highExpr (d : Expr) : Expr := .unary (.pow 2) (.binary .div (.const 1) d)

private theorem envelope {t d : Expr} {v : ℕ → ℝ} (ht : 0 ≤ t.eval v) (hd : 1 ≤ d.eval v) :
    (lowExpr t d).eval v ≤ (sinc (t.eval v)/d.eval v)^2 ∧
    (sinc (t.eval v)/d.eval v)^2 ≤ (highExpr d).eval v := by
  set z := t.eval v
  have hds : 0 < d.eval v := by linarith
  have hlow : 1 - z^2/6 ≤ sinc z := by
    rcases ht.eq_or_lt with h | h
    · simp [← h, sinc]
    · rw [sinc, ite_eq_right h.ne']
      apply (le_div_iff₀ h).2
      nlinarith [Real.sin_ge_sub_cube h.le]
  have hu : |sinc z| ≤ 1 := by
    rcases ht.eq_or_lt with h | h
    · simp [← h, sinc]
    · rw [sinc, ite_eq_right h.ne', abs_div, abs_of_nonneg h.le]
      apply (div_le_one h).2
      simpa [abs_of_nonneg h.le] using (Real.abs_sin_le_abs (x := z))
  change (max 0 (1 - z^2/6) / d.eval v)^2 ≤ _ ∧ _ ≤ (1/d.eval v)^2
  rw [div_pow, div_pow, div_pow, one_pow]
  constructor
  · apply div_le_div_of_nonneg_right _ (sq_nonneg _)
    rw [sq_le_sq, abs_of_nonneg (le_max_left _ _)]
    exact max_le (abs_nonneg _) (hlow.trans (le_abs_self _))
  · apply div_le_div_of_nonneg_right _ (sq_nonneg _)
    simpa only [one_pow] using (sq_le_sq.mpr (hu.trans_eq (abs_of_nonneg zero_le_one).symm))

/-- Compute certified squared transfer bounds for a single supplied spectrum
knot. Use direct reflection around the turnover and algebraic envelopes elsewhere. -/
def fdm? (m k : ℚ) (bits degree : Nat) : Option RationalInterval := do
  guard (0 < m ∧ 0 ≤ k)
  let t := fdmT (.const m) (.const k)
  let d := fdmD (.const m) (.const k)
  let lo ← (lowExpr t d).gridEnclosure? bits degree
  let hi ← (highExpr d).gridEnclosure? bits degree
  -- The choice affects width only; either enclosure covers the same transfer.
  if 0 < k ∧ m/16 < k ∧ k < 64*(1+m) then
    if let some mid := (fdmPowerExpr (.const m) (.const k)).gridEnclosure? bits degree then
      return ⟨max lo.lo mid.lo, min hi.hi mid.hi⟩
  return ⟨lo.lo, hi.hi⟩

theorem fdm_sound {m k : ℚ} {bits degree : ℕ} {I : RationalInterval}
    (h : fdm? m k bits degree = some I) : I.Contains (fdmPowerRatio m k) := by
  unfold fdm? at h
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  have hdom : 0 < m ∧ 0 ≤ k := by by_contra hn; simp [guard, hn] at hg
  obtain ⟨L, hl, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨R, hr, h⟩ := Option.bind_eq_some_iff.mp h
  have hm : 0 < (m : ℝ) := by exact_mod_cast hdom.1
  have hk : 0 ≤ (k : ℝ) := by exact_mod_cast hdom.2
  obtain ⟨ht,hd,he⟩ := geometry (m := .const m) (k := .const k) (v := fun _ => 0) hm hk
  obtain ⟨h1,h2⟩ := envelope ht hd
  have hl' := (Expr.contains_gridEnclosure? hl).1.trans h1
  have hr' := h2.trans (Expr.contains_gridEnclosure? hr).2
  simp only [Expr.eval] at he
  rw [← he] at hl' hr'
  split at h
  · rename_i hmid
    split at h
    · rename_i J hj
      cases h
      have hmids := Expr.contains_gridEnclosure? hj
      rw [eval_fdmPowerExpr (m22 := .const m) (k := .const k) (env := fun _ => 0) hm (by exact_mod_cast hmid.1 : (0:ℝ)<k)] at hmids
      simpa [RationalInterval.Contains] using And.intro (max_le hl' hmids.1) (le_min hr' hmids.2)
    · cases h; exact ⟨hl',hr'⟩
  · cases h; exact ⟨hl',hr'⟩

end SashimiSpectrum
