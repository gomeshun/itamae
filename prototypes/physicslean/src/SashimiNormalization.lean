import SashimiTable
import SashimiSpectrumWindow
import FloatLib.Numerics.Enclosure.Rational.Converges

/-! Common σ₈ normalization on tabulated power. The current WDM prescription
uses the top-hat variance on the full supplied spectrum at R=8/h physical Mpc.
The target may be a signed derivative or a moment on a different band. -/

namespace SashimiNormalization

open FloatLib.Numerics FloatLib.Numerics.Interval SashimiTable

/-- Normalize an already enclosed raw quantity with an independent variance
enclosure. Both the positive lower bound and the final width are checked. -/
def apply? (target tol : ℚ) (raw anchor : RationalInterval) : Option RationalInterval := do
  guard (0 < target ∧ 0 < anchor.lo)
  SashimiNumerics.within? tol ((raw.scaleNonnegative (target^2)).divPositive anchor)

theorem apply_sound {target tol : ℚ} {raw anchor I : RationalInterval} {x v : ℝ}
    (hx : raw.Contains x) (hv : anchor.Contains v) (h : apply? target tol raw anchor = some I) :
    I.Contains (SashimiSpectrum.normalize8 target v x) ∧ I.hi - I.lo ≤ tol := by
  unfold apply? at h
  by_cases hg : 0 < target ∧ 0 < anchor.lo
  · simp only [hg] at h
    obtain ⟨rfl, hw⟩ := SashimiNumerics.within_sound h
    refine ⟨?_, hw⟩
    simpa only [Rat.cast_pow, SashimiSpectrum.normalize8] using
      RationalInterval.contains_divPositive
        (RationalInterval.contains_scaleNonnegative hx (sq_nonneg target)) hv hg.2
  · simp [hg] at h

/-- Defer the width gate until after the ratio; the band may be the entire table. -/
def band? (f : Cell → Expr) (mode : Interpolation) (rows : List Knot) (a b : ℚ)
    (depth bits degree : Nat) : Option RationalInterval :=
  if accepts mode rows a b then piecewise? f a b depth bits degree (cells rows) else none

theorem band_sound {f mode rows a b depth bits degree I}
    (h : band? f mode rows a b depth bits degree = some I) :
    I.Contains (band f a b (cells rows)) := by
  unfold band? at h
  split at h
  · exact piecewise_sound h
  · contradiction

/-- Defer the width gate when normalizing a power lookup. -/
def point? (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (rows : List Knot) (k : ℚ) (bits degree : Nat) : Option RationalInterval := do
  guard (accepts mode rows k k)
  let c ← (cells rows).find? (fun c => decide (c.left.k ≤ k ∧ k ≤ c.right.k))
  (powerExpr mode placement ratio c (.const k)).gridEnclosure? bits degree

theorem point_sound {mode placement ratio rows k bits degree I}
    (h : point? mode placement ratio rows k bits degree = some I) :
    ∃ c, (cells rows).find? (fun c => decide (c.left.k ≤ k ∧ k ≤ c.right.k)) = some c ∧
      I.Contains ((powerExpr mode placement ratio c (.const k)).eval fun _ => 0) := by
  unfold point? at h
  obtain ⟨_, _, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨c, hc, h⟩ := Option.bind_eq_some_iff.mp h
  exact ⟨c, hc, Expr.contains_gridEnclosure? h⟩

/-- The spherical top-hat stays in [-1,1] and lies above 1-x²/10.
Intersecting these two analytic bounds with its sine/cosine expression
avoids spurious small-k growth when integrating an entire WMAP7 table. -/
def stableWindow (x : Expr) : Expr :=
  .binary .max (.const (-1)) (.binary .max
    (.binary .sub (.const 1) (.binary .div (x.unary (.pow 2)) (.const 10)))
    (.binary .min (.const 1) (SashimiSpectrum.topHatExpr x)))

theorem eval_stableWindow {x : Expr} {env : Nat → ℝ} (h : 0 < x.eval env) :
    (stableWindow x).eval env = SashimiSpectrum.topHat (x.eval env) := by
  obtain ⟨hl, hu⟩ := SashimiSpectrum.topHat_bounds _ h.le
  have ht := SashimiSpectrum.topHat_taylor_lower _ h.le
  change max ((-1:ℚ):ℝ) (max (((1:ℚ):ℝ) - (x.eval env)^2 / (10:ℚ))
    (min ((1:ℚ):ℝ) ((SashimiSpectrum.topHatExpr x).eval env))) = _
  push_cast
  rw [SashimiSpectrum.eval_topHatExpr h.ne', min_eq_right hu, max_eq_right ht, max_eq_right hl]

/-- The anchor uses the same interpolated power and transfer as the target.
For a sharp-k anchor the upper integration limit is the chosen cutoff. -/
def anchorMoment (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (radius : ℚ) : Cell → Expr :=
  moment 0 mode placement ratio (if radius = 0 then fun _ => Expr.const 1
    else fun k => stableWindow (k * Expr.const radius))

/-- Identity with the physical real-space spherical window (or the selected
sharp-k band), for the same transfer and base interpolation. -/
theorem eval_anchorMoment (mode : Interpolation) (placement : Placement)
    (ratio : Expr → Expr) (f : ℝ → ℝ) (c : Cell) (radius : ℚ) (k : ℝ)
    (hr : 0 ≤ radius) (hk : 0 < k)
    (h : ∀ x : Expr, (ratio x).eval (fun _ => k) = f (x.eval fun _ => k)) :
    (anchorMoment mode placement ratio radius c).eval (fun _ => k) =
      momentReal 0 mode placement f
        (fun t => if radius = 0 then 1 else SashimiSpectrum.topHat (t * radius)) c k := by
  apply eval_moment _ _ _ _ _ _ _ _ _ (h _) (h _) (h _)
  by_cases hz : radius = 0
  · simp [hz, Expr.eval]
  · simp only [hz, ↓reduceIte]
    have hr' : 0 < (radius : ℝ) := by exact_mod_cast lt_of_le_of_ne hr (Ne.symm hz)
    rw [eval_stableWindow]
    · rfl
    · change 0 < k * (radius : ℝ)
      exact mul_pos hk hr'

/-- Compute the normalization variance on its own integration band. -/
def anchor? (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (rows : List Knot) (radius a b : ℚ) (depth bits degree : Nat) : Option RationalInterval :=
  band? (anchorMoment mode placement ratio radius) mode rows a b depth bits degree

theorem anchor_sound {mode placement ratio rows radius a b depth bits degree I}
    (h : anchor? mode placement ratio rows radius a b depth bits degree = some I) :
    I.Contains (band (anchorMoment mode placement ratio radius) a b (cells rows)) :=
  band_sound h

/-- Kernel-checked normalization of a supplied functional on its own finite band. -/
def normalized? (target tol : ℚ) (raw : RationalInterval) (mode : Interpolation)
    (placement : Placement) (ratio : Expr → Expr) (rows : List Knot)
    (radius a b : ℚ) (depth bits degree : Nat) : Option RationalInterval := do
  let anchor ← anchor? mode placement ratio rows radius a b depth bits degree
  apply? target tol raw anchor

theorem normalized_sound {target tol raw mode placement ratio rows radius a b depth bits degree I}
    {x : ℝ} (hx : raw.Contains x)
    (h : normalized? target tol raw mode placement ratio rows radius a b depth bits degree = some I) :
    I.Contains (SashimiSpectrum.normalize8 target
      (band (anchorMoment mode placement ratio radius) a b (cells rows)) x) ∧ I.hi - I.lo ≤ tol := by
  obtain ⟨J, hJ, hI⟩ := Option.bind_eq_some_iff.mp h
  exact apply_sound hx (anchor_sound hJ) hI

#print axioms normalized_sound

end SashimiNormalization
