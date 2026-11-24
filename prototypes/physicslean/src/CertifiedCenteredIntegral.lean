import FloatLib.Numerics.Enclosure.Rational.Proof
import Mathlib.Analysis.Calculus.MeanValue
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-! A midpoint integral enclosure from a point evaluation and a bound on the
derivative throughout a closed cell. Constant derivative contributions cancel
regardless of their size or sign. -/
namespace CertifiedCenteredIntegral
open Set MeasureTheory FloatLib.Numerics

/-- Rational endpoints for an ordered cell, midpoint value box and whole-cell
derivative box. The radius is zero for a singleton derivative interval. -/
def cell (a b : ℚ) (F D : RationalInterval) : RationalInterval :=
  let width := b - a
  let radius := (D.hi - D.lo) * width^2 / 8
  ⟨width * F.lo - radius, width * F.hi + radius⟩

theorem integrable {f d : ℝ → ℝ} {a b : ℝ} (hab : a ≤ b)
    (hd : ∀ x ∈ Icc a b, HasDerivWithinAt f (d x) (Icc a b) x) :
    IntervalIntegrable f volume a b :=
  ContinuousOn.intervalIntegrable_of_Icc hab fun x hx => (hd x hx).continuousWithinAt

theorem absMoment {a b : ℝ} (hab : a ≤ b) :
    (∫ x in a..b, |x - (a+b)/2|) = (b-a)^2 / 4 := by
  let c := (a+b)/2
  have hac : a ≤ c := by dsimp [c]; linarith
  have hcb : c ≤ b := by dsimp [c]; linarith
  have hc : Continuous (fun x : ℝ => |x-c|) := (continuous_id.sub continuous_const).abs
  change (∫ x in a..b, |x-c|) = _
  rw [← intervalIntegral.integral_add_adjacent_intervals (hc.intervalIntegrable a c)
    (hc.intervalIntegrable c b)]
  have hleft : (∫ x in a..c, |x-c|) = ∫ x in a..c, c-x := by
    apply intervalIntegral.integral_congr
    intro x hx
    rw [uIcc_of_le hac] at hx
    change |x-c| = c-x
    rw [abs_of_nonpos (sub_nonpos.mpr hx.2), neg_sub]
  have hright : (∫ x in c..b, |x-c|) = ∫ x in c..b, x-c := by
    apply intervalIntegral.integral_congr
    intro x hx
    rw [uIcc_of_le hcb] at hx
    exact abs_of_nonneg (sub_nonneg.mpr hx.1)
  rw [hleft, hright, intervalIntegral.integral_sub (f := fun _ : ℝ => c) (g := fun x : ℝ => x)
    intervalIntegrable_const (continuous_id.intervalIntegrable a c),
    intervalIntegral.integral_sub (f := fun x : ℝ => x) (g := fun _ : ℝ => c)
      (continuous_id.intervalIntegrable c b) intervalIntegrable_const,
    intervalIntegral.integral_const, intervalIntegral.integral_const, integral_id, integral_id]
  simp only [smul_eq_mul, c]
  ring

/-- The variation in the derivative bounds the error of midpoint quadrature.
A derivative within the cell suffices, including one-sided endpoint derivatives. -/
theorem midpoint_error {f d : ℝ → ℝ} {a b l u : ℝ} (hab : a ≤ b)
    (hd : ∀ x ∈ Icc a b, HasDerivWithinAt f (d x) (Icc a b) x)
    (hdu : ∀ x ∈ Icc a b, l ≤ d x ∧ d x ≤ u) :
    |(∫ x in a..b, f x) - (b-a)*f ((a+b)/2)| ≤ (u-l)*(b-a)^2 / 8 := by
  let c := (a+b)/2
  let s := (l+u)/2
  let g := fun x => f x - s*x
  have hc : c ∈ Icc a b := by dsimp [c]; constructor <;> linarith
  have hg (x : ℝ) (hx : x ∈ Icc a b) : HasDerivWithinAt g (d x-s*1) (Icc a b) x :=
    (hd x hx).sub ((hasDerivAt_id x).const_mul s).hasDerivWithinAt
  have hb (x : ℝ) (hx : x ∈ Icc a b) : ‖d x-s*1‖ ≤ (u-l)/2 := by
    obtain ⟨h1,h2⟩ := hdu x hx
    rw [Real.norm_eq_abs, abs_le]
    dsimp [s]
    constructor <;> linarith
  have point (x : ℝ) (hx : x ∈ Icc a b) : ‖g x - g c‖ ≤ (u-l)/2 * |x-c| :=
    (convex_Icc a b).norm_image_sub_le_of_norm_hasDerivWithin_le hg hb hc hx
  have integrated := intervalIntegral.norm_integral_le_of_norm_le (μ := volume)
    (f := fun x => g x - g c) hab (Filter.Eventually.of_forall fun x hx => point x ⟨hx.1.le,hx.2⟩)
    (((continuous_id.sub continuous_const).abs.const_mul ((u-l)/2)).intervalIntegrable a b)
  rw [intervalIntegral.integral_const_mul, absMoment hab, Real.norm_eq_abs,
    intervalIntegral.integral_sub (integrable hab hg) intervalIntegrable_const] at integrated
  change |(∫ x in a..b, f x - s*x) - _| ≤ _ at integrated
  rw [intervalIntegral.integral_sub (f := f) (g := fun x => s*x) (integrable hab hd)
    ((continuous_const.mul continuous_id).intervalIntegrable a b),
    intervalIntegral.integral_const_mul s (fun x : ℝ => x), integral_id, intervalIntegral.integral_const] at integrated
  have heq : (∫ x in a..b, f x) - s*((b^2-a^2)/2) - (b-a) • g c =
      (∫ x in a..b, f x) - (b-a)*f c := by
    simp only [g, c, smul_eq_mul]
    ring
  rw [heq] at integrated
  exact integrated.trans_eq (by ring)

/-- Transfer the closed real bound to the computable rational cell enclosure. -/
theorem cell_sound {f d : ℝ → ℝ} {a b : ℚ} {F D : RationalInterval} (hab : a ≤ b)
    (hF : F.Contains (f (((a+b)/2 : ℚ) : ℝ)))
    (hd : ∀ x ∈ Icc (a:ℝ) b, HasDerivWithinAt f (d x) (Icc (a:ℝ) b) x)
    (hD : ∀ x ∈ Icc (a:ℝ) b, D.Contains (d x)) :
    (cell a b F D).Contains (∫ x in (a:ℝ)..(b:ℝ), f x) := by
  have hab' : (a:ℝ) ≤ b := by exact_mod_cast hab
  have he := abs_le.mp (midpoint_error hab' hd hD)
  obtain ⟨hlo,hhi⟩ := hF
  push_cast at hlo hhi
  have hw := sub_nonneg.mpr hab'
  have hl := mul_le_mul_of_nonneg_left hlo hw
  have hu := mul_le_mul_of_nonneg_left hhi hw
  unfold cell RationalInterval.Contains
  push_cast
  constructor <;> linarith [he.1,he.2]

#print axioms midpoint_error
#print axioms cell_sound
end CertifiedCenteredIntegral
