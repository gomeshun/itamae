import FloatLib.Numerics.Enclosure.Expression.ClosedProof
import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.Tactic.Linarith

/-! Certified rational enclosures for scalar expressions and finite integrals.
Every accepted cell encloses its entire integrand range. The width gate acts on
the cumulative integral enclosure, including the outward rounding of the cells. -/

namespace SashimiNumerics

open FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory Set

/-- Return only enclosures meeting the absolute width budget. -/
def within? (tol : ℚ) (I : RationalInterval) : Option RationalInterval :=
  if I.hi - I.lo ≤ tol then some I else none

theorem within_sound {tol : ℚ} {I J : RationalInterval} (h : within? tol I = some J) :
    I = J ∧ J.hi - J.lo ≤ tol := by
  unfold within? at h
  split at h <;> simp_all

/-- Enclose a closed FloatLib expression on a binary grid with a rational width budget. -/
def closed? (e : Expr) (bits : Nat) (tol : ℚ) (degree : Nat := 24) : Option RationalInterval :=
  (e.gridEnclosure? bits degree).bind (within? tol)

theorem closed_sound {e : Expr} {bits degree : Nat} {tol : ℚ} {I : RationalInterval}
    (h : closed? e bits tol degree = some I) :
    I.Contains (e.eval fun _ => 0) ∧ I.hi - I.lo ≤ tol := by
  obtain ⟨J, hj, hJ⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨rfl, hw⟩ := within_sound hJ
  exact ⟨Expr.contains_gridEnclosure? hj, hw⟩

/-- Darboux bounds over any real mesh, including a nonuniform mesh. -/
theorem integral_grid_bounds {f : ℝ → ℝ} {x lo hi : ℕ → ℝ} {n : ℕ}
    (hx : Monotone x) (hf : ∀ i < n, IntervalIntegrable f volume (x i) (x (i + 1)))
    (hbox : ∀ i < n, ∀ t ∈ Icc (x i) (x (i + 1)), lo i ≤ f t ∧ f t ≤ hi i) :
    (∑ i ∈ Finset.range n, (x (i + 1) - x i) * lo i) ≤ ∫ t in x 0..x n, f t ∧
    (∫ t in x 0..x n, f t) ≤ ∑ i ∈ Finset.range n, (x (i + 1) - x i) * hi i := by
  rw [← intervalIntegral.sum_integral_adjacent_intervals hf]
  constructor
  · apply Finset.sum_le_sum
    intro i hi'
    have hin := Finset.mem_range.mp hi'
    simpa using intervalIntegral.integral_mono_on (hx (Nat.le_succ i))
      (intervalIntegrable_const (c := lo i)) (hf i hin) (fun t ht => (hbox i hin t ht).1)
  · apply Finset.sum_le_sum
    intro i hi'
    have hin := Finset.mem_range.mp hi'
    simpa using intervalIntegral.integral_mono_on (hx (Nat.le_succ i))
      (hf i hin) (intervalIntegrable_const (c := hi i)) (fun t ht => (hbox i hin t ht).2)

/-- At depth zero enclose the entire range. At each higher depth bisect the cell. -/
def rectangles? (range : ℚ → ℚ → Option RationalInterval) (a b : ℚ) :
    Nat → Option RationalInterval
  | 0 => (range a b).map (·.scaleNonnegative (b - a))
  | n + 1 => do
    let mid := (a + b) / 2
    let left ← rectangles? range a mid n
    let right ← rectangles? range mid b n
    pure (left.add right)

theorem rectangles_sound {f : ℝ → ℝ} {range : ℚ → ℚ → Option RationalInterval}
    (hbox : ∀ a b I, range a b = some I → ∀ t ∈ Icc (a : ℝ) (b : ℝ), I.Contains (f t))
    {n : Nat} {a b : ℚ} {I : RationalInterval} (hab : a ≤ b)
    (hf : IntervalIntegrable f volume a b) (h : rectangles? range a b n = some I) :
    I.Contains (∫ t in (a : ℝ)..(b : ℝ), f t) := by
  induction n generalizing a b I with
  | zero =>
    obtain ⟨J, hj, rfl⟩ := Option.map_eq_some_iff.mp h
    have hab' : (a : ℝ) ≤ b := by exact_mod_cast hab
    obtain ⟨hl, hu⟩ := And.intro
      (intervalIntegral.integral_mono_on hab' (intervalIntegrable_const (c := (J.lo : ℝ)))
        hf (fun t ht => (hbox _ _ _ hj t ht).1))
      (intervalIntegral.integral_mono_on hab' hf (intervalIntegrable_const (c := (J.hi : ℝ)))
        (fun t ht => (hbox _ _ _ hj t ht).2))
    simpa [RationalInterval.Contains, RationalInterval.scaleNonnegative_def] using And.intro hl hu
  | succ n ih =>
    let m := (a + b) / 2
    have ham : a ≤ m := by dsimp [m]; linarith
    have hmb : m ≤ b := by dsimp [m]; linarith
    have ham' : (a : ℝ) ≤ m := by exact_mod_cast ham
    have hmb' : (m : ℝ) ≤ b := by exact_mod_cast hmb
    have hfa : IntervalIntegrable f volume a m := hf.mono_set (by
      simpa [uIcc_of_le ham', uIcc_of_le (ham'.trans hmb')] using Icc_subset_Icc le_rfl hmb')
    have hfb : IntervalIntegrable f volume m b := hf.mono_set (by
      simpa [uIcc_of_le hmb', uIcc_of_le (ham'.trans hmb')] using Icc_subset_Icc ham' le_rfl)
    obtain ⟨J, hj, hrest⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨K, hk, hI⟩ := Option.bind_eq_some_iff.mp hrest
    cases hI
    rw [← intervalIntegral.integral_add_adjacent_intervals hfa hfb]
    exact RationalInterval.contains_add (ih ham hfa hj) (ih hmb hfb hk)

/-- Evaluate over real variable 0; FloatLib checks the complete input interval. -/
def exprRange (e : Expr) (bits degree : Nat) (a b : ℚ) : Option RationalInterval :=
  e.gridEnclosureWith? [⟨a, b⟩] bits degree

theorem exprRange_sound (e : Expr) (bits degree : Nat) (a b : ℚ) (I : RationalInterval)
    (h : exprRange e bits degree a b = some I) (t : ℝ) (ht : t ∈ Icc (a : ℝ) (b : ℝ)) :
    I.Contains (e.eval fun _ => t) := by
  apply Expr.contains_gridEnclosureWith? (h := h)
  intro i J hJ
  cases i with
  | zero => simp at hJ; subst J; exact ht
  | succ i => simp at hJ

/-- Integrate over 2^depth panels. An undefined range, reversed domain or missed
width budget returns none. Precision is in binary fractional bits. -/
def integrate? (e : Expr) (a b : ℚ) (depth bits : Nat) (tol : ℚ) (degree : Nat := 24) :
    Option RationalInterval :=
  if a ≤ b then (rectangles? (exprRange e bits degree) a b depth).bind (within? tol)
  else none

theorem integrate_sound {e : Expr} {a b tol : ℚ} {depth bits degree : Nat} {I : RationalInterval}
    (hf : IntervalIntegrable (fun t => e.eval fun _ => t) volume a b)
    (h : integrate? e a b depth bits tol degree = some I) :
    I.Contains (∫ t in (a : ℝ)..(b : ℝ), e.eval fun _ => t) ∧ I.hi - I.lo ≤ tol := by
  unfold integrate? at h
  split at h
  · obtain ⟨J, hj, hJ⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨rfl, hw⟩ := within_sound hJ
    exact ⟨rectangles_sound (exprRange_sound e bits degree) ‹a ≤ b› hf hj, hw⟩
  · contradiction

#print axioms integral_grid_bounds
#print axioms integrate_sound

end SashimiNumerics
