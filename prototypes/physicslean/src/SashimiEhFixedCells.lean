import SashimiEhFixedWeight
import SashimiTighterWindow

/-! Variance and signed derivative on a common fixed mesh. -/
namespace SashimiEhFixedCells
open FloatLib.Numerics MeasureTheory Set SashimiEhSweep SashimiEisensteinIntegral

def total (window : ℚ → ℚ → RationalInterval) (bits : Nat) (a : ℚ) : List Cell → RationalInterval
  | [] => .point 0
  | q::qs => (rounded bits (q.integral.mul (window a q.right))).add (total window bits q.right qs)

/-- A nonnegative density allows a signed window. Its cell integral and its
pointwise range play separate roles in the proof. -/
theorem total_sound {w g : ℝ → ℝ} {window : ℚ → ℚ → RationalInterval} {a qs bits}
    (hw : Measurable w) (hfg : Measurable (fun x => w x * g x)) (hn : ∀ x, 0 ≤ w x)
    (hg : Good w a qs)
    (hW : ∀ (l u : ℚ), 0 < l → l ≤ u → ∀ x ∈ Icc (l:ℝ) u, (window l u).Contains (g x)) :
    IntervalIntegrable (fun x => w x * g x) volume a (lastEdge a qs) ∧
      (total window bits a qs).Contains (∫ k in (a:ℝ)..(lastEdge a qs:ℝ), w k * g k) := by
  induction qs generalizing a with
  | nil => simp [lastEdge,total,RationalInterval.point,RationalInterval.Contains]
  | cons q qs ih =>
    obtain ⟨hv,ht⟩ := hg
    have hJ := hW a q.right hv.pos hv.ordered
    have hi := (cellIntegral_sound hw hv.ordered hv.range).1
    have hp := (cellIntegral_sound (I := q.pointwise.mul (window a q.right)) hfg hv.ordered
      (fun x hx => RationalInterval.contains_mul (hv.range x hx) (hJ x hx))).1
    have hA := CertifiedWeightedIntegral.weighted_sound (by exact_mod_cast hv.ordered) hi hp
      (fun x _ => hn x) hv.mass hJ
    obtain ⟨hj,hB⟩ := ih ht
    exact ⟨hp.trans hj, (intervalIntegral.integral_add_adjacent_intervals hp hj) ▸
      RationalInterval.contains_add (round_sound bits hA) hB⟩

private theorem argument {l u : ℚ} {R : RationalInterval} {k r : ℝ}
    (hl : 0 < l) (hR : R.Contains r) (h0 : 0 ≤ R.lo) (hk : k ∈ Icc (l:ℝ) u) :
    k*r ∈ Icc ((l*R.lo:ℚ):ℝ) ((u*R.hi:ℚ):ℝ) := by
  have hrp : 0 ≤ (R.lo:ℝ) := by exact_mod_cast h0
  have hlp : 0 ≤ (l:ℝ) := by exact_mod_cast hl.le
  constructor <;> push_cast
  · exact mul_le_mul hk.1 hR.1 hrp (hlp.trans hk.1)
  · exact mul_le_mul hk.2 hR.2 (hrp.trans hR.1) ((hlp.trans hk.1).trans hk.2)

private theorem weight_nonneg (c : Cosmology) (x : ℝ) : 0 ≤ weightReal c x := by
  unfold weightReal
  positivity

def variance (R : RationalInterval) (bits degree : Nat) :=
  total (fun l u => SashimiTighterWindow.windowSq (l*R.lo) (u*R.hi) bits degree) bits

def slope (R : RationalInterval) (bits degree : Nat) :=
  total (fun l u => SashimiTighterWindow.weight (l*R.lo) (u*R.hi) bits degree) bits

theorem variance_sound {c a qs R bits degree r} (hg : Good (weightReal c) a qs)
    (hR : R.Contains r) (hp : 0 < R.lo) :
    (variance R bits degree a qs).Contains (raw c a (lastEdge a qs) r) := by
  have hm := SashimiAnalyticWindow.measurableIntegrand (weighted c)
    (fun i => ((coefficients c).toList[i]?).getD 0) r
  apply (total_sound (g := fun k => SashimiSpectrum.topHat (k*r)^2)
    (SashimiEhDerivativeBounds.measurable_weight c) hm (weight_nonneg c) hg ?_).2
  intro l u hl _ k hk
  exact SashimiTighterWindow.windowSq_sound (mul_nonneg hl.le hp.le) (argument hl hR hp.le hk).1
    (argument hl hR hp.le hk).2 bits degree

theorem slope_sound {c a qs R bits degree r} (hg : Good (weightReal c) a qs)
    (hR : R.Contains r) (hp : 0 < R.lo) :
    (slope R bits degree a qs).Contains (SashimiEhDerivativeBounds.rawSlope c r a (lastEdge a qs)) := by
  have hw := SashimiEhDerivativeBounds.measurable_weight c
  have hm := hw.mul (SashimiEhDerivativeBounds.measurable_shape.comp (measurable_id.mul_const r))
  apply (total_sound (g := fun k => SashimiEhDerivativeBounds.shape (k*r))
    hw hm (weight_nonneg c) hg ?_).2
  intro l u hl _ k hk
  exact SashimiTighterWindow.weight_sound (mul_pos hl hp) (argument hl hR hp.le hk) bits degree

/-- The same denominator and mass radius as the existing finite-band real model. -/
theorem normalized_variance {a : ℚ} {qs R I m tol bits degree} (ha : 0 < a)
    (hab : a ≤ lastEdge a qs) (hg : Good (weightReal matched) a qs)
    (hr : radiusBox? matched m bits degree = some R) (hp : 0 < R.lo)
    (h : SashimiNormalization.apply? matched.sigma8 tol (variance R bits degree a qs)
      (variance (.point 8) bits degree a qs) = some I) :
    I.Contains (SashimiEhVariance.normalizedS m a (lastEdge a qs)) ∧ I.hi-I.lo ≤ tol := by
  rw [← SashimiEhSweep.normalized_eq ha hab m]
  apply SashimiNormalization.apply_sound (variance_sound hg (radiusBox_sound hr) hp) _ h
  apply variance_sound hg (R := RationalInterval.point 8)
  · norm_num [RationalInterval.Contains,RationalInterval.point]
  · decide

theorem normalized_slope {a : ℚ} {qs R I m tol bits degree} (ha : 0 < a) (hm : 0 < m)
    (hab : a ≤ lastEdge a qs) (hg : Good (weightReal matched) a qs)
    (hr : radiusBox? matched m bits degree = some R) (hp : 0 < R.lo)
    (h : SashimiNormalization.apply? matched.sigma8 tol ((slope R bits degree a qs).scaleNonnegative (1/m))
      (variance (.point 8) bits degree a qs) = some I) :
    I.Contains (deriv (fun t => SashimiEhVariance.normalizedS t a (lastEdge a qs)) m) ∧
      I.hi-I.lo ≤ tol := by
  have hd := RationalInterval.contains_scaleNonnegative (slope_sound (bits := bits) (degree := degree) hg (radiusBox_sound hr) hp)
    (one_div_nonneg.mpr hm.le)
  have hA := variance_sound (bits := bits) (degree := degree) hg (R := RationalInterval.point 8)
    (r := 8) (by norm_num [RationalInterval.Contains,RationalInterval.point]) (by decide)
  obtain ⟨hI,hwidth⟩ := SashimiNormalization.apply_sound hd hA h
  push_cast at hI
  rw [SashimiEhDerivativeBounds.normalized_eq hm ha hab] at hI
  rw [(SashimiEhContinuity.matched_derivative (by exact_mod_cast hm)
    (by exact_mod_cast ha) (by exact_mod_cast hab)).deriv]
  exact ⟨hI,hwidth⟩

structure MassResult where
  radius : RationalInterval
  numerator : RationalInterval
  variance : RationalInterval
  derivative : RationalInterval
deriving Repr, DecidableEq

/-- Reuse the same finite-band anchor and cell integrals for each mass. -/
def mass? (m tol a : ℚ) (cells : List Cell) (anchor : RationalInterval) (bits degree : Nat) :
    Option MassResult := do
  let R ← radiusBox? matched m bits degree
  guard (0 < m ∧ 0 < R.lo)
  let N := variance R bits degree a cells
  let S ← SashimiNormalization.apply? matched.sigma8 tol N anchor
  let D ← SashimiNormalization.apply? matched.sigma8 tol
    ((slope R bits degree a cells).scaleNonnegative (1/m)) anchor
  return ⟨R,N,S,D⟩

/-- Standalone finite-band evaluator on arbitrary prescribed rational edges. -/
def evaluate? (m tol a : ℚ) (ends : List ℚ) (bits degree : Nat) : Option MassResult := do
  let boxes ← boxes? matched bits degree
  let cells ← build? matched boxes bits degree a ends
  mass? m tol a cells (variance (.point 8) bits degree a cells) bits degree

theorem evaluate_sound {m tol a ends bits degree result}
    (ha : 0 < a) (hab : a ≤ ends.getLastD a)
    (h : evaluate? m tol a ends bits degree = some result) :
    result.variance.Contains (SashimiEhVariance.normalizedS m a (ends.getLastD a)) ∧
    result.derivative.Contains (deriv (fun t => SashimiEhVariance.normalizedS t a (ends.getLastD a)) m) ∧
    result.variance.hi-result.variance.lo ≤ tol ∧ result.derivative.hi-result.derivative.lo ≤ tol := by
  obtain ⟨C,hC,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨cs,hcs,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨hg,he⟩ := build_sound hC hcs
  rw [← he] at hab ⊢
  unfold mass? at h
  obtain ⟨R,hr,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨_,hpos,h⟩ := Option.bind_eq_some_iff.mp h
  have hp : 0 < m ∧ 0 < R.lo := by by_contra hn; simp [guard,hn] at hpos
  obtain ⟨S,hS,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨D,hD,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  obtain ⟨hV,hvw⟩ := normalized_variance ha hab hg hr hp.2 hS
  obtain ⟨hL,hlw⟩ := normalized_slope ha hp.1 hab hg hr hp.2 hD
  exact ⟨hV,hL,hvw,hlw⟩

#print axioms evaluate_sound
#print axioms normalized_variance
#print axioms normalized_slope
#print axioms variance_sound
#print axioms slope_sound
end SashimiEhFixedCells
