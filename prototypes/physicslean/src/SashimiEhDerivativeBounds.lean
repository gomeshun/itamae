import SashimiEhSweep
import SashimiEhContinuity

/-! Whole-cell bounds for the signed physical-mass derivative of the spherical
window. The power and the σ₈ anchor can be shared with the variance sweep. -/
namespace SashimiEhDerivativeBounds
open Real Set MeasureTheory FloatLib.Numerics FloatLib.Numerics.Interval SashimiSpectrum

/-- The logarithmic-mass window derivative, before division by physical mass. -/
noncomputable def shape (x : ℝ) := 2 * (topHat x * (sin x / x - topHat x))

theorem shape_eq {x : ℝ} (hx : 0 < x) :
    shape x = 2 * topHat x * topHatPrime x * x / 3 := by
  simp only [shape, topHat, topHatPrime, ite_eq_right hx.ne']
  field_simp
  ring

/-- Intersect an analytic estimate with FloatLib trigonometric evaluation in
the oscillatory middle band. The analytic range works at either extreme. -/
private def refine (I : RationalInterval) (e : Expr) (l u : ℚ) (bits degree : Nat) :
    RationalInterval :=
  if 1/4 ≤ l ∧ u ≤ 48 then
    match SashimiNumerics.exprRange e bits degree l u with
    | some J => ⟨max I.lo J.lo, min I.hi J.hi⟩
    | none => I
  else I

private theorem refine_sound {I : RationalInterval} {y x : ℝ} (h : I.Contains y)
    {e : Expr} {l u : ℚ} (hx : x ∈ Icc (l:ℝ) u)
    (he : e.eval (fun _ => x) = y) (bits degree : Nat) :
    (refine I e l u bits degree).Contains y := by
  unfold refine
  split
  · split
    · rename_i J hj
      have hh := SashimiNumerics.exprRange_sound _ _ _ _ _ _ hj x hx
      rw [he] at hh
      constructor <;> push_cast
      · exact max_le h.1 hh.1
      · exact le_min h.2 hh.2
    · exact h
  · exact h

private def windowBox (l u : ℚ) : RationalInterval :=
  let bound := min 1 (3*(1+u)/l^3)
  ⟨max (-bound) (1-u^2/10), bound⟩

private def sincBox (l u : ℚ) : RationalInterval :=
  let bound := min 1 (1/l)
  ⟨max (-bound) (1-u^2/6), bound⟩

private theorem windowBox_sound {l u : ℚ} {x : ℝ} (hl : 0 < l)
    (hx : x ∈ Icc (l:ℝ) u) : (windowBox l u).Contains (topHat x) := by
  have hl' : 0 < (l:ℝ) := by exact_mod_cast hl
  have hp := hl'.trans_le hx.1
  have hu := hp.le.trans hx.2
  obtain ⟨hlo, hhi⟩ := topHat_bounds x hp.le
  have hdec : |topHat x| ≤ 3*(1+(u:ℝ))/(l:ℝ)^3 := (topHat_decay x hp).trans (by gcongr <;> first | exact hx.1 | exact hx.2)
  have hb : |topHat x| ≤ min 1 (3*(1+(u:ℝ))/(l:ℝ)^3) := le_min (abs_le.mpr ⟨hlo,hhi⟩) hdec
  have ht : 1-(u:ℝ)^2/10 ≤ topHat x := by
    have h2 := pow_le_pow_left₀ hp.le hx.2 2
    linarith [topHat_taylor_lower x hp.le]
  unfold windowBox
  constructor <;> push_cast
  · exact max_le (abs_le.mp hb).1 ht
  · exact (abs_le.mp hb).2

private theorem sincBox_sound {l u : ℚ} {x : ℝ} (hl : 0 < l)
    (hx : x ∈ Icc (l:ℝ) u) : (sincBox l u).Contains (sin x / x) := by
  have hl' : 0 < (l:ℝ) := by exact_mod_cast hl
  have hp := hl'.trans_le hx.1
  have habs : |sin x / x| ≤ min 1 (1/(l:ℝ)) := by
    rw [abs_div, abs_of_pos hp]
    apply le_min
    · apply (div_le_one hp).2
      simpa [abs_of_pos hp] using (Real.abs_sin_le_abs : |sin x| ≤ |x|)
    · exact (div_le_div_of_nonneg_right (Real.abs_sin_le_one x) hp.le).trans
        (one_div_le_one_div_of_le hl' hx.1)
  have htx : 1-x^2/6 ≤ sin x / x := (le_div_iff₀ hp).2 (by nlinarith [Real.sin_ge_sub_cube hp.le])
  have htu : 1-(u:ℝ)^2/6 ≤ sin x / x := by
    have h2 := pow_le_pow_left₀ hp.le hx.2 2
    linarith
  unfold sincBox
  constructor <;> push_cast
  · exact max_le (abs_le.mp habs).1 htu
  · exact (abs_le.mp habs).2

/-- Enclose 2 W(x) (sinc(x) − W(x)) for every x in a positive rational box. -/
def weight (l u : ℚ) (bits degree : Nat) : RationalInterval :=
  let W := refine (windowBox l u) (topHatExpr (.var 0)) l u bits degree
  let S := refine (sincBox l u) (.binary .div (.unary .sin (.var 0)) (.var 0)) l u bits degree
  (W.mul (S.sub W)).scaleNonnegative 2

theorem weight_sound {l u : ℚ} {x : ℝ} (hl : 0 < l) (hx : x ∈ Icc (l:ℝ) u) (bits degree : Nat) :
    (weight l u bits degree).Contains (shape x) := by
  have hp : 0 < x := (by exact_mod_cast hl : 0 < (l:ℝ)).trans_le hx.1
  have hw := refine_sound (windowBox_sound hl hx) hx
    (eval_topHatExpr (show (Expr.var 0).eval (fun _ => x) ≠ 0 from hp.ne')) bits degree
  have hs := refine_sound (sincBox_sound hl hx) hx
    (e := Expr.binary .div (.unary .sin (.var 0)) (.var 0)) rfl bits degree
  exact_mod_cast RationalInterval.contains_scaleNonnegative
    (RationalInterval.contains_mul hw (RationalInterval.contains_sub hs hw)) (by norm_num : (0:ℚ)≤2)


open SashimiEhSweep SashimiEisensteinIntegral SashimiEhVariance

/-- The integral prior to division by M and normalization against R=8. -/
noncomputable def rawSlope (c : Cosmology) (r : ℝ) (a b : ℚ) :=
  ∫ k in (a:ℝ)..(b:ℝ), weightReal c k * shape (k*r)

private theorem measurable_weight (c : Cosmology) : Measurable (weightReal c) := by
  have hh := SashimiAnalyticWindow.measurableIntegrand (weighted c)
    (fun i => ((coefficients c).toList[i]?).getD 0) 0
  change Measurable (fun k => weightReal c k * topHat (k*0)^2) at hh
  simpa [topHat] using hh

private theorem measurable_shape : Measurable shape := by
  have hw : Measurable topHat := by
    unfold topHat
    exact Measurable.ite (measurableSet_singleton 0) measurable_const (by fun_prop)
  unfold shape
  fun_prop

def panel (R : RationalInterval) (bits degree : Nat) (l : ℚ) (q : Bin) : RationalInterval :=
  SashimiEhSweep.rounded bits
    ((q.power.mul (weight (l*R.lo) (q.right*R.hi) bits degree)).scaleNonnegative (q.right-l))

def sum (R : RationalInterval) (bits degree : Nat) := totalBy (panel R bits degree)

theorem sum_append (R : RationalInterval) (bits degree : Nat) (a : ℚ) (xs ys : List Bin) :
    sum R bits degree a (xs ++ ys) =
      (sum R bits degree a xs).add (sum R bits degree (last a xs) ys) := by
  induction xs generalizing a with
  | nil =>
    simp [sum, totalBy, last, RationalInterval.point, RationalInterval.add_def]
  | cons x xs ih =>
    change (panel R bits degree a x).add (sum R bits degree x.right (xs ++ ys)) = _
    rw [ih]
    change RationalInterval.add _ (RationalInterval.add _ _) =
      RationalInterval.add (RationalInterval.add _ _) _
    simp [RationalInterval.add_def, add_assoc, last, sum]

theorem sum_sound {c a bins R r bits degree} (hgood : Good (weightReal c) a bins)
    (hR : R.Contains r) (hpos : 0 < R.lo) :
    (sum R bits degree a bins).Contains (rawSlope c r a (last a bins)) := by
  let f := fun k => weightReal c k * shape (k*r)
  have hf : Measurable f := (measurable_weight c).mul (measurable_shape.comp (measurable_id.mul_const r))
  apply (totalBy_sound (f := f) (piece := panel R bits degree) hgood ?_).2
  intro l q hl hu hp
  have hi := cellIntegral_sound (I := q.power.mul (weight (l*R.lo) (q.right*R.hi) bits degree)) hf hu (by
    intro k hk
    apply RationalInterval.contains_mul (hp k hk)
    apply weight_sound (mul_pos hl hpos)
    have hrp : 0 ≤ (R.lo:ℝ) := by exact_mod_cast hpos.le
    have hlp : 0 ≤ (l:ℝ) := by exact_mod_cast hl.le
    constructor <;> push_cast
    · exact mul_le_mul hk.1 hR.1 hrp (hlp.trans hk.1)
    · exact mul_le_mul hk.2 hR.2 (hrp.trans hR.1) ((hlp.trans hk.1).trans hk.2))
  exact ⟨hi.1, SashimiEhSweep.round_sound _ hi.2⟩

private theorem ordered {w a bins} (h : Good w a bins) : a ≤ last a bins := by
  induction bins generalizing a with
  | nil => exact le_rfl
  | cons b more ih => exact h.2.1.trans (ih h.2.2.2)

theorem normalized_eq {m a b : ℚ} (hm : 0 < m) (ha : 0 < a) (hab : a ≤ b) :
    normalize8 matched.sigma8 (raw matched a b 8)
        ((1 / (m:ℝ)) * rawSlope matched (radius matched m) a b) =
      normalizedDerivative m a b := by
  have hm' : 0 < (m:ℝ) := by exact_mod_cast hm
  have hR : 0 < massRadius massDensity (m:ℝ) := by
    unfold massRadius massDensity
    positivity
  have hi : rawSlope matched (radius matched m) a b =
      ((2*π^2)*(m:ℝ)) * (∫ k in (a:ℝ)..(b:ℝ), SashimiMassVariance.kernel
          .topHatDerivative realPower massDensity m 1 k) := by
    unfold rawSlope
    rw [← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr
    intro k hk
    rw [uIcc_of_le (by exact_mod_cast hab : (a:ℝ) ≤ b)] at hk
    have hk' := (by exact_mod_cast ha : (0:ℝ) < a).trans_le hk.1
    dsimp only
    rw [radius_eq, shape_eq (mul_pos hk' hR)]
    unfold weightReal
    rw [tilt_eq matched hk', Real.rpow_add hk']
    have hc : coefficients matched = SashimiEisensteinHu.prepareReal 0.674 0.315 0.049 2.7255 := by
      norm_num [coefficients,matched]
    rw [hc]
    simp only [SashimiMassVariance.kernel,realPower,matched,Real.rpow_two]
    norm_num
    field_simp [hm'.ne']
  rw [hi, raw_eq_scaled ha hab]
  unfold normalize8 normalizedDerivative
  norm_num [matched]
  field_simp [hm'.ne']

/-- Reuse a certified power mesh and its R=8 enclosure. A successful result has
a determined sign and the requested absolute width, in inverse solar masses. -/
def derivative? (m tol : ℚ) (a : ℚ) (bins : List Bin) (anchor : RationalInterval)
    (bits degree : Nat) : Option RationalInterval := do
  let R ← radiusBox? matched m bits degree
  guard (0 < m ∧ 0 < a ∧ 0 < R.lo)
  let I ← SashimiNormalization.apply? matched.sigma8 tol
    ((sum R bits degree a bins).scaleNonnegative (1/m)) anchor
  guard (I.hi < 0 ∨ 0 < I.lo)
  return I

theorem derivative_sound {m tol a bins anchor bits degree I}
    (hg : Good (weightReal matched) a bins) (ha : anchor.Contains (raw matched a (last a bins) 8))
    (h : derivative? m tol a bins anchor bits degree = some I) :
    I.Contains (normalizedDerivative m a (last a bins)) ∧ I.hi-I.lo ≤ tol ∧ (I.hi<0 ∨ 0<I.lo) := by
  unfold derivative? at h
  obtain ⟨R,hR,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨_,hpos,h⟩ := Option.bind_eq_some_iff.mp h
  have hp : 0 < m ∧ 0 < a ∧ 0 < R.lo := by by_contra hh; simp [guard,hh] at hpos
  obtain ⟨J,hJ,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨_,hsg,hh⟩ := Option.bind_eq_some_iff.mp h
  cases hh
  have hs : I.hi<0 ∨ 0<I.lo := by by_contra hn; simp [guard,hn] at hsg
  have hsum := sum_sound (bits := bits) (degree := degree) hg (radiusBox_sound hR) hp.2.2
  have hv := RationalInterval.contains_scaleNonnegative hsum (one_div_nonneg.mpr hp.1.le)
  obtain ⟨hN,hw⟩ := SashimiNormalization.apply_sound hv ha hJ
  push_cast at hN
  rw [normalized_eq hp.1 hp.2.1 (ordered hg)] at hN
  exact ⟨hN,hw,hs⟩

/-- Compose a replayed signed sum with the shared checked EH spectrum and
denominator. The conclusion refers to the derivative of the variance itself. -/
theorem replayed_deriv {m a b tol : ℚ} {bins R A N I bits degree}
    (hm : 0 < m) (ha : 0 < a) (hab : a ≤ b)
    (hg : Good (weightReal matched) a bins) (he : last a bins = b)
    (hr : radiusBox? matched m bits degree = some R) (hp : 0 < R.lo)
    (hn : sum R bits degree a bins = N) (hd : A.Contains (raw matched a b 8))
    (hI : SashimiNormalization.apply? matched.sigma8 tol (N.scaleNonnegative (1 / m)) A = some I) :
    I.Contains (deriv (fun t => normalizedS t a b) m) := by
  have hsum := sum_sound (bits := bits) (degree := degree) hg (radiusBox_sound hr) hp
  rw [hn, he] at hsum
  have hv := RationalInterval.contains_scaleNonnegative hsum (one_div_nonneg.mpr hm.le)
  obtain ⟨hN,_⟩ := SashimiNormalization.apply_sound hv hd hI
  push_cast at hN
  rw [normalized_eq hm ha hab] at hN
  rw [(SashimiEhContinuity.matched_derivative (by exact_mod_cast hm)
    (by exact_mod_cast ha) (by exact_mod_cast hab)).deriv]
  exact hN

/-- Standalone checker for replay or one-off bands. In a sweep, the
transfer bins and anchor can instead be computed once for all masses. -/
def fromMesh? (m tol a : ℚ) (ends : List ℚ) (bits degree : Nat) : Option RationalInterval := do
  let C ← boxes? matched bits degree
  let bins ← build? matched C bits degree a ends
  let A := SashimiEhSweep.total (.point 8) bits degree a bins
  derivative? m tol a bins A bits degree

theorem fromMesh_sound {m tol a ends bits degree I}
    (h : fromMesh? m tol a ends bits degree = some I) :
    I.Contains (normalizedDerivative m a (ends.getLastD a)) ∧ I.hi-I.lo ≤ tol ∧ (I.hi<0 ∨ 0<I.lo) := by
  obtain ⟨C,hC,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨bins,hb,h⟩ := Option.bind_eq_some_iff.mp h
  have hg := build_sound hC hb
  have hA := total_sound (bits := bits) (degree := degree) hg
    (R := RationalInterval.point 8) (r := 8)
    (by norm_num [RationalInterval.Contains,RationalInterval.point])
    (by norm_num [RationalInterval.point])
  simpa only [build_last hb] using derivative_sound hg hA h

/-- An accepted standalone result encloses the mathematical derivative with
respect to physical mass of the normalized variance on this wavenumber band. -/
theorem checked_deriv {m tol a : ℚ} {ends : List ℚ} {bits degree I}
    (hm : 0 < m) (ha : 0 < a) (hab : a ≤ ends.getLastD a)
    (h : fromMesh? m tol a ends bits degree = some I) :
    I.Contains (deriv (fun M => normalizedS M a (ends.getLastD a)) (m : ℝ)) := by
  rw [(SashimiEhContinuity.matched_derivative
    (show 0 < (m : ℝ) by exact_mod_cast hm)
    (show 0 < (a : ℝ) by exact_mod_cast ha)
    (show (a : ℝ) ≤ (ends.getLastD a : ℚ) by exact_mod_cast hab)).deriv]
  exact (fromMesh_sound h).1

#print axioms replayed_deriv
#print axioms checked_deriv
#print axioms derivative_sound
#print axioms fromMesh_sound
end SashimiEhDerivativeBounds
