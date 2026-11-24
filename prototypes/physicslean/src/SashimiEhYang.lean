import SashimiEhFullGrid
import SashimiEpsKernel
import SashimiEhSweep

/-! The two variances in the Yang density use the same spectrum normalization.
Subtract their raw full integrals before dividing by that common positive anchor. -/

namespace SashimiEhYang

open SashimiEisensteinIntegral SashimiEhFull SashimiEhFullGrid FloatLib.Numerics

/-- Extend the child, parent and radius-8 raw integrals with their individual
tail allowances and propagate their shared normalization into the variance difference. -/
def varianceDifference (rc rp : ℚ) (child parent anchor : RationalInterval) : RationalInterval :=
  (((pad child (allowance rc)).sub (pad parent (allowance rp))).scaleNonnegative
    (matched.sigma8 ^ 2)).divPositive (pad anchor (allowance 8))

theorem varianceDifference_sound {mc mp rc rp : ℚ} {child parent anchor : RationalInterval}
    (hc : (rc:ℝ) ≤ radius matched mc) (hp : (rp:ℝ) ≤ radius matched mp)
    (hcRaw : child.Contains (raw matched lowEdge highEdge (radius matched mc)))
    (hpRaw : parent.Contains (raw matched lowEdge highEdge (radius matched mp)))
    (haRaw : anchor.Contains (raw matched lowEdge highEdge 8))
    (h : 0 < rc ∧ 1 ≤ highEdge * rc ∧ 0 < rp ∧ 1 ≤ highEdge * rp ∧ 0 < anchor.lo) :
    (varianceDifference rc rp child parent anchor).Contains (S mc - S mp) := by
  obtain ⟨hcr,hck,hpr,hpk,hpos⟩ := h
  have c := fullRaw_sound hc hcr hck hcRaw
  have p := fullRaw_sound hp hpr hpk hpRaw
  have a := fullRaw_sound (rmin := 8) (by norm_num) (by norm_num)
    (by norm_num [highEdge]) haRaw
  have result := RationalInterval.contains_divPositive
    (RationalInterval.contains_scaleNonnegative (RationalInterval.contains_sub c p)
      (sq_nonneg matched.sigma8)) a hpos
  simpa only [varianceDifference, S, SashimiSpectrum.normalize8, Rat.cast_pow,
    mul_sub, sub_div] using result

/-- Apply the guarded EPS enclosure to a previously certified positive difference. -/
theorem density_sound {d1 d0 ds out : RationalInterval} {b1 b0 s : ℝ} {bits degree : Nat}
    (hc : d1.Contains b1) (hp : d0.Contains b0) (hv : ds.Contains s)
    (hh : EpsKernel.enclose? d1 d0 ds (.point 0) bits degree = some out) :
    out.Contains (EpsKernel.density (b1-b0) s) := by
  simpa using (EpsKernel.sound hc hp hv
    (show (RationalInterval.point 0).Contains (0:ℝ) by norm_num [RationalInterval.Contains, RationalInterval.point]) hh).2.2

#print axioms varianceDifference_sound
#print axioms density_sound

open SashimiEhSweep SashimiSpectrum MeasureTheory Set

private def window (r : RationalInterval) (l u : ℚ) (bits degree : Nat) : RationalInterval :=
  SashimiTableWindow.windowSq (l*r.lo) (u*r.hi) bits degree

private theorem window_sound {r : RationalInterval} {l u : ℚ} {x k : ℝ} {bits degree : Nat}
    (hl : 0 < l) (hr : 0 ≤ r.lo) (hx : r.Contains x) (hk : k ∈ Icc (l:ℝ) u) :
    (window r l u bits degree).Contains ((topHat (k*x))^2) := by
  have hr' : 0 ≤ (r.lo : ℝ) := by exact_mod_cast hr
  have hl' : 0 ≤ (l : ℝ) := by exact_mod_cast hl.le
  apply SashimiTableWindow.windowSq_sound (mul_nonneg hl.le hr)
  · push_cast; exact mul_le_mul hk.1 hx.1 hr' (hl'.trans hk.1)
  · push_cast; exact mul_le_mul hk.2 hx.2 (hr'.trans hx.1) ((hl'.trans hk.1).trans hk.2)

/-- Reuse one spectrum enclosure while subdividing the two windows in k.
The signed difference is formed before multiplication by their common power. -/
def pairedPanel (c p power : RationalInterval) (bits degree : Nat) :
    Nat → ℚ → ℚ → RationalInterval
  | 0, l, u =>
    let dw := (window c l u bits degree).sub (window p l u bits degree)
    rounded bits ((power.mul dw).scaleNonnegative (u-l))
  | d+1, l, u =>
    let mid := (l+u)/2
    (pairedPanel c p power bits degree d l mid).add (pairedPanel c p power bits degree d mid u)

def pairedTotal (c p : RationalInterval) (bits degree depth : Nat) :=
  totalBy (fun l b => pairedPanel c p b.power bits degree depth l b.right)

/-- The physical integrand whose signed quadrature is enclosed by paired panels. -/
private noncomputable def profile (r k : ℝ) := weightReal matched k * (topHat (k*r))^2

private theorem measurableProfile (r : ℝ) : Measurable (profile r) :=
  SashimiAnalyticWindow.measurableIntegrand (weighted matched)
    (fun i => ((coefficients matched).toList[i]?).getD 0) r

private theorem panel_sound {c p power : RationalInterval} {bits degree : Nat} {x y : ℝ}
    (hc : c.Contains x) (hp : p.Contains y) (hc0 : 0 ≤ c.lo) (hp0 : 0 ≤ p.lo)
    (depth : Nat) {l u : ℚ} (hl : 0 < l) (hlu : l ≤ u)
    (hw : ∀ k ∈ Icc (l:ℝ) u, power.Contains (weightReal matched k)) :
    IntervalIntegrable (fun k => profile x k - profile y k) volume l u ∧
    (pairedPanel c p power bits degree depth l u).Contains
      (∫ k in (l:ℝ)..u, profile x k - profile y k) := by
  induction depth generalizing l u with
  | zero =>
    have I := cellIntegral_sound (I := power.mul ((window c l u bits degree).sub (window p l u bits degree)))
      ((measurableProfile x).sub (measurableProfile y)) hlu (by
        intro k hk
        simpa only [Pi.sub_apply,profile,mul_sub] using (RationalInterval.contains_mul (hw k hk)
          (RationalInterval.contains_sub (window_sound hl hc0 hc hk) (window_sound hl hp0 hp hk))))
    exact ⟨I.1,round_sound _ I.2⟩
  | succ d ih =>
    let mid : ℚ := (l+u)/2
    have lm : l ≤ mid := by dsimp [mid]; linarith
    have mu : mid ≤ u := by dsimp [mid]; linarith
    obtain ⟨hA,A⟩ := ih hl lm (fun k hk => hw k ⟨hk.1,hk.2.trans (by exact_mod_cast mu)⟩)
    obtain ⟨hB,B⟩ := ih (hl.trans_le lm) mu (fun k hk => hw k ⟨le_trans (by exact_mod_cast lm) hk.1,hk.2⟩)
    exact ⟨hA.trans hB,(intervalIntegral.integral_add_adjacent_intervals hA hB) ▸ RationalInterval.contains_add A B⟩

theorem paired_append (c p : RationalInterval) (bits degree depth : Nat) (a : ℚ) (xs ys : List Bin) :
    pairedTotal c p bits degree depth a (xs++ys) =
      (pairedTotal c p bits degree depth a xs).add (pairedTotal c p bits degree depth (last a xs) ys) := by
  induction xs generalizing a with
  | nil => simp [pairedTotal,totalBy,last,RationalInterval.point,RationalInterval.add_def]
  | cons x xs ih =>
    change (pairedPanel c p x.power bits degree depth a x.right).add (pairedTotal c p bits degree depth x.right (xs++ys)) = _
    rw [ih]
    change RationalInterval.add _ (RationalInterval.add _ _) = RationalInterval.add (RationalInterval.add _ _) _
    simp [RationalInterval.add_def,add_assoc,last,pairedTotal]

/-- For any actual radii within the boxes, the paired sum bounds the difference
of the two finite-band continuum EH integrals. -/
theorem pairedTotal_sound {a b : ℚ} {bins : List Bin} {bits degree depth : Nat} {c p : RationalInterval} {x y : ℝ}
    (good : Good (weightReal matched) a bins) (hEnd : last a bins = b)
    (hLo : a = lowEdge) (hHi : b = highEdge)
    (hc : c.Contains x) (hp : p.Contains y) (hc0 : 0 ≤ c.lo) (hp0 : 0 ≤ p.lo) :
    (pairedTotal c p bits degree depth a bins).Contains (raw matched a b x - raw matched a b y) := by
  have h := (totalBy_sound good (fun l q hl hu hw =>
    panel_sound (bits := bits) (degree := degree) hc hp hc0 hp0 depth hl hu hw)).2
  rw [hEnd] at h
  have hx := matched_integrable x ((by exact_mod_cast hc0 : (0:ℝ) ≤ c.lo).trans hc.1)
  have hy := matched_integrable y ((by exact_mod_cast hp0 : (0:ℝ) ≤ p.lo).trans hp.1)
  rw [← hLo, ← hHi] at hx hy
  change IntervalIntegrable (profile x) volume a b at hx
  change IntervalIntegrable (profile y) volume a b at hy
  rw [intervalIntegral.integral_sub hx hy] at h
  exact h

/-- The finite difference can have either sign. The parent tail contributes to
the lower end and the child tail to the upper end. The common anchor stays positive. -/
def pairedDifference (rc rp : ℚ) (J anchor : RationalInterval) : RationalInterval :=
  ((⟨J.lo-allowance rp, J.hi+allowance rc⟩ : RationalInterval).scaleNonnegative
    (matched.sigma8^2)).divPositive (pad anchor (allowance 8))

theorem pairedDifference_sound {rc rp : ℚ} {x y : ℝ} {J a : RationalInterval}
    (hrc : (rc:ℝ) ≤ x) (hrp : (rp:ℝ) ≤ y)
    (hJ : J.Contains (raw matched lowEdge highEdge x - raw matched lowEdge highEdge y))
    (ha : a.Contains (raw matched lowEdge highEdge 8))
    (h : 0 < rc ∧ 1 ≤ highEdge * rc ∧ 0 < rp ∧ 1 ≤ highEdge * rp ∧ 0 < a.lo) :
    (pairedDifference rc rp J a).Contains
      (normalize8 matched.sigma8 (fullRaw matched 8) (fullRaw matched x) -
       normalize8 matched.sigma8 (fullRaw matched 8) (fullRaw matched y)) := by
  obtain ⟨hcr,hck,hpr,hpk,hpos⟩ := h
  have C := fullRaw_tail hrc hcr hck
  have P := fullRaw_tail hrp hpr hpk
  have A := fullRaw_sound (rmin := 8) (by norm_num) (by norm_num) (by norm_num [highEdge]) ha
  have hj : (⟨J.lo-allowance rp,J.hi+allowance rc⟩ : RationalInterval).Contains
      (fullRaw matched x - fullRaw matched y) := by
    constructor <;> dsimp only <;> push_cast <;> linarith [hJ.1,hJ.2,C.1,C.2,P.1,P.2]
  have result := RationalInterval.contains_divPositive
    (RationalInterval.contains_scaleNonnegative hj (sq_nonneg matched.sigma8)) A hpos
  simpa only [pairedDifference,normalize8,Rat.cast_pow,mul_sub,sub_div] using result

#print axioms pairedTotal_sound
#print axioms pairedDifference_sound

/-- Check the Yang density on consecutive subdivisions of the variance difference.
Each subinterval uses one instance of the same difference in both the exponential
and the inverse power. Their hull covers the whole input. -/
def densitySplit? (b1 b0 : RationalInterval) (bits degree : Nat) :
    Nat → RationalInterval → Option RationalInterval
  | 0, d => EpsKernel.enclose? b1 b0 d (.point 0) bits degree
  | n+1, d => do
    let m := (d.lo + d.hi)/2
    let l ← densitySplit? b1 b0 bits degree n ⟨d.lo,m⟩
    let r ← densitySplit? b1 b0 bits degree n ⟨m,d.hi⟩
    return ⟨min l.lo r.lo,max l.hi r.hi⟩

theorem densitySplit_sound {b1 b0 : RationalInterval} {x y s : ℝ} {bits degree depth : Nat}
    {d I : RationalInterval} (h1 : b1.Contains x) (h0 : b0.Contains y) (hd : d.Contains s)
    (done_ : densitySplit? b1 b0 bits degree depth d = some I) :
    I.Contains (EpsKernel.density (x-y) s) := by
  induction depth generalizing d I with
  | zero => exact density_sound h1 h0 hd done_
  | succ n ih =>
    let m := (d.lo+d.hi)/2
    obtain ⟨l,hl,h⟩ := Option.bind_eq_some_iff.mp done_
    obtain ⟨r,hr,e⟩ := Option.bind_eq_some_iff.mp h
    cases e
    by_cases hm : s ≤ (m:ℝ)
    · have hv := ih (d := ⟨d.lo,m⟩) (I := l) ⟨hd.1,hm⟩ hl
      constructor
      · change ((min l.lo r.lo : ℚ):ℝ) ≤ _
        rw [Rat.cast_min]
        exact (min_le_left ..).trans hv.1
      · change _ ≤ ((max l.hi r.hi : ℚ):ℝ)
        rw [Rat.cast_max]
        exact hv.2.trans (le_max_left ..)
    · have hv := ih (d := ⟨m,d.hi⟩) (I := r) ⟨le_of_not_ge hm,hd.2⟩ hr
      constructor
      · change ((min l.lo r.lo : ℚ):ℝ) ≤ _
        rw [Rat.cast_min]
        exact (min_le_right ..).trans hv.1
      · change _ ≤ ((max l.hi r.hi : ℚ):ℝ)
        rw [Rat.cast_max]
        exact hv.2.trans (le_max_right ..)

#print axioms densitySplit_sound
end SashimiEhYang
