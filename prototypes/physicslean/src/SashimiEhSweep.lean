import SashimiEhVariance
import SashimiEisensteinIntegral

/-! Shared positive-k power enclosures for analytic EH98 on a finite band.
A cell's transfer range is certified once and reused for every true mass
radius. The quadrature is a sum of range-times-width, without interpolation. -/
namespace SashimiEhSweep
open SashimiEisensteinHu SashimiEisensteinIntegral SashimiEhVariance
open SashimiSpectrum FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory Set Real

structure Bin where
  right : ℚ
  power : RationalInterval
deriving Repr, DecidableEq

def tiltExpr (c : Cosmology) : Expr := .unary .exp (.const (2+c.ns) * .unary .log (.var 0))

def weight? (c : Cosmology) (boxes : List RationalInterval) (l u : ℚ) (bits degree : Nat) :
    Option RationalInterval := do
  guard (0 < l ∧ l ≤ u)
  let T ← rangeFast? boxes ⟨l,u⟩ bits degree
  let P ← SashimiNumerics.exprRange (tiltExpr c) bits degree l u
  -- The mathematical square is nonnegative even when its interval crosses zero.
  let Q := T.mul T
  return P.mul ⟨max 0 Q.lo,Q.hi⟩

def build? (c : Cosmology) (boxes : List RationalInterval) (bits degree : Nat) (a : ℚ) :
    List ℚ → Option (List Bin)
  | [] => some []
  | b::bs => do
    let P ← weight? c boxes a b bits degree
    return ⟨b,P⟩ :: (← build? c boxes bits degree b bs)

def rounded (bits : Nat) (I : RationalInterval) : RationalInterval :=
  let s : ℚ := 2^bits
  ⟨⌊I.lo*s⌋/s,⌈I.hi*s⌉/s⟩

def panel (R : RationalInterval) (bits degree : Nat) (a : ℚ) (q : Bin) : RationalInterval :=
  let W := SashimiTableWindow.windowSq (a*R.lo) (q.right*R.hi) bits degree
  rounded bits ((q.power.mul W).scaleNonnegative (q.right-a))

def totalBy (piece : ℚ → Bin → RationalInterval) (a : ℚ) : List Bin → RationalInterval
  | [] => .point 0
  | q::qs => (piece a q).add (totalBy piece q.right qs)

def total (R : RationalInterval) (bits degree : Nat) := totalBy (panel R bits degree)

def normalize? (c : Cosmology) (m tol : ℚ) (bins : List Bin) (a : ℚ) (A : RationalInterval)
    (bits degree : Nat) : Option RationalInterval := do
  let R ← radiusBox? c m bits degree
  guard (0 ≤ R.lo)
  SashimiNormalization.apply? c.sigma8 tol (total R bits degree a bins) A

noncomputable def weightReal (c : Cosmology) (k : ℝ) :=
  exp (((2+c.ns : ℚ):ℝ)*log k)*(transferReal (coefficients c) k)^2

def Good (f : ℝ → ℝ) (a : ℚ) : List Bin → Prop
  | [] => True
  | q::qs => 0<a ∧ a≤q.right ∧ (∀ k∈Icc (a:ℝ) q.right, q.power.Contains (f k)) ∧ Good f q.right qs

def last (a : ℚ) : List Bin → ℚ
  | [] => a
  | q::qs => last q.right qs

theorem weight_sound {c boxes cb cd bits degree l u P}
    (hb : boxes? c cb cd = some boxes)
    (h : weight? c boxes l u bits degree = some P) :
    0<l ∧ l≤u ∧ ∀ k ∈ Icc (l:ℝ) u, P.Contains (weightReal c k) := by
  unfold weight? at h
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  have hd : 0 < l ∧ l ≤ u := by by_contra hn; simp [guard, hn] at hg
  obtain ⟨T,hT,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨V,hV,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  refine ⟨hd.1,hd.2,?_⟩
  intro k hk
  have tr := rangeFast_sound (a := coefficients c) (K := ⟨l,u⟩) hk
    (fun _ _ hi => coefficientBoxes_sound hb hi) hT
  have sqt := RationalInterval.contains_mul tr tr
  rw [← sq] at sqt
  have ht := SashimiNumerics.exprRange_sound (tiltExpr c) bits degree _ _ _ hV k hk
  apply RationalInterval.contains_mul ht
  exact ⟨by simpa only [Rat.cast_max, Rat.cast_zero] using max_le (sq_nonneg _) sqt.1, sqt.2⟩

theorem build_sound {c boxes cb cd bits degree a ends bins}
    (hb : boxes? c cb cd = some boxes) (h : build? c boxes bits degree a ends = some bins) :
    Good (weightReal c) a bins := by
  induction ends generalizing a bins with
  | nil => cases h; trivial
  | cons b bs ih =>
    obtain ⟨P,hP,h'⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨rest,hr,hI⟩ := Option.bind_eq_some_iff.mp h'
    cases hI
    obtain ⟨ha,hab,hv⟩ := weight_sound hb hP
    exact ⟨ha,hab,hv,ih hr⟩

theorem build_last {c boxes bits degree a ends bins}
    (h : build? c boxes bits degree a ends = some bins) : last a bins = ends.getLastD a := by
  induction ends generalizing a bins with
  | nil => cases h; rfl
  | cons b bs ih =>
    obtain ⟨_,_,h'⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨rest,hr,hI⟩ := Option.bind_eq_some_iff.mp h'
    cases hI
    exact (List.getLastD_cons ..).symm ▸ ih hr

theorem cellIntegral_sound {f : ℝ → ℝ} {l u : ℚ} {I : RationalInterval} (hm : Measurable f)
    (hlu : l ≤ u) (hb : ∀ k∈Icc (l:ℝ) u, I.Contains (f k)) :
    IntervalIntegrable f volume l u ∧ (I.scaleNonnegative (u-l)).Contains (∫ k in (l:ℝ)..u,f k) := by
  let range := fun a b : ℚ => if a=l ∧ b=u then some I else none
  have hb' : ∀ a b J, range a b = some J → ∀ k∈Icc (a:ℝ) b, J.Contains (f k) := by
    intro a b J hj
    dsimp [range] at hj
    split at hj
    · rename_i hv; obtain ⟨rfl,rfl⟩ := hv; cases hj; exact hb
    · contradiction
  have hv : SashimiNumerics.rectangles? range l u 0 = some (I.scaleNonnegative (u-l)) := by
    simp [range,SashimiNumerics.rectangles?]
  have hi := SashimiNumerics.rectangles_integrable hm hb' hlu hv
  exact ⟨hi,SashimiNumerics.rectangles_sound hb' hlu hi hv⟩

theorem totalBy_sound {w f : ℝ → ℝ} {a bins} {piece : ℚ → Bin → RationalInterval} (hgood : Good w a bins)
    (hp : ∀ (l : ℚ) (q : Bin), 0 < l → l ≤ q.right → (∀ k ∈ Icc (l:ℝ) q.right, q.power.Contains (w k)) →
      IntervalIntegrable f volume l q.right ∧ (piece l q).Contains (∫ k in (l:ℝ)..q.right, f k)) :
    IntervalIntegrable f volume a (last a bins) ∧
      (totalBy piece a bins).Contains (∫ k in (a:ℝ)..(last a bins:ℝ), f k) := by
  induction bins generalizing a with
  | nil => simp [last,totalBy,RationalInterval.point,RationalInterval.Contains]
  | cons q qs ih =>
    obtain ⟨hl,hu,hp',hrest⟩ := hgood
    obtain ⟨hi,hI⟩ := hp _ _ hl hu hp'
    obtain ⟨hj,hJ⟩ := ih hrest
    exact ⟨hi.trans hj, (intervalIntegral.integral_add_adjacent_intervals hi hj) ▸ RationalInterval.contains_add hI hJ⟩

theorem round_sound (bits : Nat) {I : RationalInterval} {x : ℝ} (h : I.Contains x) :
    (rounded bits I).Contains x := by
  let s : ℚ := 2^bits
  have hs : 0 < s := by dsimp [s]; positivity
  have hl : (⌊I.lo*s⌋ : ℚ)/s ≤ I.lo := (div_le_iff₀ hs).2 (Int.floor_le _)
  have hu : I.hi ≤ (⌈I.hi*s⌉ : ℚ)/s := (le_div_iff₀ hs).2 (Int.le_ceil _)
  exact ⟨(by exact_mod_cast hl : (((rounded bits I).lo : ℚ) : ℝ) ≤ I.lo).trans h.1,
    h.2.trans (by exact_mod_cast hu : (I.hi : ℝ) ≤ ((rounded bits I).hi : ℝ))⟩

theorem total_sound {c a bins} {R : RationalInterval} {r bits degree} (hgood : Good (weightReal c) a bins)
    (hR : R.Contains r) (hR0 : 0 ≤ R.lo) :
    (total R bits degree a bins).Contains (SashimiEisensteinIntegral.raw c a (last a bins) r) := by
  let f := fun k => weightReal c k * (topHat (k*r))^2
  have hf : Measurable f := SashimiAnalyticWindow.measurableIntegrand (weighted c)
    (fun i => ((coefficients c).toList[i]?).getD 0) r
  apply (totalBy_sound (f := f) (piece := panel R bits degree) hgood ?_).2
  intro l q hl hu hp
  let W := SashimiTableWindow.windowSq (l*R.lo) (q.right*R.hi) bits degree
  have hi := cellIntegral_sound (I := q.power.mul W) hf hu (by
    intro k hk
    apply RationalInterval.contains_mul (hp k hk)
    have hr0 : 0 ≤ (R.lo : ℝ) := by exact_mod_cast hR0
    have hl0 : 0 ≤ (l : ℝ) := by exact_mod_cast hl.le
    apply SashimiTableWindow.windowSq_sound (mul_nonneg hl.le hR0)
    · push_cast; exact mul_le_mul hk.1 hR.1 hr0 (hl0.trans hk.1)
    · push_cast; exact mul_le_mul hk.2 hR.2 (hr0.trans hR.1) ((hl0.trans hk.1).trans hk.2))
  exact ⟨hi.1,round_sound _ hi.2⟩

/-- Production analytic band approaches positive wavenumber from above, so
the real-powers convention at nonpositive k is irrelevant. -/
theorem raw_eq_scaled {a b : ℚ} (ha : 0 < a) (hab : a ≤ b) (r : ℝ) :
    SashimiEisensteinIntegral.raw matched a b r = (2*π^2)*rawAtRadius r a b := by
  unfold rawAtRadius SashimiEisensteinIntegral.raw
  rw [← intervalIntegral.integral_const_mul]
  apply intervalIntegral.integral_congr
  intro k hk
  rw [uIcc_of_le (by exact_mod_cast hab : (a:ℝ) ≤ b)] at hk
  have hp := (by exact_mod_cast ha : (0:ℝ) < a).trans_le hk.1
  dsimp only
  rw [tilt_eq matched hp, Real.rpow_add hp]
  have hc : coefficients matched = prepareReal 0.674 0.315 0.049 2.7255 := by
    norm_num [coefficients, matched]
  rw [hc]
  simp only [matched, realPower, Real.rpow_two]
  norm_num
  field_simp

theorem radius_eq (m : ℚ) :
    SashimiEisensteinIntegral.radius matched m = massRadius massDensity m := by
  rw [radius_identity]
  unfold SashimiEisensteinIntegral.radius matched
  push_cast
  ring_nf

theorem normalized_eq {a b : ℚ} (ha : 0 < a) (hab : a ≤ b) (m : ℚ) :
    normalize8 matched.sigma8 (SashimiEisensteinIntegral.raw matched a b 8)
      (SashimiEisensteinIntegral.raw matched a b (SashimiEisensteinIntegral.radius matched m)) =
        normalizedS m a b := by
  rw [raw_eq_scaled ha hab, raw_eq_scaled ha hab, radius_eq]
  unfold normalize8 normalizedS
  norm_num [matched]; field_simp

/-- An output from the cached evaluator refines the same S as the real mass
specification used for the high-resolution comparison. -/
theorem normalized_sound {m tol a bs bits degree A I}
    (ha : 0 < a) (hab : a ≤ last a bs)
    (hg : Good (weightReal matched) a bs)
    (hA : A = total (.point 8) bits degree a bs)
    (h : normalize? matched m tol bs a A bits degree = some I) :
    I.Contains (normalizedS m a (last a bs)) ∧ I.hi - I.lo ≤ tol := by
  unfold normalize? at h
  obtain ⟨R,hR,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨_,hR0,hap⟩ := Option.bind_eq_some_iff.mp h
  have h0 : 0 ≤ R.lo := by by_contra hn; simp [guard, hn] at hR0
  have hn := total_sound (bits := bits) (degree := degree) hg (radiusBox_sound hR) h0
  have han : A.Contains (SashimiEisensteinIntegral.raw matched a (last a bs) 8) := by
    rw [hA]
    apply total_sound hg (R := RationalInterval.point 8)
    · simpa using RationalInterval.contains_point (8:ℚ)
    · norm_num [RationalInterval.point]
  rw [← normalized_eq ha hab m]
  exact SashimiNormalization.apply_sound hn han hap

#print axioms normalized_sound
end SashimiEhSweep
