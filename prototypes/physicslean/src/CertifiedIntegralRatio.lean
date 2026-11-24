import SashimiEhFullGrid
import SashimiEhSweep

/-! Enclose a ratio of integrals by testing the sign of N − qD. Each cell
uses its common weight only once. A proposed endpoint may come from any
optimizer; the final sign tests and tail allowances use rational arithmetic. -/
namespace CertifiedIntegralRatio
open FloatLib.Numerics SashimiEhSweep SashimiEisensteinIntegral SashimiEhFullGrid
open SashimiSpectrum MeasureTheory Set

structure Cell where
  weight : RationalInterval
  signal : RationalInterval
  anchor : RationalInterval

/-- The weight includes cell width. The signal and anchor range over the cell. -/
def contribution (q : ℚ) (c : Cell) : RationalInterval :=
  c.weight.mul (c.signal.sub (c.anchor.scaleNonnegative q))

def window (r : RationalInterval) (l u : ℚ) (bits degree : Nat) : RationalInterval :=
  rounded bits (SashimiTableWindow.windowSq (l*r.lo) (u*r.hi) bits degree)

def cell (r : RationalInterval) (bits degree : Nat) (l : ℚ) (b : Bin) : Cell :=
  ⟨rounded bits (b.power.scaleNonnegative (b.right-l)),
   window r l b.right bits degree, window (.point 8) l b.right bits degree⟩

def residual (r : RationalInterval) (bits degree : Nat) (q : ℚ) :=
  totalBy (fun l b => contribution q (cell r bits degree l b))

/-- Bounds on a signed cell integral from a bound on its integrand times
the cell width. Works with arbitrary measurable integrands and positive widths. -/
theorem scaledCell {f : ℝ → ℝ} {l u : ℚ} {J : RationalInterval}
    (hlu : l < u) (hm : Measurable f)
    (hJ : ∀ k ∈ Icc (l:ℝ) u, J.Contains (((u-l:ℚ):ℝ)*f k)) :
    IntervalIntegrable f volume l u ∧ J.Contains (∫ k in (l:ℝ)..u, f k) := by
  have hpos : 0 < u-l := sub_pos.mpr hlu
  have hreal : (((u-l:ℚ):ℝ)) ≠ 0 := by exact_mod_cast hpos.ne'
  have h := cellIntegral_sound (I := J.scaleNonnegative (u-l)⁻¹) hm hlu.le (by
    intro k hk
    have hh := RationalInterval.contains_scaleNonnegative (hJ k hk) (inv_nonneg.mpr hpos.le)
    simpa only [Rat.cast_inv,← mul_assoc,inv_mul_cancel₀ hreal,one_mul] using hh)
  obtain ⟨hi,hv⟩ := h
  have he : (J.scaleNonnegative (u-l)⁻¹).scaleNonnegative (u-l) = J := by
    cases J
    simp [RationalInterval.scaleNonnegative, hpos.ne']
  exact ⟨hi,he ▸ hv⟩

/-- All dependence on the particular window is supplied by the two pointwise
range proofs. The common weight may vary arbitrarily within its enclosure. -/
theorem weightedCell {w n d : ℝ → ℝ} {l u q : ℚ} {P X A : RationalInterval} {bits : Nat}
    (hlu : l < u) (hq : 0 ≤ q) (hm : Measurable (fun k => w k * (n k - q*d k)))
    (hw : ∀ k ∈ Icc (l:ℝ) u, P.Contains (w k))
    (hn : ∀ k ∈ Icc (l:ℝ) u, X.Contains (n k))
    (hd : ∀ k ∈ Icc (l:ℝ) u, A.Contains (d k)) :
    IntervalIntegrable (fun k => w k * (n k - q*d k)) volume l u ∧
    (contribution q ⟨rounded bits (P.scaleNonnegative (u-l)), X, A⟩).Contains
      (∫ k in (l:ℝ)..u,w k * (n k - q*d k)) := by
  apply scaledCell hlu hm
  intro k hk
  have W := round_sound bits (RationalInterval.contains_scaleNonnegative (hw k hk)
    (sub_nonneg.mpr hlu.le))
  simpa only [contribution,mul_assoc] using RationalInterval.contains_mul W
    (RationalInterval.contains_sub (hn k hk) (RationalInterval.contains_scaleNonnegative (hd k hk) hq))

theorem window_sound {r : RationalInterval} {l u : ℚ} {x k : ℝ} {bits degree : Nat}
    (hl : 0 < l) (hr : 0 ≤ r.lo) (hx : r.Contains x) (hk : k ∈ Icc (l:ℝ) u) :
    (window r l u bits degree).Contains ((topHat (k*x))^2) := by
  have hr' : 0 ≤ (r.lo : ℝ) := by exact_mod_cast hr
  have hl' : 0 ≤ (l : ℝ) := by exact_mod_cast hl.le
  apply round_sound
  apply SashimiTableWindow.windowSq_sound (mul_nonneg hl.le hr)
  · push_cast; exact mul_le_mul hk.1 hx.1 hr' (hl'.trans hk.1)
  · push_cast; exact mul_le_mul hk.2 hx.2 (hr'.trans hx.1) ((hl'.trans hk.1).trans hk.2)

/-- For each bin, the same spectral function weights numerator and denominator. -/
theorem residual_sound {r : RationalInterval} {x : ℝ} {bins : List Bin} {q : ℚ} {bits degree : Nat}
    (hg : Good (weightReal matched) lowEdge bins) (he : last lowEdge bins = highEdge)
    (hx : r.Contains x) (hr : 0 ≤ r.lo) (hq : 0 ≤ q) :
    (residual r bits degree q lowEdge bins).Contains
      (raw matched lowEdge highEdge x - q*raw matched lowEdge highEdge 8) := by
  let w := weightReal matched
  let n := fun k : ℝ => (topHat (k*x))^2
  let d := fun k : ℝ => (topHat (k*8))^2
  let f := fun k : ℝ => w k * (n k - q*d k)
  have hf : f = fun k => w k*n k - q*(w k*d k) := by funext k; dsimp [f]; ring
  have meas (z : ℝ) : Measurable (fun k => w k * (topHat (k*z))^2) :=
    SashimiAnalyticWindow.measurableIntegrand (weighted matched)
      (fun i => ((coefficients matched).toList[i]?).getD 0) z
  have hmf : Measurable f := by
    rw [hf]; exact (meas x).sub (measurable_const.mul (meas 8))
  have h := (totalBy_sound (f := f) (piece := fun l b => contribution q (cell r bits degree l b)) hg (by
    intro l b hl hu hw
    rcases lt_or_eq_of_le hu with hlt | eq
    · exact weightedCell hlt hq hmf hw (fun _ hk => window_sound hl hr hx hk)
        (fun _ hk => window_sound hl (r := .point 8) (by norm_num [RationalInterval.point])
          (by norm_num [RationalInterval.point,RationalInterval.Contains]) hk)
    · have heq := eq.symm
      simp [heq,cell,contribution,rounded,CertifiedDyadicInterval.rounded,RationalInterval.scaleNonnegative,
        RationalInterval.mul,RationalInterval.Contains])).2
  rw [he,hf, intervalIntegral.integral_sub] at h
  · rwa [intervalIntegral.integral_const_mul] at h
  · exact matched_integrable x ((by exact_mod_cast hr : (0:ℝ) ≤ r.lo).trans hx.1)
  · exact (matched_integrable 8 (by norm_num)).const_mul _

/-- Assemble replayed chunks without asking the kernel to re-evaluate the full partition. -/
theorem residual_append (r : RationalInterval) (bits degree : Nat) (q a : ℚ) (xs ys : List Bin) :
    residual r bits degree q a (xs ++ ys) =
      (residual r bits degree q a xs).add (residual r bits degree q (last a xs) ys) := by
  induction xs generalizing a with
  | nil => simp [residual,totalBy,last,RationalInterval.point,RationalInterval.add_def]
  | cons x xs ih =>
    change (contribution q (cell r bits degree a x)).add (residual r bits degree q x.right (xs ++ ys)) = _
    rw [ih]
    change RationalInterval.add _ (RationalInterval.add _ _) = RationalInterval.add (RationalInterval.add _ _) _
    simp [RationalInterval.add_def,add_assoc,last,residual]

/-- A signed residual on a finite band is padded by separately bounded positive
tails. A positive denominator converts both signs into normalized endpoints. -/
theorem ratio_bounds {n d nf df : ℝ} {l u tn td : ℚ} {L U : RationalInterval}
    (hd : 0 < d) (hl : 0 ≤ l) (hu : 0 ≤ u)
    (hn : nf ≤ n ∧ n ≤ nf + tn) (ha : df ≤ d ∧ d ≤ df + td)
    (hL : L.Contains (nf-l*df)) (hU : U.Contains (nf-u*df))
    (lower : l*td ≤ L.lo) (upper : U.hi+tn ≤ 0) : (⟨l,u⟩:RationalInterval).Contains (n/d) := by
  have h1 : (l:ℝ)*td ≤ L.lo := by exact_mod_cast lower
  have h2 : (U.hi:ℝ)+tn ≤ 0 := by exact_mod_cast upper
  have hl0 : (0:ℝ) ≤ l := by exact_mod_cast hl
  have hu0 : (0:ℝ) ≤ u := by exact_mod_cast hu
  constructor
  · apply (le_div_iff₀ hd).2
    linarith [hL.1,mul_le_mul_of_nonneg_left ha.2 hl0,hn.1]
  · apply (div_le_iff₀ hd).2
    linarith [hU.2,mul_le_mul_of_nonneg_left ha.1 hu0,hn.2]

/-- End-to-end normalized, improper EH98 variance from shared-weight residuals. -/
theorem full_sound {x : ℝ} {r l u : ℚ} {L U A : RationalInterval}
    (hr : (r:ℝ) ≤ x) (hp : 0 < r) (hk : 1 ≤ highEdge*r)
    (hA : A.Contains (raw matched lowEdge highEdge 8)) (ha : 0 < A.lo)
    (hL : L.Contains (raw matched lowEdge highEdge x-l*raw matched lowEdge highEdge 8))
    (hU : U.Contains (raw matched lowEdge highEdge x-u*raw matched lowEdge highEdge 8))
    (h : 0 ≤ l ∧ 0 ≤ u ∧ l*allowance 8 ≤ L.lo ∧ U.hi+allowance r ≤ 0) :
    ((⟨l,u⟩:RationalInterval).scaleNonnegative (matched.sigma8^2)).Contains
      (normalize8 matched.sigma8 (SashimiEhFull.fullRaw matched 8) (SashimiEhFull.fullRaw matched x)) := by
  have A8 := fullRaw_tail (r := (8:ℝ)) (rmin := 8) (by norm_num) (by norm_num) (by norm_num [highEdge])
  have hx := fullRaw_tail hr hp hk
  have hd : 0 < SashimiEhFull.fullRaw matched 8 :=
    ((by exact_mod_cast ha : (0:ℝ) < A.lo).trans_le hA.1).trans_le A8.1
  have hv := ratio_bounds hd h.1 h.2.1 hx A8 hL hU h.2.2.1 h.2.2.2
  simpa only [normalize8,Rat.cast_pow,mul_div_assoc] using
    (RationalInterval.contains_scaleNonnegative hv (sq_nonneg matched.sigma8))

end CertifiedIntegralRatio
