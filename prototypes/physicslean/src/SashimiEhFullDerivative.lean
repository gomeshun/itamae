import SashimiEhFullCalculus

/-! Signed derivative of the improper, independently σ₈-normalized EH variance
as a function of physical mass in solar masses. -/
namespace SashimiEhFullDerivative
open MeasureTheory Real FloatLib.Numerics SashimiEisensteinIntegral SashimiEhSweep
open SashimiEhFull SashimiEhFullGrid SashimiEhDerivativeBounds SashimiEhSlopeTails SashimiEhVariance

/-- Real-mass extension of the mass-grid specification. -/
noncomputable def fullS (m : ℝ) : ℝ :=
  SashimiSpectrum.normalize8 matched.sigma8 (fullRaw matched 8)
    (fullRaw matched (SashimiSpectrum.massRadius massDensity m))

theorem fullS_eq_grid (m : ℚ) : fullS m = SashimiEhFullGrid.S m := by
  rw [fullS, SashimiEhFullGrid.S, radius_eq]

theorem hasDerivAt_fullS {m rmin : ℝ} (hm : 0<m) (hr : 0<rmin)
    (hle : rmin < SashimiSpectrum.massRadius massDensity m) (hKR : 1≤(highEdge:ℝ)*rmin) :
    HasDerivAt fullS (SashimiSpectrum.normalize8 matched.sigma8 (fullRaw matched 8)
      ((1/m)*fullSlope (SashimiSpectrum.massRadius massDensity m))) m := by
  have hd := SashimiEhFullCalculus.hasDerivAt_fullRaw hm hr hle hKR
  apply ((hd.const_mul ((matched.sigma8:ℝ)^2)).div_const (fullRaw matched 8)).congr_deriv
  unfold SashimiSpectrum.normalize8
  ring

def signedBudget (r : ℚ) : ℚ := 36*lowEdge^3 + 84*highConstant^2/(r^3*highEdge)

def signedPad (i : RationalInterval) (q : ℚ) : RationalInterval := ⟨i.lo-q, i.hi+q⟩

/-- Half the positive lower radius bound leaves a neighborhood for uniform
dominated differentiation. Both integral tails and the normalization are
checked before the proposed outward endpoints are accepted. -/
def check (m : ℚ) (R N A I : RationalInterval) : Bool :=
  let r := R.lo/2
  let raw := (signedPad N (signedBudget r)).scaleNonnegative (1/m)
  let anchor := pad A (allowance 8)
  let J := (raw.scaleNonnegative (matched.sigma8^2)).divPositive anchor
  decide (0<m ∧ 0<r ∧ 1≤highEdge*r ∧ 0<anchor.lo ∧ I.lo≤J.lo ∧ J.hi≤I.hi)

theorem sound {m : ℚ} {R N A I : RationalInterval} (hR : R.Contains (radius matched m))
    (hn : N.Contains (rawSlope matched (radius matched m) lowEdge highEdge))
    (ha : A.Contains (raw matched lowEdge highEdge 8)) (hc : check m R N A I = true) :
    I.Contains (deriv fullS (m:ℝ)) := by
  let r := R.lo / 2
  let P := signedPad N (signedBudget r)
  let B := pad A (allowance 8)
  let J := ((P.scaleNonnegative (1/m)).scaleNonnegative (matched.sigma8^2)).divPositive B
  change decide (0<m ∧ 0<r ∧ 1≤highEdge*r ∧ 0<B.lo ∧ I.lo≤J.lo ∧ J.hi≤I.hi) = true at hc
  obtain ⟨hm,hr,hkr,hb,hl,hu⟩ := of_decide_eq_true hc
  have hrr : r < R.lo := by dsimp [r] at hr ⊢; linarith
  have hp : (0:ℝ)<r := by exact_mod_cast hr
  have hle := (by exact_mod_cast hrr : (r:ℝ)<R.lo).trans_le hR.1
  have hk : 1 ≤ (highEdge:ℝ)*(r:ℝ) := by exact_mod_cast hkr
  have ht := abs_le.mp (tails hp hle.le hk)
  have hP : P.Contains (fullSlope (radius matched m)) := by
    change ((N.lo-signedBudget r : ℚ):ℝ) ≤ _ ∧ _ ≤ ((N.hi+signedBudget r : ℚ):ℝ)
    unfold signedBudget
    push_cast
    constructor <;> linarith [hn.1,hn.2,ht.1,ht.2]
  have hA := fullRaw_sound (rmin:=8) (by norm_num) (by norm_num) (by norm_num [highEdge]) ha
  have hv := RationalInterval.contains_scaleNonnegative hP (one_div_nonneg.mpr hm.le)
  have hJ := RationalInterval.contains_divPositive
    (RationalInterval.contains_scaleNonnegative hv (sq_nonneg matched.sigma8)) hA hb
  have he := (hasDerivAt_fullS (m:=(m:ℝ)) (by exact_mod_cast hm) hp (by rwa [← radius_eq]) hk).deriv
  rw [he, ← radius_eq]
  have hnormalized : J.Contains (SashimiSpectrum.normalize8 matched.sigma8 (fullRaw matched 8)
        ((1/(m:ℝ))*fullSlope (radius matched m))) := by
    simpa only [SashimiSpectrum.normalize8, Rat.cast_pow, Rat.cast_div, Rat.cast_one] using hJ
  exact ⟨(by exact_mod_cast hl : (I.lo:ℝ)≤J.lo).trans hnormalized.1,
    hnormalized.2.trans (by exact_mod_cast hu)⟩

theorem replayed {m : ℚ} {R N A I bits degree bins}
    (hg : Good (weightReal matched) lowEdge bins) (he : last lowEdge bins = highEdge)
    (hr : radiusBox? matched m bits degree = some R) (hp : 0 < R.lo)
    (hn : sum R bits degree lowEdge bins = N) (ha : A.Contains (raw matched lowEdge highEdge 8))
    (hc : check m R N A I = true) : I.Contains (deriv fullS (m:ℝ)) := by
  have h := sum_sound (bits:=bits) (degree:=degree) hg (radiusBox_sound hr) hp
  rw [hn,he] at h
  exact sound (radiusBox_sound hr) h ha hc

#print axioms hasDerivAt_fullS
#print axioms replayed
end SashimiEhFullDerivative
