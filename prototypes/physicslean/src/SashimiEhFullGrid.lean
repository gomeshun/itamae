import SashimiEhFull
import SashimiEisensteinTail

/-! Improper EH98 variance and its independently integrated σ₈ normalization.
The radii are in Mpc/h and physical input masses are in solar masses. -/
namespace SashimiEhFullGrid
open MeasureTheory Set FloatLib.Numerics SashimiEisensteinHu SashimiEisensteinIntegral SashimiEhFull

def lowEdge : ℚ := 1/10^8
def highEdge : ℚ := 10^10
def highConstant : ℚ := 1/10^9

noncomputable def S (m : ℚ) : ℝ :=
  SashimiSpectrum.normalize8 matched.sigma8 (fullRaw matched 8) (fullRaw matched (radius matched m))

private theorem bandIntegrable {r : ℝ} (hr : 0 ≤ r)
    (hT : ∀ k : ℝ, 0 < k → |transferReal (coefficients matched) k| ≤ 3) :
    IntervalIntegrable (spectrumWeight matched r) volume lowEdge highEdge := by
  have hL : (0:ℝ) < lowEdge := by norm_num [lowEdge]
  have hK : (1:ℝ) ≤ highEdge := by norm_num [highEdge]
  have hLK : (lowEdge:ℝ) ≤ highEdge := by norm_num [lowEdge,highEdge]
  rw [intervalIntegrable_iff_integrableOn_Ioc_of_le hLK]
  have hi : IntegrableOn (fun _ : ℝ => (highEdge:ℝ)^3*9) (Ioc (lowEdge:ℝ) highEdge) :=
    integrableOn_const (by simp [Real.volume_Ioc])
  apply hi.mono' (measurableWeight matched r).aestronglyMeasurable
  filter_upwards [ae_restrict_mem (μ := volume) (measurableSet_Ioc (a := (lowEdge:ℝ)) (b := (highEdge:ℝ)))] with k hk
  have hkp : 0 < k := hL.trans hk.1
  have hns : (0:ℝ) ≤ matched.ns := by norm_num [matched]
  have hns1 : (matched.ns:ℝ) ≤ 1 := by norm_num [matched]
  have ht : (transferReal (coefficients matched) k)^2 ≤ (3:ℝ)^2 :=
    sq_le_sq.mpr ((hT k hkp).trans (le_abs_self _))
  have hw : SashimiSpectrum.topHat (k*r)^2 ≤ (1:ℝ)^2 :=
    sq_le_sq.mpr ((abs_le.mpr (SashimiSpectrum.topHat_bounds _ (mul_nonneg hkp.le hr))).trans (le_abs_self _))
  have hp : k^(matched.ns:ℝ) ≤ (highEdge:ℝ) :=
    (Real.rpow_le_rpow hkp.le hk.2 hns).trans
      (by simpa using Real.rpow_le_rpow_of_exponent_le hK hns1)
  unfold spectrumWeight SashimiVarianceTails.weight
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  calc
    _ ≤ (highEdge:ℝ)^2*highEdge*3^2*1^2 := by gcongr; exact hk.2
    _ = _ := by ring

def allowance (Rmin : ℚ) := budget lowEdge highEdge Rmin 3 highConstant

/-- Exact arithmetic checker for padded numerator and denominator, normalized
variance, and an outward-rounded square root. -/
def check (rmin : ℚ) (raw anchor variance sigma : RationalInterval) : Bool :=
  let n := pad raw (allowance rmin)
  let a := pad anchor (allowance 8)
  let v := (n.scaleNonnegative (matched.sigma8^2)).divPositive a
  decide (0<rmin ∧ 1≤highEdge*rmin ∧ 0<a.lo ∧ variance.lo≤v.lo ∧ v.hi≤variance.hi ∧
    0≤sigma.lo ∧ sigma.lo^2 ≤ variance.lo ∧ variance.hi ≤ sigma.hi^2 ∧ 0≤sigma.hi)

private theorem raw_sound {r : ℝ} {rmin : ℚ} {J : RationalInterval}
    (hlo : ∀ k : ℝ, 0 < k → |transferReal (coefficients matched) k| ≤ 3)
    (hhi : ∀ k ∈ Ioi (highEdge:ℝ), |transferReal (coefficients matched) k| ≤ (highConstant:ℝ)/k)
    (hr : (rmin:ℝ)≤r) (hp : 0<rmin) (hk : 1≤highEdge*rmin)
    (hj : J.Contains (raw matched lowEdge highEdge r)) :
    (pad J (allowance rmin)).Contains (fullRaw matched r) := by
  apply pad_sound (by norm_num [matched]) (by norm_num [lowEdge]) (by norm_num [lowEdge])
    (by norm_num [highEdge]) hp hr hk (fun k h => by simpa using hlo k h.1) hhi
    (bandIntegrable ((by exact_mod_cast hp.le : (0:ℝ) ≤ rmin).trans hr) hlo) hj

theorem check_with_envelopes {m rmin : ℚ} {j a v sigma : RationalInterval}
    (hlo : ∀ k : ℝ, 0 < k → |transferReal (coefficients matched) k| ≤ 3)
    (hhi : ∀ k ∈ Ioi (highEdge:ℝ), |transferReal (coefficients matched) k| ≤ (highConstant:ℝ)/k)
    (hr : (rmin:ℝ) ≤ radius matched m)
    (hn : j.Contains (raw matched lowEdge highEdge (radius matched m)))
    (ha : a.Contains (raw matched lowEdge highEdge 8)) (hc : check rmin j a v sigma = true) :
    v.Contains (S m) ∧ sigma.Contains (Real.sqrt (S m)) := by
  let N := pad j (allowance rmin)
  let A := pad a (allowance 8)
  let V := (N.scaleNonnegative (matched.sigma8^2)).divPositive A
  change decide (_ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _) = true at hc
  obtain ⟨hm,hK,hA,hl,hu,_,hs1,hs2,hsp⟩ := of_decide_eq_true hc
  have hN := raw_sound hlo hhi hr hm hK hn
  have hAn := raw_sound (rmin:=8) hlo hhi (by norm_num) (by norm_num) (by norm_num [highEdge]) ha
  have hv : V.Contains (S m) := by
    simpa only [S, SashimiSpectrum.normalize8,Rat.cast_pow,mul_comm] using
      RationalInterval.contains_divPositive
        (RationalInterval.contains_scaleNonnegative hN (sq_nonneg matched.sigma8)) hAn hA
  have hV : v.Contains (S m) := ⟨(show (v.lo:ℝ)≤V.lo by exact_mod_cast hl).trans hv.1,
    hv.2.trans (by exact_mod_cast hu)⟩
  exact ⟨hV, Real.le_sqrt_of_sq_le ((by exact_mod_cast hs1 : (sigma.lo:ℝ)^2≤v.lo).trans hV.1),
    Real.sqrt_le_iff.mpr ⟨by exact_mod_cast hsp, hV.2.trans (by exact_mod_cast hs2)⟩⟩

private theorem matchedCoefficients : coefficients matched = fiducial := by
  norm_num [coefficients, matched, fiducial]

private theorem matchedLow (k : ℝ) (hk : 0 < k) :
    |transferReal (coefficients matched) k| ≤ 3 := by
  rw [matchedCoefficients]
  exact low_envelope fiducial_bounds hk

private theorem matchedHigh (k : ℝ) (hk : k ∈ Ioi (highEdge:ℝ)) :
    |transferReal (coefficients matched) k| ≤ (highConstant:ℝ)/k := by
  rw [matchedCoefficients]
  change (highEdge:ℝ) < k at hk
  have hk' : (10:ℝ)^10 ≤ k := by
    simpa only [highEdge, Rat.cast_pow, Rat.cast_ofNat] using (le_of_lt hk : (highEdge:ℝ) ≤ k)
  simpa only [highConstant, Rat.cast_div, Rat.cast_one, Rat.cast_pow, Rat.cast_ofNat] using
    high_envelopeBig fiducial_bounds hk'

/-- Padded unnormalized integral at any certified radius. Multiple radii can
share a single padded denominator when bounding a difference of variances. -/
theorem fullRaw_sound {r : ℝ} {rmin : ℚ} {j : RationalInterval}
    (hr : (rmin:ℝ) ≤ r) (hp : 0<rmin) (hk : 1 ≤ highEdge*rmin)
    (hj : j.Contains (raw matched lowEdge highEdge r)) :
    (pad j (allowance rmin)).Contains (fullRaw matched r) :=
  raw_sound matchedLow matchedHigh hr hp hk hj

theorem matched_integrable (r : ℝ) (hr : 0≤r) :
    IntervalIntegrable (fun k => Real.exp (((2+matched.ns:ℚ):ℝ)*Real.log k)*
      (transferReal (coefficients matched) k)^2*SashimiSpectrum.topHat (k*r)^2) volume lowEdge highEdge := by
  apply (bandIntegrable hr matchedLow).congr
  intro k hk
  rw [uIoc_of_le (by norm_num [lowEdge, highEdge] : (lowEdge:ℝ) ≤ highEdge)] at hk
  exact weight_eq matched r ((by norm_num [lowEdge] : (0:ℝ)<lowEdge).trans hk.1)

/-- The numerical tail allowance also bounds the *additive* remainder. Keeping
the raw difference intact avoids two independent enclosures of common power. -/
theorem fullRaw_tail {r : ℝ} {rmin : ℚ} (hr : (rmin:ℝ) ≤ r) (hp : 0<rmin) (hk : 1≤highEdge*rmin) :
    raw matched lowEdge highEdge r ≤ fullRaw matched r ∧
      fullRaw matched r ≤ raw matched lowEdge highEdge r + (allowance rmin:ℝ) := by
  let x := raw matched lowEdge highEdge r
  let y := fullRaw matched r
  have H (J : RationalInterval) (hj : J.Contains x) := fullRaw_sound hr hp hk hj
  constructor
  · by_contra! he
    obtain ⟨q, hq, hx⟩ := exists_rat_btwn he
    have hj : (⟨q,(⌈x⌉:ℤ)⟩:RationalInterval).Contains x :=
      ⟨hx.le, by simpa using Int.le_ceil x⟩
    exact hq.not_ge (H _ hj).1
  · by_contra! he
    obtain ⟨q,hx,hq⟩ := exists_rat_btwn (show x < y - (allowance rmin:ℝ) by dsimp [x,y]; linarith)
    have hj : (⟨(⌊x⌋:ℤ),q⟩:RationalInterval).Contains x :=
      ⟨by simpa using Int.floor_le x, hx.le⟩
    have hu := (H _ hj).2
    change y ≤ (((q + allowance rmin : ℚ)) : ℝ) at hu
    push_cast at hu
    linarith

/-- All tail envelopes are discharged for the matched cosmology. Only the
finite-band quadrature, actual radius and rational endpoint certificate remain. -/
theorem check_sound {m rmin : ℚ} {j a v sigma : RationalInterval}
    (hr : (rmin:ℝ) ≤ radius matched m)
    (hn : j.Contains (raw matched lowEdge highEdge (radius matched m)))
    (ha : a.Contains (raw matched lowEdge highEdge 8)) (hc : check rmin j a v sigma = true) :
    v.Contains (S m) ∧ sigma.Contains (Real.sqrt (S m)) :=
  check_with_envelopes matchedLow matchedHigh hr hn ha hc

#print axioms fullRaw_sound
#print axioms check_sound

end SashimiEhFullGrid
