import SashimiEhFullGrid
import SashimiEhDerivativeDecay

/-! Absolute tail allowances for the signed logarithmic mass derivative of the
matched EH98 top-hat variance, before its common σ₈ normalization. -/
namespace SashimiEhSlopeTails
open MeasureTheory Real Set SashimiEisensteinHu SashimiEisensteinIntegral
open SashimiEhFull SashimiEhFullGrid SashimiEhDerivativeBounds SashimiEhSweep

noncomputable def slope (r k : ℝ) := weightReal matched k * shape (k*r)
noncomputable def fullSlope (r : ℝ) := ∫ k in Ioi (0:ℝ), slope r k

theorem measurableSlope (r : ℝ) : Measurable (slope r) := by
  have hw : Measurable SashimiSpectrum.topHat := by
    unfold SashimiSpectrum.topHat
    exact Measurable.ite (measurableSet_singleton 0) measurable_const (by fun_prop)
  have h := SashimiAnalyticWindow.measurableIntegrand (weighted matched)
    (fun i => ((coefficients matched).toList[i]?).getD 0) 0
  change Measurable (fun k => weightReal matched k * SashimiSpectrum.topHat (k*0)^2) at h
  have hp : Measurable (weightReal matched) := by simpa [SashimiSpectrum.topHat] using h
  unfold slope shape
  fun_prop

private theorem matched_coeff : coefficients matched = fiducial := by
  norm_num [coefficients, matched, fiducial]

private theorem smallT {k : ℝ} (hk : 0 < k) : |transferReal (coefficients matched) k| ≤ 3 := by
  rw [matched_coeff]
  exact low_envelope fiducial_bounds hk

private theorem largeT {k : ℝ} (hk : (highEdge:ℝ) ≤ k) :
    |transferReal (coefficients matched) k| ≤ (highConstant:ℝ)/k := by
  rw [matched_coeff]
  have h : (10:ℝ)^10 ≤ k := by simpa only [highEdge, Rat.cast_pow, Rat.cast_ofNat] using hk
  simpa only [highConstant, Rat.cast_div, Rat.cast_one, Rat.cast_pow, Rat.cast_ofNat] using
    high_envelopeBig fiducial_bounds h

private theorem power_eq {k : ℝ} (hk : 0<k) :
    weightReal matched k = k^2*k^(matched.ns:ℝ)*(transferReal (coefficients matched) k)^2 := by
  rw [weightReal, tilt_eq matched hk, Real.rpow_add hk]
  norm_num

private theorem abs_slope {r k : ℝ} (hk : 0 < k) :
    |slope r k| = k^2*k^(matched.ns:ℝ)*(transferReal (coefficients matched) k)^2*|shape (k*r)| := by
  rw [slope, abs_mul, power_eq hk, abs_of_nonneg (by positivity)]

theorem midBound {r k K : ℝ} (hr : 0<r) (hk : 0<k) (hu : k≤K) (hK : 1≤K) :
    |slope r k| ≤ 36*K^3 := by
  have hn : (0:ℝ) ≤ matched.ns := by norm_num [matched]
  have hn1 : (matched.ns:ℝ) ≤ 1 := by norm_num [matched]
  have hp : k^(matched.ns:ℝ) ≤ K :=
    (Real.rpow_le_rpow hk.le hu hn).trans (by simpa using Real.rpow_le_rpow_of_exponent_le hK hn1)
  have ht : (transferReal (coefficients matched) k)^2 ≤ (3:ℝ)^2 :=
    sq_le_sq.mpr ((smallT hk).trans (le_abs_self _))
  rw [abs_slope hk]
  calc
    _ ≤ K^2*K*3^2*4 := by gcongr; exact abs_shape_le_four (mul_pos hk hr)
    _ = _ := by ring

theorem highBound {r rmin K k : ℝ} (hr : 0<rmin) (hle : rmin≤r)
    (hK : (highEdge:ℝ)≤K) (hKR : 1≤K*rmin) (hk : K<k) :
    |slope r k| ≤ (84*(highConstant:ℝ)^2/rmin^3)*k^(-2:ℝ) := by
  have hk1 : 1≤k := (by norm_num [highEdge] : (1:ℝ)≤highEdge).trans (hK.trans hk.le)
  have hkp := zero_lt_one.trans_le hk1
  have hx : 1≤k*r := hKR.trans ((mul_le_mul_of_nonneg_right hk.le hr.le).trans
    (mul_le_mul_of_nonneg_left hle hkp.le))
  have hp : k^(matched.ns:ℝ) ≤ k := by
    simpa using Real.rpow_le_rpow_of_exponent_le hk1 (by norm_num [matched] : (matched.ns:ℝ)≤1)
  have ht : (transferReal (coefficients matched) k)^2 ≤ ((highConstant:ℝ)/k)^2 :=
    sq_le_sq.mpr ((largeT (hK.trans hk.le)).trans (le_abs_self _))
  have hw : |shape (k*r)| ≤ 84/(k*rmin)^3 := (abs_shape_decay hx).trans
    (div_le_div_of_nonneg_left (by norm_num) (by positivity) (by gcongr))
  rw [abs_slope hkp]
  calc
    _ ≤ k^2*k*((highConstant:ℝ)/k)^2*(84/(k*rmin)^3) := by gcongr
    _ = _ := by rw [Real.rpow_neg hkp.le]; norm_num; field_simp

private theorem boundedIntegral {f g : ℝ → ℝ} {s : Set ℝ} (hs : MeasurableSet s)
    (hf : Measurable f) (hg : IntegrableOn g s) (hb : ∀ x∈s, |f x| ≤ g x) :
    IntegrableOn f s ∧ |∫ x in s, f x| ≤ ∫ x in s, g x := by
  have hi : IntegrableOn f s := by
    apply hg.mono' hf.aestronglyMeasurable
    filter_upwards [ae_restrict_mem hs] with x hx
    exact hb x hx
  rw [← Real.norm_eq_abs]
  refine ⟨hi, (norm_integral_le_integral_norm (μ := volume.restrict s) f).trans ?_⟩
  exact setIntegral_mono_on hi.norm hg hs hb

theorem fullRaw_integrable {r rmin : ℝ} (hr : 0<rmin) (hle : rmin≤r)
    (hKR : 1≤(highEdge:ℝ)*rmin) : IntegrableOn (spectrumWeight matched r) (Ioi 0) := by
  let K : ℝ := highEdge
  have hK : 1 ≤ K := by norm_num [K, highEdge]
  have hkp := zero_lt_one.trans_le hK
  have hhi := SashimiVarianceTails.high (measurableWeight matched r)
    (by norm_num [matched]) (hr.trans_le hle) hK (hKR.trans (mul_le_mul_of_nonneg_left hle hkp.le))
    (fun k hk => largeT hk.le)
  rw [← Ioc_union_Ioi_eq_Ioi hkp.le]
  refine IntegrableOn.union ?_ hhi.1
  have hc : IntegrableOn (fun _ : ℝ => K^3*9) (Ioc 0 K) :=
    integrableOn_const (by simp [Real.volume_Ioc])
  apply hc.mono' (measurableWeight matched r).aestronglyMeasurable
  filter_upwards [ae_restrict_mem (μ := volume) (measurableSet_Ioc (a := (0:ℝ)) (b := K))] with k hk
  have h0 := hk.1
  have h1 := hk.2
  have hp : k^(matched.ns:ℝ) ≤ K :=
    (Real.rpow_le_rpow h0.le h1 (by norm_num [matched])).trans
      (by simpa using Real.rpow_le_rpow_of_exponent_le hK (by norm_num [matched] : (matched.ns:ℝ)≤1))
  have ht : (transferReal (coefficients matched) k)^2 ≤ (3:ℝ)^2 :=
    sq_le_sq.mpr ((smallT h0).trans (le_abs_self _))
  have hw : SashimiSpectrum.topHat (k*r)^2 ≤ (1:ℝ)^2 :=
    sq_le_sq.mpr (by simpa using abs_le.mpr (SashimiSpectrum.topHat_bounds _ (mul_nonneg h0.le (hr.le.trans hle))))
  unfold spectrumWeight SashimiVarianceTails.weight
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  calc
    _ ≤ K^2*K*3^2*1^2 := by gcongr
    _ = _ := by ring

/-- The signed tails can move either endpoint. The derivative's lower-frequency
allowance scales as L³ and its higher-frequency allowance as 1/K. -/
theorem tails {r rmin : ℝ} (hr : 0 < rmin) (hle : rmin ≤ r) (hKR : 1 ≤ (highEdge:ℝ)*rmin) :
    |fullSlope r - rawSlope matched r lowEdge highEdge| ≤
      36*(lowEdge:ℝ)^3 + 84*(highConstant:ℝ)^2/(rmin^3*(highEdge:ℝ)) := by
  let L : ℝ := lowEdge
  let K : ℝ := highEdge
  have hL : 0 < L := by norm_num [L,lowEdge]
  have hL1 : L ≤ 1 := by norm_num [L,lowEdge]
  have hK : 1 ≤ K := by norm_num [K,highEdge]
  have hKp := zero_lt_one.trans_le hK
  have hLK := hL1.trans hK
  have hconst (p q c : ℝ) : IntegrableOn (fun _ : ℝ => c) (Ioc p q) :=
    integrableOn_const (by simp [Real.volume_Ioc])
  have ⟨hi0, h0⟩ := boundedIntegral measurableSet_Ioc (measurableSlope r) (hconst 0 L (36*L^2)) (by
    intro k hk
    have hk0 := hk.1
    have hkL := hk.2
    have hn := Real.rpow_le_one hk0.le (hkL.trans hL1) (by norm_num [matched] : 0 ≤ (matched.ns:ℝ))
    have ht : (transferReal (coefficients matched) k)^2 ≤ (3:ℝ)^2 :=
      sq_le_sq.mpr ((smallT hk0).trans (le_abs_self _))
    rw [abs_slope hk0]
    calc
      _ ≤ L^2*1*3^2*4 := by gcongr; exact abs_shape_le_four (mul_pos hk0 (hr.trans_le hle))
      _ = _ := by ring)
  let C := 84*(highConstant:ℝ)^2/rmin^3
  have hmajor : IntegrableOn (fun k : ℝ => C*k^(-2:ℝ)) (Ioi K) :=
    (integrableOn_Ioi_rpow_of_lt (by norm_num : (-2:ℝ)< -1) hKp).const_mul _
  have ⟨hiK,hhigh⟩ := boundedIntegral measurableSet_Ioi (measurableSlope r) hmajor
    (fun k hk => highBound hr hle le_rfl hKR hk)
  have himid := (boundedIntegral measurableSet_Ioc (measurableSlope r) (hconst L K (36*K^3))
    (fun k hk => midBound (hr.trans_le hle) (hL.trans hk.1) hk.2 hK)).1
  have hiL : IntegrableOn (slope r) (Ioi L) := by
    rw [← Ioc_union_Ioi_eq_Ioi hLK]
    exact himid.union hiK
  have hi : IntegrableOn (slope r) (Ioi 0) := by
    rw [← Ioc_union_Ioi_eq_Ioi hL.le]
    exact hi0.union hiL
  have he0 := intervalIntegral.integral_Ioi_sub_Ioi hi hL.le
  have he1 := intervalIntegral.integral_Ioi_sub_Ioi hiL hLK
  rw [intervalIntegral.integral_of_le hL.le] at he0
  have hconstint : (∫ _ in Ioc (0:ℝ) L, 36*L^2) = 36*L^3 := by
    rw [← intervalIntegral.integral_of_le hL.le, intervalIntegral.integral_const]; simp; ring
  rw [hconstint] at h0
  have hhighint : (∫ k in Ioi K, C*k^(-2:ℝ)) = 84*(highConstant:ℝ)^2/(rmin^3*K) := by
    rw [integral_const_mul, integral_Ioi_rpow_of_lt (by norm_num) hKp]; norm_num [C, hKp.le, Real.rpow_neg hKp.le]; ring
  rw [hhighint] at hhigh
  change |(∫ k in Ioi (0:ℝ), slope r k) - (∫ k in L..K, slope r k)| ≤ _
  have he : (∫ k in Ioi (0:ℝ), slope r k) - (∫ k in L..K, slope r k) =
      (∫ k in Ioc (0:ℝ) L, slope r k) + ∫ k in Ioi K, slope r k := by linarith
  rw [he]
  exact (abs_add_le ..).trans (add_le_add h0 hhigh)

#print axioms tails
end SashimiEhSlopeTails
