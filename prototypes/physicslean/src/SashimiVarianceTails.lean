import SashimiSpectrumWindow
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals

/-! Missing low- and high-wavenumber contributions to top-hat variance.
These bounds are for the unnormalised integral before the common factor 1/(2π²).
They can be applied separately to numerator and σ₈ anchor, with an actual transfer envelope. -/
namespace SashimiVarianceTails
open Real MeasureTheory Set SashimiSpectrum

noncomputable def weight (T : ℝ → ℝ) (ns R k : ℝ) : ℝ :=
  k^2 * k^ns * (T k)^2 * topHat (k*R)^2

private theorem nonneg {T : ℝ → ℝ} {ns R k : ℝ} (hk : 0 ≤ k) : 0 ≤ weight T ns R k := by
  unfold weight
  positivity

private theorem dominated {s : Set ℝ} {f g : ℝ → ℝ} (hs : MeasurableSet s)
    (hf : Measurable f) (hg : IntegrableOn g s)
    (hb : ∀ k ∈ s, 0 ≤ f k ∧ f k ≤ g k) : IntegrableOn f s := by
  apply hg.mono' hf.aestronglyMeasurable
  filter_upwards [ae_restrict_mem hs] with k hk
  simpa only [Real.norm_eq_abs, abs_of_nonneg (hb k hk).1] using (hb k hk).2

/-- On 0 < k ≤ L ≤ 1, a uniform transfer envelope gives an additive variance budget. -/
theorem low {T : ℝ → ℝ} {ns R L C : ℝ} (hm : Measurable (weight T ns R))
    (hn : 0 ≤ ns) (hR : 0 ≤ R) (hL : 0 < L) (hL1 : L ≤ 1)
    (henv : ∀ k ∈ Ioc 0 L, |T k| ≤ C) :
    IntegrableOn (weight T ns R) (Ioc 0 L) ∧
      (∫ k in (0:ℝ)..L, weight T ns R k) ≤ C^2 * L^3 := by
  have bound (k : ℝ) (hk : k ∈ Ioc 0 L) : weight T ns R k ≤ L^2*C^2 := by
    have ht := henv k hk
    have ht2 : (T k)^2 ≤ C^2 := sq_le_sq.mpr (ht.trans (le_abs_self C))
    have hw := topHat_bounds (k*R) (mul_nonneg hk.1.le hR)
    have hw2 : topHat (k*R)^2 ≤ 1^2 := sq_le_sq.mpr (by simpa using abs_le.mpr hw)
    unfold weight
    calc
      _ ≤ L^2 * 1 * C^2 * 1^2 := by
        have hk0 := hk.1.le
        have hkL := hk.2
        have hns := Real.rpow_le_one hk0 (hkL.trans hL1) hn
        have hp := Real.rpow_nonneg hk0 ns
        gcongr
      _ = _ := by ring
  have hconst : IntegrableOn (fun _ : ℝ => L^2*C^2) (Ioc 0 L) :=
    integrableOn_const (by simp [Real.volume_Ioc])
  have hi : IntegrableOn (weight T ns R) (Ioc 0 L) :=
    dominated measurableSet_Ioc hm hconst (fun k hk => ⟨nonneg hk.1.le, bound k hk⟩)
  refine ⟨hi, ?_⟩
  calc
    _ ≤ ∫ _ in (0:ℝ)..L, L^2*C^2 := by
      simpa only [intervalIntegral.integral_of_le hL.le] using
        setIntegral_mono_on hi hconst measurableSet_Ioc bound
    _ = _ := by rw [intervalIntegral.integral_const]; simp; ring

/-- For k>K, the top-hat envelope contributes four inverse powers after squaring.
The input C/k decay can be conservative. K and KR must both be at least one. -/
theorem high {T : ℝ → ℝ} {ns R K C : ℝ} (hm : Measurable (weight T ns R))
    (hn : ns ≤ 1) (hR : 0 < R) (hK : 1 ≤ K) (hKR : 1 ≤ K*R)
    (henv : ∀ k ∈ Ioi K, |T k| ≤ C/k) :
    IntegrableOn (weight T ns R) (Ioi K) ∧
      (∫ k in Ioi K, weight T ns R k) ≤ 18*C^2/(R^4*K^2) := by
  have hKp : 0 < K := zero_lt_one.trans_le hK
  let A := 36*C^2/R^4
  have bound (k : ℝ) (hk : k ∈ Ioi K) : weight T ns R k ≤ A * k^(-3:ℝ) := by
    have hkp := hKp.trans hk
    have hx : 0 < k*R := mul_pos hkp hR
    have hx1 : 1 ≤ k*R := hKR.trans (mul_le_mul_of_nonneg_right hk.le hR.le)
    have hw : |topHat (k*R)| ≤ 6/(k*R)^2 := by
      apply (topHat_decay _ hx).trans
      apply (div_le_div_iff₀ (by positivity) (by positivity)).2
      nlinarith [mul_nonneg (sub_nonneg.mpr hx1) (sq_nonneg (k*R))]
    have ht2 : (T k)^2 ≤ (C/k)^2 := sq_le_sq.mpr ((henv k hk).trans (le_abs_self _))
    have hw2 : topHat (k*R)^2 ≤ (6/(k*R)^2)^2 := sq_le_sq.mpr (hw.trans (le_abs_self _))
    have hn' : k^ns ≤ k := by simpa using Real.rpow_le_rpow_of_exponent_le (hK.trans hk.le) hn
    unfold weight
    calc
      _ ≤ k^2 * k * (C/k)^2 * (6/(k*R)^2)^2 := by gcongr
      _ = A * k^(-3:ℝ) := by
        rw [Real.rpow_neg hkp.le]; norm_num [A]; field_simp; ring
  have major : IntegrableOn (fun k : ℝ => A*k^(-3:ℝ)) (Ioi K) :=
    (integrableOn_Ioi_rpow_of_lt (by norm_num : (-3:ℝ) < -1) hKp).const_mul A
  have hi := dominated measurableSet_Ioi hm major (fun k hk => ⟨nonneg (hKp.trans hk).le, bound k hk⟩)
  refine ⟨hi, (setIntegral_mono_on hi major measurableSet_Ioi bound).trans_eq ?_⟩
  rw [integral_const_mul, integral_Ioi_rpow_of_lt (by norm_num) hKp]
  norm_num [A, hKp.ne']; ring

/-- Combine a finite-band certificate with both tail contributions. All terms on the
right are nonnegative; numerator and anchor may be covered at different radii. -/
theorem extend {T : ℝ → ℝ} {ns R L K C₀ C₁ : ℝ}
    (hm : Measurable (weight T ns R)) (hn : 0 ≤ ns ∧ ns ≤ 1) (hR : 0 < R)
    (hL : 0 < L) (hL1 : L ≤ 1) (hK : 1 ≤ K) (hKR : 1 ≤ K*R)
    (lo : ∀ k ∈ Ioc 0 L, |T k| ≤ C₀) (hi : ∀ k ∈ Ioi K, |T k| ≤ C₁/k)
    (hband : IntervalIntegrable (weight T ns R) volume L K) :
    let band := ∫ k in L..K, weight T ns R k
    let full := ∫ k in Ioi (0:ℝ), weight T ns R k
    band ≤ full ∧ full ≤ band + C₀^2*L^3 + 18*C₁^2/(R^4*K^2) := by
  dsimp only
  obtain ⟨hlo, hlow⟩ := low hm hn.1 hR.le hL hL1 lo
  obtain ⟨hhi, hhigh⟩ := high hm hn.2 hR hK hKR hi
  have hLK := hL1.trans hK
  have htail : IntegrableOn (weight T ns R) (Ioi L) := by
    rw [← Ioc_union_Ioi_eq_Ioi hLK]
    exact hband.1.union hhi
  have hall : IntegrableOn (weight T ns R) (Ioi (0:ℝ)) := by
    rw [← Ioc_union_Ioi_eq_Ioi hL.le]
    exact hlo.union htail
  have eq₁ := intervalIntegral.integral_Ioi_sub_Ioi hall hL.le
  have eq₂ := intervalIntegral.integral_Ioi_sub_Ioi htail hLK
  have pos₁ : 0 ≤ ∫ k in (0:ℝ)..L, weight T ns R k := by
    rw [intervalIntegral.integral_of_le hL.le]
    exact setIntegral_nonneg measurableSet_Ioc (fun k hk => nonneg hk.1.le)
  have pos₂ : 0 ≤ ∫ k in Ioi K, weight T ns R k :=
    setIntegral_nonneg measurableSet_Ioi (fun k hk => nonneg (zero_le_one.trans (hK.trans hk.le)))
  constructor <;> linarith

end SashimiVarianceTails
