import CertifiedCenteredIntegral
import FloatLib.Numerics.Enclosure.Rational.Converges

/-! Scalar integration rules for enclosed, signed integrands. The weighted rule
uses the exact integral of a nonnegative weight. The affine-band rule integrates
pointwise witnesses that may depend on the point and need no regularity. -/
namespace CertifiedWeightedIntegral
open Set MeasureTheory FloatLib.Numerics

private theorem weighted_value {w g : ℝ → ℝ} {a b : ℝ} {J : RationalInterval}
    (hab : a ≤ b) (hw : IntervalIntegrable w volume a b)
    (hwg : IntervalIntegrable (fun x => w x * g x) volume a b)
    (hw0 : ∀ x ∈ Icc a b, 0 ≤ w x) (hJ : ∀ x ∈ Icc a b, J.Contains (g x)) :
    ∃ v, J.Contains v ∧ (∫ x in a..b, w x * g x) = (∫ x in a..b, w x) * v := by
  let W := ∫ x in a..b, w x
  let I := ∫ x in a..b, w x * g x
  have hW : 0 ≤ W := intervalIntegral.integral_nonneg hab hw0
  have hlo : W * (J.lo:ℝ) ≤ I := by
    rw [show W * (J.lo:ℝ) = ∫ x in a..b, w x * (J.lo:ℝ) by
      rw [intervalIntegral.integral_mul_const]]
    exact intervalIntegral.integral_mono_on hab (hw.mul_const _) hwg
      fun x hx => mul_le_mul_of_nonneg_left (hJ x hx).1 (hw0 x hx)
  have hhi : I ≤ W * (J.hi:ℝ) := by
    rw [show W * (J.hi:ℝ) = ∫ x in a..b, w x * (J.hi:ℝ) by
      rw [intervalIntegral.integral_mul_const]]
    exact intervalIntegral.integral_mono_on hab hwg (hw.mul_const _)
      fun x hx => mul_le_mul_of_nonneg_left (hJ x hx).2 (hw0 x hx)
  rcases hW.eq_or_lt with h | h
  · refine ⟨g a, hJ a ⟨le_rfl,hab⟩, ?_⟩
    change I = W * g a
    rw [← h] at hlo hhi ⊢
    simp only [zero_mul] at hlo hhi ⊢
    exact le_antisymm hhi hlo
  · refine ⟨I/W, ⟨(le_div_iff₀ h).mpr (by simpa [mul_comm] using hlo),
      (div_le_iff₀ h).mpr (by simpa [mul_comm] using hhi)⟩, ?_⟩
    change I = W * (I/W)
    field_simp

/-- Integral-average rule for general signed g and a nonnegative weight. In
particular, a weight with integral zero yields integral zero for every bounded g. -/
theorem weighted_sound {w g : ℝ → ℝ} {a b : ℝ} {A J : RationalInterval}
    (hab : a ≤ b) (hw : IntervalIntegrable w volume a b)
    (hwg : IntervalIntegrable (fun x => w x * g x) volume a b)
    (hw0 : ∀ x ∈ Icc a b, 0 ≤ w x)
    (hA : A.Contains (∫ x in a..b, w x)) (hJ : ∀ x ∈ Icc a b, J.Contains (g x)) :
    (A.mul J).Contains (∫ x in a..b, w x * g x) := by
  obtain ⟨v,hv,he⟩ := weighted_value hab hw hwg hw0 hJ
  rw [he]
  exact RationalInterval.contains_mul hA hv

/-- Exact second moment of Lebesgue measure on an interval with rational endpoints. -/
def squareMoment (a b : ℚ) := (b^3-a^3)/3

theorem squareMoment_eq (a b : ℚ) :
    (∫ k in (a:ℝ)..(b:ℝ), k^2) = (squareMoment a b : ℝ) := by
  rw [integral_pow]
  norm_num [squareMoment]

theorem square_weight_sound {g : ℝ → ℝ} {a b : ℚ} {P : RationalInterval}
    (hab : a ≤ b) (hg : IntervalIntegrable (fun k => k^2 * g k) volume (a:ℝ) b)
    (hP : ∀ k ∈ Icc (a:ℝ) b, P.Contains (g k)) :
    (P.scaleNonnegative (squareMoment a b)).Contains (∫ k in (a:ℝ)..(b:ℝ), k^2 * g k) := by
  have hab' : (a:ℝ) ≤ b := by exact_mod_cast hab
  have hw : IntervalIntegrable (fun k : ℝ => k^2) volume (a:ℝ) b :=
    (continuous_id.pow 2).intervalIntegrable _ _
  obtain ⟨v,hv,he⟩ := weighted_value (w := fun k : ℝ => k^2) hab' hw hg
    (fun k _ => sq_nonneg k) hP
  have hpos : 0 ≤ squareMoment a b := by
    exact_mod_cast (squareMoment_eq a b) ▸ (intervalIntegral.integral_nonneg hab' fun k _ => sq_nonneg k)
  rw [he,squareMoment_eq]
  exact RationalInterval.contains_scaleNonnegative hv hpos

/-- The two witnesses need not be measurable functions of k. Existence at
each point suffices. Measurability of the actual integrand is still required. -/
theorem affine_band {f : ℝ → ℝ} {a b : ℚ} {C D : RationalInterval}
    (hab : a ≤ b) (hf : Measurable f)
    (hband : ∀ k ∈ Icc (a:ℝ) b, ∃ v : ℝ, C.Contains v ∧ ∃ d : ℝ, D.Contains d ∧
      f k = v + d * (k - (((a+b)/2 : ℚ):ℝ))) :
    IntervalIntegrable f volume (a:ℝ) b ∧
      (CertifiedCenteredIntegral.cell a b C D).Contains (∫ k in (a:ℝ)..(b:ℝ), f k) := by
  let m := ((a:ℝ)+(b:ℝ))/2
  let s := ((D.lo:ℝ)+(D.hi:ℝ))/2
  let e := ((D.hi:ℝ)-(D.lo:ℝ))/2
  let q := fun k : ℝ => s*(k-m)
  let t := fun k : ℝ => e*|k-m|
  let L := fun k : ℝ => ((C.lo:ℝ)+q k) - t k
  let U := fun k : ℝ => ((C.hi:ℝ)+q k) + t k
  have hab' : (a:ℝ) ≤ b := by exact_mod_cast hab
  have range (x : ℝ) (hx : x ∈ Icc (a:ℝ) b) : L x ≤ f x ∧ f x ≤ U x := by
    obtain ⟨v,hv,d,hd,he⟩ := hband x hx
    push_cast at he
    have h : |d-s| ≤ e := by
      unfold s e
      exact abs_le.mpr ⟨by linarith [hd.1],by linarith [hd.2]⟩
    have hb : |(d-s)*(x-m)| ≤ t x := by
      rw [abs_mul]
      exact mul_le_mul_of_nonneg_right h (abs_nonneg _)
    have eq : f x = v + q x + (d-s)*(x-m) := by rw [he]; dsimp [q,m]; ring
    rw [eq]
    dsimp only [L,U]
    exact ⟨by linarith [hv.1,(abs_le.mp hb).1],by linarith [hv.2,(abs_le.mp hb).2]⟩
  have cq : Continuous q := continuous_const.mul (continuous_id.sub continuous_const)
  have ct : Continuous t := continuous_const.mul (continuous_id.sub continuous_const).abs
  have cl : Continuous L := (continuous_const.add cq).sub ct
  have cu : Continuous U := (continuous_const.add cq).add ct
  have il : IntervalIntegrable L volume (a:ℝ) b := cl.intervalIntegrable _ _
  have iu : IntervalIntegrable U volume (a:ℝ) b := cu.intervalIntegrable _ _
  have hi : IntervalIntegrable f volume (a:ℝ) b := by
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le hab']
    apply integrable_of_le_of_le hf.aestronglyMeasurable ?_ ?_ il.1 iu.1
    · filter_upwards [ae_restrict_mem (μ := volume) (measurableSet_Ioc (a := (a:ℝ)) (b := (b:ℝ)))] with k hk
      exact (range k ⟨hk.1.le, hk.2⟩).1
    · filter_upwards [ae_restrict_mem (μ := volume) (measurableSet_Ioc (a := (a:ℝ)) (b := (b:ℝ)))] with k hk
      exact (range k ⟨hk.1.le, hk.2⟩).2
  have iq : (∫ k in (a:ℝ)..(b:ℝ), q k) = 0 := by
    change (∫ k in (a:ℝ)..(b:ℝ), s*(k-m)) = 0
    rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_sub
      (f := fun k : ℝ => k) (g := fun _ : ℝ => m) (continuous_id.intervalIntegrable _ _) intervalIntegrable_const,
      integral_id, intervalIntegral.integral_const]
    dsimp [m]; ring
  have it : (∫ k in (a:ℝ)..(b:ℝ), t k) = e*((b:ℝ)-(a:ℝ))^2/4 := by
    change (∫ k in (a:ℝ)..(b:ℝ), e*|k-m|) = _
    rw [intervalIntegral.integral_const_mul, CertifiedCenteredIntegral.absMoment hab']
    ring
  have lo := intervalIntegral.integral_mono_on hab' il hi (fun k hk => (range k hk).1)
  have hu := intervalIntegral.integral_mono_on hab' hi iu (fun k hk => (range k hk).2)
  change (∫ k in (a:ℝ)..(b:ℝ), ((C.lo:ℝ)+q k)-t k) ≤ _ at lo
  change _ ≤ ∫ k in (a:ℝ)..(b:ℝ), ((C.hi:ℝ)+q k)+t k at hu
  rw [intervalIntegral.integral_sub (f := fun k => (C.lo:ℝ)+q k) (g := t)
    ((continuous_const.add cq).intervalIntegrable _ _) (ct.intervalIntegrable _ _),
    intervalIntegral.integral_add (f := fun _ => (C.lo:ℝ)) (g := q)
    intervalIntegrable_const (cq.intervalIntegrable _ _),intervalIntegral.integral_const,iq,it] at lo
  rw [intervalIntegral.integral_add (f := fun k => (C.hi:ℝ)+q k) (g := t)
    ((continuous_const.add cq).intervalIntegrable _ _) (ct.intervalIntegrable _ _),
    intervalIntegral.integral_add (f := fun _ => (C.hi:ℝ)) (g := q)
    intervalIntegrable_const (cq.intervalIntegrable _ _),intervalIntegral.integral_const,iq,it] at hu
  refine ⟨hi,?_⟩
  unfold CertifiedCenteredIntegral.cell RationalInterval.Contains
  push_cast
  simp only [smul_eq_mul, add_zero, e] at lo hu
  constructor <;> linarith

theorem affine_band_sound {f : ℝ → ℝ} {a b : ℚ} {C D : RationalInterval}
    (hab : a ≤ b) (hf : Measurable f)
    (hband : ∀ k ∈ Icc (a:ℝ) b, ∃ v : ℝ, C.Contains v ∧ ∃ d : ℝ, D.Contains d ∧
      f k = v + d * (k - (((a+b)/2 : ℚ):ℝ))) :
    (CertifiedCenteredIntegral.cell a b C D).Contains (∫ k in (a:ℝ)..(b:ℝ), f k) :=
  (affine_band hab hf hband).2

#print axioms weighted_sound
#print axioms square_weight_sound
#print axioms affine_band_sound
end CertifiedWeightedIntegral
