import SashimiNumerics
import HzCertificate
import Physlib.Cosmology.FLRW.Age
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.SpecialFunctions.Arsinh
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

/-! Scalar cosmic age and elapsed time in the flat matter–Λ background.
The improper integral extends to infinite redshift. A finite lookback window
uses an exact difference of the same antiderivative. -/

namespace SashimiTime
open MeasureTheory Set Real Filter FloatLib.Numerics FloatLib.Numerics.Interval
open scoped Topology

/-- The decimal km/s/Mpc to Gyr⁻¹ factor used by ITAMAE. -/
def gyrFactor : ℚ := 10227121650537077 / 10000000000000000000

/-- The time integrand in the chosen inverse-Hubble units, including conversion factor c. -/
noncomputable def rate (H m c t : ℝ) : ℝ := 1 / ((1 + t) * HzCertificate.modelH H m t * c)

/-- Cosmic age with the same improper domain as ITAMAE's native implementation. -/
noncomputable def age (H m c z : ℝ) : ℝ := ∫ t in Ioi z, rate H m c t

/-- The elapsed time between two redshifts, including the lookback case a = 0. -/
noncomputable def elapsed (H m c a b : ℝ) : ℝ := ∫ t in a..b, rate H m c t

noncomputable def ageForm (D m z : ℝ) : ℝ :=
  if m = 1 then 2 / (3 * D) / sqrt (m * (1 + z)^3)
  else 2 / (3 * D * sqrt (1 - m)) * arsinh (sqrt (1 - m) / sqrt (m * (1 + z)^3))

theorem age_deriv {H m c t : ℝ} (hH : 0 < H) (hm : 0 < m) (hm₁ : m ≤ 1)
    (hc : 0 < c) (ht : 0 ≤ t) :
    HasDerivAt (ageForm (H * c) m) (-rate H m c t) t := by
  let v := 1 + t
  let w := sqrt (m * v^3)
  let s := sqrt (1 - m)
  let j := sqrt ((1 - m) + m * v^3)
  have hv : 0 < v := by dsimp [v]; linarith
  have hw : 0 < w := Real.sqrt_pos.2 (mul_pos hm (pow_pos hv _))
  have hj : 0 < j := Real.sqrt_pos.2 (add_pos_of_nonneg_of_pos (sub_nonneg.mpr hm₁)
    (mul_pos hm (pow_pos hv _)))
  have hw2 : w^2 = m * v^3 := Real.sq_sqrt (mul_nonneg hm.le (pow_nonneg hv.le _))
  have hder : HasDerivAt (fun u : ℝ => sqrt (m * (1 + u)^3))
      (m * (3 * v^2) / (2 * w)) t := by
    have hd := (((hasDerivAt_id t).const_add 1).pow 3).const_mul m
    simpa [v, w] using hd.sqrt (ne_of_gt (mul_pos hm (pow_pos hv 3)))
  have hrate : rate H m c t = 1 / (v * (H * j) * c) := by
    rw [rate, HzCertificate.modelH_eq H m t ht]
  by_cases h1 : m = 1
  · have hej : j = w := by simp [j, w, h1]
    unfold ageForm
    simp only [ite_eq_left h1]
    have hd := (hasDerivAt_const t (2 / (3 * (H * c)))).div hder hw.ne'
    apply hd.congr_deriv
    change (0*w - 2/(3*(H*c)) * (m*(3*v^2)/(2*w))) / w^2 = _
    rw [hrate, hej]
    field_simp
    nlinarith [hw2]
  · have hs : 0 < s := Real.sqrt_pos.2 (sub_pos.mpr (lt_of_le_of_ne hm₁ h1))
    have hs2 : s^2 = 1 - m := Real.sq_sqrt (sub_nonneg.mpr hm₁)
    have hrw : sqrt (1 + (s / w)^2) = j / w := by
      have hratio : 1 + (s / w)^2 = ((1 - m) + m * v^3) / w^2 := by
        field_simp
        nlinarith [hw2, hs2]
      rw [hratio, Real.sqrt_div (add_nonneg (sub_nonneg.mpr hm₁) (mul_nonneg hm.le (pow_nonneg hv.le _))),
        Real.sqrt_sq hw.le]
    have hd := (((hasDerivAt_const t s).div hder hw.ne').arsinh).const_mul
      (2 / (3 * (H * c) * s))
    unfold ageForm
    simp only [ite_eq_right h1]
    apply hd.congr_deriv
    change 2/(3*(H*c)*s) * ((sqrt (1+(s/w)^2))⁻¹ * ((0*w - s*(m*(3*v^2)/(2*w))) / w^2)) = _
    rw [hrate, hrw]
    field_simp
    nlinarith [hw2]

private theorem age_limit {D m : ℝ} (hm : 0 < m) :
    Tendsto (ageForm D m) atTop (𝓝 0) := by
  have hv : Tendsto (fun x : ℝ => 1 + x) atTop atTop := tendsto_atTop_add_const_left _ _ tendsto_id
  have hp := (tendsto_pow_atTop (α := ℝ) (by decide : 3 ≠ 0)).comp hv
  have hw := Real.tendsto_sqrt_atTop.comp (hp.const_mul_atTop hm)
  by_cases h1 : m = 1
  · unfold ageForm
    simpa only [ite_eq_left h1, Function.comp_def] using
      (tendsto_const_nhds (x := 2/(3*D))).div_atTop hw
  · have ha := (Real.continuous_arsinh.tendsto 0).comp
      ((tendsto_const_nhds (x := sqrt (1-m))).div_atTop hw)
    unfold ageForm
    simpa [h1, Function.comp_def] using
      (tendsto_const_nhds (x := 2 / (3*D*sqrt (1-m)))).mul ha

private theorem rate_pos {H m c t : ℝ} (hH : 0 < H) (hm : 0 < m) (hc : 0 < c)
    (ht : 0 ≤ t) : 0 < rate H m c t := by
  unfold rate
  have h := HzCertificate.modelH_pos H m t hH hm.le ht
  positivity

theorem age_eq {H m c z : ℝ} (hH : 0 < H) (hm : 0 < m) (hm₁ : m ≤ 1)
    (hc : 0 < c) (hz : 0 ≤ z) :
    age H m c z = ageForm (H * c) m z := by
  have hder : ∀ t ∈ Ici z, HasDerivAt (ageForm (H*c) m) (-rate H m c t) t :=
    fun t ht => age_deriv hH hm hm₁ hc (hz.trans ht)
  have hn : ∀ t ∈ Ioi z, -rate H m c t ≤ 0 :=
    fun t ht => neg_nonpos.mpr (rate_pos hH hm hc (hz.trans ht.le)).le
  have h := integral_Ioi_of_hasDerivAt_of_nonpos' hder hn (age_limit hm)
  rw [integral_neg] at h
  simpa [age] using congrArg Neg.neg h

theorem elapsed_eq {H m c a b : ℝ} (hH : 0 < H) (hm : 0 < m) (hm₁ : m ≤ 1)
    (hc : 0 < c) (ha : 0 ≤ a) (hab : a ≤ b) :
    elapsed H m c a b = age H m c a - age H m c b := by
  have hf : ∀ t ∈ Ici a, HasDerivAt (ageForm (H*c) m) (-rate H m c t) t :=
    fun t ht => age_deriv hH hm hm₁ hc (ha.trans ht)
  have hi := integrableOn_Ioi_deriv_of_nonpos' hf (fun x hx => neg_nonpos.mpr
    (rate_pos hH hm hc (ha.trans hx.le)).le) (age_limit hm)
  have hi' : IntegrableOn (rate H m c) (Ioi a) := by simpa using hi.neg
  exact (intervalIntegral.integral_Ioi_sub_Ioi hi' hab).symm

/-- Closed expression in the redshift, with rational radicand formed before
rounding. FloatLib evaluates arsinh through its defining real log and sqrt. -/
def ageExpr (H m c z : ℚ) : Expr :=
  let root : Expr := .unary .sqrt (.const (m * (1 + z)^3))
  if m = 1 then .binary .div (.const (2 / (3 * (H * c)))) root
  else
    let s : Expr := .unary .sqrt (.const (1 - m))
    let u : Expr := .binary .div s root
    let ash : Expr := .unary .log (u + .unary .sqrt (.const 1 + .unary (.pow 2) u))
    .binary .mul (.binary .div (.const (2 / (3 * (H * c)))) s) ash

theorem ageExpr_eq (H m c z : ℚ) :
    (ageExpr H m c z).eval (fun _ => 0) = ageForm ((H : ℝ)*c) m z := by
  by_cases hm : m = 1
  · subst m
    simp [ageExpr, ageForm, Expr.eval, UnaryOp.eval, BinaryOp.eval]
  · have hm' : (m : ℝ) ≠ 1 := by exact_mod_cast hm
    simp [ageExpr, ageForm, hm, hm', Expr.eval, UnaryOp.eval, BinaryOp.eval, Real.arsinh]
    ring_nf
    simp

/-- Return cosmic age at a for none, and elapsed time from a to b for some b.
The domain excludes zero matter density, for which the infinite-age formula diverges. -/
private def positiveEnclose (H m c a : ℚ) (until? : Option ℚ) (tol : ℚ)
    (bits : Nat := 64) (degree : Nat := 24) : Option RationalInterval := do
  guard (0 < H ∧ 0 < m ∧ m ≤ 1 ∧ 0 < c ∧ 0 ≤ a)
  match until? with
  | none => SashimiNumerics.closed? (ageExpr H m c a) bits tol degree
  | some b =>
    guard (a ≤ b)
    SashimiNumerics.closed? (.binary .sub (ageExpr H m c a) (ageExpr H m c b)) bits tol degree

/-- A returned enclosure satisfies its width budget and contains the precise
finite or improper integral built directly from Physlib's Hubble prediction. -/
private theorem positiveEnclose_sound {H m c a tol : ℚ} {until? : Option ℚ} {bits degree : Nat}
    {I : RationalInterval} (h : positiveEnclose H m c a until? tol bits degree = some I) :
    I.Contains (match until? with
      | none => age H m c a
      | some b => elapsed H m c a b) ∧ I.hi - I.lo ≤ tol := by
  by_cases hv : 0 < H ∧ 0 < m ∧ m ≤ 1 ∧ 0 < c ∧ 0 ≤ a
  · obtain ⟨hH, hm, hm₁, hc, ha⟩ := hv
    have hH' : 0 < (H : ℝ) := by exact_mod_cast hH
    have hm' : 0 < (m : ℝ) := by exact_mod_cast hm
    have hm₁' : (m : ℝ) ≤ 1 := by exact_mod_cast hm₁
    have hc' : 0 < (c : ℝ) := by exact_mod_cast hc
    have ha' : 0 ≤ (a : ℝ) := by exact_mod_cast ha
    cases until? with
    | none =>
      simp [positiveEnclose, hH, hm, hm₁, hc, ha] at h
      simpa [ageExpr_eq, age_eq hH' hm' hm₁' hc' ha'] using SashimiNumerics.closed_sound h
    | some b =>
      by_cases hab : a ≤ b
      · have hab' : (a : ℝ) ≤ b := by exact_mod_cast hab
        simp [positiveEnclose, hH, hm, hm₁, hc, ha, hab] at h
        obtain ⟨hx, hw⟩ := SashimiNumerics.closed_sound h
        refine ⟨?_, hw⟩
        change I.Contains ((ageExpr H m c a).eval _ - (ageExpr H m c b).eval _) at hx
        simpa only [ageExpr_eq, elapsed_eq hH' hm' hm₁' hc' ha' hab',
          age_eq hH' hm' hm₁' hc' ha', age_eq hH' hm' hm₁' hc' (ha'.trans hab')] using hx
      · simp [positiveEnclose, hH, hm, hm₁, hc, ha, hab] at h
  · simp [positiveEnclose, hv] at h

#print axioms age_eq
#print axioms elapsed_eq

/-- The time from scale factor zero to x, in inverse Hubble units. This
representation extends the age continuously to x = 0. -/
private noncomputable def scaleForm (m x : ℝ) : ℝ :=
  if m = 1 then 2/3 * sqrt (x^3 / m)
  else 2/(3*sqrt (1-m)) * arsinh (sqrt (1-m) * sqrt (x^3/m))

private theorem scale_eq {m x : ℝ} (hm : 0 < m) (hx : 0 < x) :
    scaleForm m x = ageForm 1 m (x⁻¹-1) := by
  have h1 : 1 + (x⁻¹-1) = x⁻¹ := by ring
  have he : sqrt (x^3 / m) = (sqrt (m * x⁻¹ ^ 3))⁻¹ := by
    rw [← Real.sqrt_inv]
    congr 1
    field_simp
  unfold scaleForm ageForm
  rw [h1, he]
  split <;> simp [div_eq_mul_inv]

/-- At the reference redshift the redshift-domain age equals Physlib's
scale-factor age integral for exactly the same reduced Hubble function.
H and c are combined here so both integrals use the same time units. -/
theorem age_zero_eq_physlib {H m c : ℝ} (hH : 0 < H) (hm : 0 < m) (hm₁ : m ≤ 1)
    (hc : 0 < c) :
    age H m c 0 = Cosmology.FLRW.FriedmannEquation.age (H*c)
      (Cosmology.FLRW.FriedmannEquation.reducedHubble (1-m) m 0 0) := by
  let E := Cosmology.FLRW.FriedmannEquation.reducedHubble (1-m) m 0 0
  let f := fun x : ℝ => 1/(x * E x)
  have hcont : Continuous (scaleForm m) := by
    have hu : Continuous (fun x : ℝ => sqrt (x^3/m)) :=
      Real.continuous_sqrt.comp ((continuous_id.pow 3).div_const m)
    unfold scaleForm
    by_cases hh : m = 1
    · simp only [ite_eq_left hh]
      exact continuous_const.mul hu
    · simp only [ite_eq_right hh]
      exact continuous_const.mul (Real.continuous_arsinh.comp (continuous_const.mul hu))
  have hd (x : ℝ) (hx : x ∈ Ioo (0 : ℝ) 1) : HasDerivAt (scaleForm m) (f x) x := by
    have hx0 := hx.1
    have hz : 0 ≤ x⁻¹ - 1 := sub_nonneg.mpr ((one_le_inv₀ hx0).mpr hx.2.le)
    have hd := (age_deriv (H := 1) (m := m) (c := 1) zero_lt_one hm hm₁ zero_lt_one hz).comp x
      (((hasDerivAt_id x).inv hx0.ne').sub_const 1)
    have hmatch : (fun u : ℝ => ageForm (1*1) m (u⁻¹-1)) =ᶠ[𝓝 x] scaleForm m := by
      filter_upwards [lt_mem_nhds hx0] with u hu
      simpa using (scale_eq hm hu).symm
    apply (hd.congr_of_eventuallyEq hmatch.symm).congr_deriv
    have hsum : 1 + (x⁻¹-1) = x⁻¹ := by ring
    simp [f, rate, HzCertificate.modelH, hsum, E, div_eq_mul_inv]
    field_simp
  have hpos (x : ℝ) (hx : x ∈ Ioo (0 : ℝ) 1) : 0 ≤ f x := by
    unfold f E Cosmology.FLRW.FriedmannEquation.reducedHubble
    exact one_div_nonneg.mpr (mul_nonneg hx.1.le (Real.sqrt_nonneg _))
  have hInt : IntervalIntegrable f volume 0 1 :=
    intervalIntegral.intervalIntegrable_deriv_of_nonneg hcont.continuousOn
      (by simpa using hd) (by simpa using hpos)
  have hFTC := intervalIntegral.integral_eq_sub_of_hasDerivAt_of_le zero_le_one
    hcont.continuousOn hd hInt
  have hz0 : scaleForm m 0 = 0 := by
    unfold scaleForm
    split <;> simp
  have hz1 : scaleForm m 1 = ageForm 1 m 0 := by simpa using scale_eq hm zero_lt_one
  rw [hz0, hz1, sub_zero] at hFTC
  rw [age_eq hH hm hm₁ hc le_rfl]
  change ageForm (H*c) m 0 = 1/(H*c) * ∫ x in (0:ℝ)..1, f x
  rw [hFTC]
  unfold ageForm
  split <;> ring

#print axioms age_zero_eq_physlib

/-- In pure Λ the finite time window has a logarithmic antiderivative. -/
theorem elapsed_zero {H c a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) :
    elapsed H 0 c a b = (1/(H*c)) * log ((1+b)/(1+a)) := by
  have hpos : 0 < 1+a := by linarith
  have h0 : (0:ℝ) ∉ uIcc (1+a) (1+b) := by
    rw [uIcc_of_le (by linarith)]
    exact fun h => (not_le.mpr hpos) h.1
  unfold elapsed
  calc
    _ = ∫ t in a..b, (1/(H*c)) * (1+t)⁻¹ := by
      apply intervalIntegral.integral_congr
      intro t ht
      rw [uIcc_of_le hab] at ht
      rw [rate, HzCertificate.modelH_eq H 0 t (ha.trans ht.1)]
      simp [div_eq_mul_inv]
      ring
    _ = _ := by
      rw [intervalIntegral.integral_const_mul,
        intervalIntegral.integral_comp_add_left (fun t : ℝ => t⁻¹) 1, integral_inv h0]

private def lambdaWindow (H c a b : ℚ) : Expr :=
  .binary .mul (.const (1/(H*c))) (.unary .log (.const ((1+b)/(1+a))))

/-- An infinite age needs positive matter density. A finite pure-Λ window can
be integrated directly, even though the two infinite-age terms would diverge. -/
def enclose (H m c a : ℚ) (until? : Option ℚ) (tol : ℚ)
    (bits : Nat := 64) (degree : Nat := 24) : Option RationalInterval :=
  if m = 0 then do
    let b ← until?
    guard (0 < H ∧ 0 < c ∧ 0 ≤ a ∧ a ≤ b)
    SashimiNumerics.closed? (lambdaWindow H c a b) bits tol degree
  else positiveEnclose H m c a until? tol bits degree

theorem enclose_sound {H m c a tol : ℚ} {until? : Option ℚ} {bits degree : Nat}
    {I : RationalInterval} (h : enclose H m c a until? tol bits degree = some I) :
    I.Contains (match until? with
      | none => age H m c a
      | some b => elapsed H m c a b) ∧ I.hi - I.lo ≤ tol := by
  unfold enclose at h
  split at h
  · rename_i hm
    subst m
    cases until? with
    | none => contradiction
    | some b =>
      by_cases hv : 0 < H ∧ 0 < c ∧ 0 ≤ a ∧ a ≤ b
      · simp [hv] at h
        obtain ⟨hc, hw⟩ := SashimiNumerics.closed_sound h
        refine ⟨?_, hw⟩
        have ha' : 0 ≤ (a : ℝ) := by exact_mod_cast hv.2.2.1
        have hab' : (a : ℝ) ≤ b := by exact_mod_cast hv.2.2.2
        simpa [lambdaWindow, Expr.eval, BinaryOp.eval, UnaryOp.eval, elapsed_zero ha' hab'] using hc
      · simp [hv] at h
  · cases until? <;> exact positiveEnclose_sound h

#print axioms enclose_sound

set_option maxRecDepth 20000 in
set_option maxHeartbeats 8000000 in
/-- Reference age in Gyr for ITAMAE’s default parameters, checked end-to-end. -/
theorem native_reference_age : RationalInterval.Contains
    ⟨254495609656065729025 / 18446744073709551616, 254495609657592531779 / 18446744073709551616⟩
    (age (337/5 : ℚ) (63/200 : ℚ) gyrFactor (0 : ℚ)) := by
  exact (enclose_sound (H := 337/5) (m := 63/200) (c := gyrFactor) (a := 0)
    (until? := none) (tol := 1/100000000) (bits := 64) (degree := 24) (by decide +kernel)).1

#print axioms native_reference_age

end SashimiTime
