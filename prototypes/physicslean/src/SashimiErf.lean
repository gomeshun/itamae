import SashimiNumerics
import Mathlib.Analysis.SpecialFunctions.Gaussian.GaussianIntegral
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-! Gaussian error function. Integrate the exponential Taylor polynomial
exactly over rational endpoints and bound its remainder on the entire interval.
Saturate only after proving monotonicity and the normalized Gaussian integral. -/
namespace SashimiErf
open FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory Set Real Finset

noncomputable def erf (x : ℝ) : ℝ := 2 / sqrt π * ∫ t in (0:ℝ)..x, exp (-(t^2))

private theorem eint (a b : ℝ) :
    IntervalIntegrable (fun t : ℝ => exp (-(t^2))) volume a b := by
  apply Continuous.intervalIntegrable
  fun_prop

theorem erf_mono : Monotone erf := by
  intro a b hab
  have h := intervalIntegral.integral_add_adjacent_intervals (eint 0 a) (eint a b)
  have hg := intervalIntegral.integral_nonneg_of_forall (μ := volume) hab (fun t : ℝ => (exp_pos (-(t^2))).le)
  exact mul_le_mul_of_nonneg_left (by linarith) (by positivity)

theorem erf_neg (x : ℝ) : erf (-x) = -erf x := by
  have h := intervalIntegral.integral_comp_neg (a := (0:ℝ)) (b := x)
    (fun t : ℝ => exp (-(t^2)))
  simp only [neg_zero, neg_sq] at h
  unfold erf
  rw [intervalIntegral.integral_symm (a := -x) (b := 0), ← h]
  ring

theorem erf_le_one (x : ℝ) : erf x ≤ 1 := by
  have hall : Integrable (fun t : ℝ => exp (-(t^2))) := by
    simpa using integrable_exp_neg_mul_sq (b := (1:ℝ)) zero_lt_one
  have h := intervalIntegral.integral_interval_add_Ioi (a := (0:ℝ)) (b := x)
    hall.integrableOn hall.integrableOn
  have hfull : (∫ t in Ioi (0:ℝ), exp (-(t^2))) = sqrt π / 2 := by
    simpa using integral_gaussian_Ioi 1
  have htail : 0 ≤ (∫ t in Ioi x, exp (-(t^2))) := MeasureTheory.integral_nonneg (fun _ => (exp_pos _).le)
  rw [hfull] at h
  unfold erf
  have hs : 0 < sqrt π := Real.sqrt_pos.2 Real.pi_pos
  calc
    _ ≤ 2 / sqrt π * (sqrt π / 2) := mul_le_mul_of_nonneg_left (by linarith) (by positivity)
    _ = 1 := by field_simp

private noncomputable def term (i : ℕ) (x : ℝ) : ℝ := (-(x^2))^i / i.factorial

/-- Primitive of the first n terms at a nonnegative rational endpoint. -/
def sumAt (x : ℚ) (n : Nat) : ℚ :=
  ∑ i ∈ range n, (-1:ℚ)^i / i.factorial * (x^(2*i+1)/(2*i+1))

private theorem sum_integral (x : ℚ) (n : Nat) :
    (∫ t in (0:ℝ)..x, ∑ i ∈ range n, term i t) = (sumAt x n : ℝ) := by
  have hint (i : Nat) : IntervalIntegrable (term i) volume 0 x := by
    apply Continuous.intervalIntegrable; unfold term; fun_prop
  rw [intervalIntegral.integral_finsetSum (fun i _ => hint i)]
  unfold sumAt
  push_cast
  apply sum_congr rfl
  intro i hi
  have hfun : term i = fun t : ℝ => ((-1:ℝ)^i/i.factorial) * t^(2*i) := by
    funext t
    unfold term
    rw [show -(t^2) = (-1:ℝ) * t^2 by ring, mul_pow, pow_mul]
    ring
  rw [hfun, intervalIntegral.integral_const_mul, integral_pow]
  push_cast
  simp

private theorem exp_error (n : Nat) (x : ℝ) (hx : |x|/(n+1) ≤ 1/2) :
    |exp x - ∑ i ∈ range n, x^i/i.factorial| ≤ |x|^n/n.factorial*2 := by
  have hc : ‖(x:ℂ)‖/(n.succ : ℝ) ≤ 1/2 := by simpa using hx
  convert Complex.exp_bound' hc using 1 <;> norm_cast

def radius (x : ℚ) (n : Nat) : ℚ := (x^2)^n/n.factorial*2*x

private theorem integral_error {x : ℚ} (hx : 0 ≤ x) (n : Nat)
    (hn : 2*x^2 ≤ (n:ℚ)+1) :
    (RationalInterval.around (sumAt x n) (radius x n)).Contains
      (∫ t in (0:ℝ)..x, exp (-(t^2))) := by
  apply RationalInterval.contains_around
  rw [← sum_integral]
  have hxr : (0:ℝ) ≤ x := by exact_mod_cast hx
  have hn' : 2*(x:ℝ)^2 ≤ n+1 := by exact_mod_cast hn
  have hpoly : IntervalIntegrable (fun t => ∑ i ∈ range n, term i t) volume 0 x := by
    apply Continuous.intervalIntegrable
    unfold term
    fun_prop
  rw [← intervalIntegral.integral_sub (eint 0 x) hpoly]
  have hle := intervalIntegral.norm_integral_le_of_norm_le_const (a := (0:ℝ)) (b := (x:ℝ))
    (f := fun t => exp (-(t^2)) - ∑ i ∈ range n, term i t) (C := ((x:ℝ)^2)^n/n.factorial*2) (by
      intro t ht
      rw [uIoc_of_le hxr] at ht
      have ht2 := pow_le_pow_left₀ ht.1.le ht.2 2
      have hd : |-(t^2)|/((n:ℝ)+1) ≤ 1/2 := by
        rw [abs_neg, abs_of_nonneg (sq_nonneg _)]
        apply (div_le_iff₀ (by positivity : 0 < (n:ℝ)+1)).mpr
        linarith
      have h := exp_error n (-(t^2)) hd
      simp only [Real.norm_eq_abs, term, abs_neg, abs_of_nonneg (sq_nonneg t)] at h ⊢
      exact h.trans (by gcongr))
  simpa [radius, abs_of_nonneg hxr, mul_comm] using hle

private def normalization : Expr := .binary .mul
  (.binary .div (.const 2) (.unary .sqrt (.binary .mul (.const 4) (.unary .atan (.const 1)))))
  (.var 0)

/-- Only sum up to five. For larger x the lower bound at five and the Gaussian
integral upper bound one avoid huge cancelling terms. -/
private def positive? (x : ℚ) (bits degree : Nat) : Option RationalInterval := do
  let y := min x 5
  let n := max degree (⌈4*y^2⌉₊ + 16)
  if 2*y^2 ≤ (n:ℚ)+1 then
    let J := RationalInterval.around (sumAt y n) (radius y n)
    let I ← normalization.gridEnclosureWith? [J] bits degree
    if x ≤ 5 then return I else return ⟨I.lo, 1⟩
  else none

private theorem positive_sound {x : ℚ} {bits degree : Nat} {I : RationalInterval}
    (hx : 0 ≤ x) (h : positive? x bits degree = some I) : I.Contains (erf x) := by
  let y : ℚ := min x 5
  let n := max degree (⌈4*y^2⌉₊ + 16)
  let J := RationalInterval.around (sumAt y n) (radius y n)
  by_cases hn : 2*y^2 ≤ (n:ℚ)+1
  · have hJ : J.Contains (∫ t in (0:ℝ)..y, exp (-(t^2))) :=
      integral_error (le_min hx (by norm_num)) n hn
    have hbind : (normalization.gridEnclosureWith? [J] bits degree).bind
        (fun K => if x ≤ 5 then some K else some ⟨K.lo, 1⟩) = some I := by
      change (if 2*y^2 ≤ (n:ℚ)+1 then _ else none) = some I at h
      rwa [ite_eq_left hn] at h
    obtain ⟨K, hk, hi⟩ := Option.bind_eq_some_iff.mp hbind
    have hk' := Expr.contains_gridEnclosureWith?
      (values := fun _ => ∫ t in (0:ℝ)..y, exp (-(t^2))) (h := hk) (by
        intro i U hU
        cases i with
        | zero => simp at hU; subst U; exact hJ
        | succ i => simp at hU)
    have hK : K.Contains (erf y) := by
      change K.Contains (2 / sqrt (4 * Real.arctan 1) * _) at hk'
      rw [Real.arctan_one, show 4*(π/4)=π by ring] at hk'
      exact hk'
    split at hi
    · cases hi
      simpa [y, min_eq_left ‹x ≤ 5›] using hK
    · cases hi
      exact ⟨hK.1.trans (erf_mono (by exact_mod_cast min_le_left x (5:ℚ))), by simpa using erf_le_one x⟩
  · change (if 2*y^2 ≤ (n:ℚ)+1 then _ else none) = some I at h
    simp [hn] at h

private def point? (x : ℚ) (bits degree : Nat) : Option RationalInterval :=
  if 0 ≤ x then positive? x bits degree else (positive? (-x) bits degree).map RationalInterval.neg

private theorem point_sound {x : ℚ} {bits degree : Nat} {I : RationalInterval}
    (h : point? x bits degree = some I) : I.Contains (erf x) := by
  unfold point? at h
  split at h
  · exact positive_sound ‹0 ≤ x› h
  · obtain ⟨J, hj, rfl⟩ := Option.map_eq_some_iff.mp h
    simpa [erf_neg] using RationalInterval.contains_neg (positive_sound (by linarith : 0 ≤ -x) hj)

/-- Enclose erf throughout an input interval by monotonicity. Endpoints and
all Taylor sums are rational. Elementary normalization uses FloatLib's grid. -/
def erf? (input : RationalInterval) (bits degree : Nat) : Option RationalInterval := do
  if input.lo ≤ input.hi then
    let L ← point? input.lo bits degree
    let U ← point? input.hi bits degree
    return ⟨L.lo, U.hi⟩
  else none

theorem erf_sound {input I : RationalInterval} {x : ℝ} {bits degree : Nat}
    (hx : input.Contains x) (h : erf? input bits degree = some I) : I.Contains (erf x) := by
  have hv := RationalInterval.lo_le_hi hx
  unfold erf? at h
  rw [ite_eq_left hv] at h
  obtain ⟨L, hL, hrest⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨U, hU, hI⟩ := Option.bind_eq_some_iff.mp hrest
  cases hI
  exact ⟨(point_sound hL).1.trans (erf_mono hx.1), (erf_mono hx.2).trans (point_sound hU).2⟩

#print axioms erf_sound

end SashimiErf
