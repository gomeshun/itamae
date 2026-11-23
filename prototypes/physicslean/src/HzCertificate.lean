import Physlib.Cosmology.FLRW.ReducedHubble
import Mathlib.Tactic.NormNum

/-!
# Certified H(z) for flat matter + Lambda FLRW cosmology

This layer connects the Physlib `reducedHubble` definition (the exact real-valued
specification) with rational-square certificate checking.  The physical theorem
relating `reducedHubble` to Friedmann solutions is already in Physlib; we use its
reduced Hubble definition but do not instantiate a full spacetime solution here.

The numeric witness is a lower/upper rational enclosure for a nonnegative real
whose *square* is a known rational number.
-/

namespace HzCertificate

open Cosmology.FLRW.FriedmannEquation

/-- Physical interpretation: `H₀` in chosen Hubble units, `Ωm` dimensionless,
    `z` dimensionless. Flat ΛCDM, no radiation, ΩΛ=1-Ωm. -/
noncomputable def modelH (H₀ Ωm z : ℝ) : ℝ :=
  H₀ * reducedHubble (1 - Ωm) Ωm 0 0 (1 / (1 + z))

/-- Convert Physlib's inverse-scale-factor power into the standard redshift form. -/
private theorem inverseScale_pow (z : ℝ) :
    (1 / (1 + z) : ℝ) ^ (-3 : ℝ) = (1 + z) ^ 3 := by
  rw [show (-3 : ℝ) = -((3 : ℕ) : ℝ) by norm_num,
      Real.rpow_neg_eq_inv_rpow, Real.rpow_natCast]
  simp [one_div]

/-- Algebraic identity for the total real operations used by Physlib and FloatLib. -/
theorem modelH_eq_unconditional (H₀ Ωm z : ℝ) :
    modelH H₀ Ωm z = H₀ * Real.sqrt ((1 - Ωm) + Ωm * (1 + z) ^ 3) := by
  unfold modelH reducedHubble
  simp only [zero_mul, add_zero, inverseScale_pow]

-- Keep the original named redshift hypothesis for callers of the physical-domain API.
set_option linter.unusedVariables false in
/-- Redshift form on the physical redshift domain. -/
theorem modelH_eq (H₀ Ωm z : ℝ) (hz : 0 ≤ z) :
    modelH H₀ Ωm z = H₀ * Real.sqrt ((1 - Ωm) + Ωm * (1 + z) ^ 3) :=
  modelH_eq_unconditional H₀ Ωm z

/-- Nonnegativity of the expanding-model prediction. -/
theorem modelH_nonneg (H₀ Ωm z : ℝ) (hH : 0 ≤ H₀) (hz : 0 ≤ z) :
    0 ≤ modelH H₀ Ωm z := by
  rw [modelH_eq H₀ Ωm z hz]
  exact mul_nonneg hH (Real.sqrt_nonneg _)

/-- For positive H₀ and nonnegative matter density and redshift, H(z) is strictly positive. -/
theorem modelH_pos (H₀ Ωm z : ℝ) (hH : 0 < H₀) (hΩ : 0 ≤ Ωm) (hz : 0 ≤ z) :
    0 < modelH H₀ Ωm z := by
  rw [modelH_eq H₀ Ωm z hz]
  apply mul_pos hH
  apply Real.sqrt_pos.2
  have hc : (1 : ℝ) ≤ (1 + z) ^ 3 := one_le_pow₀ (le_add_of_nonneg_right hz)
  calc
    (0 : ℝ) < 1 := zero_lt_one
    _ = (1 - Ωm) + Ωm * 1 := by simp
    _ ≤ _ := add_le_add le_rfl (mul_le_mul_of_nonneg_left hc hΩ)

/-- The square of H(z) is a polynomial expression over the input reals. -/
theorem modelH_sq (H₀ Ωm z : ℝ) (hz : 0 ≤ z)
    (hr : 0 ≤ (1 - Ωm) + Ωm * (1 + z) ^ 3) :
    modelH H₀ Ωm z ^ 2 = H₀ ^ 2 * ((1 - Ωm) + Ωm * (1 + z) ^ 3) := by
  rw [modelH_eq H₀ Ωm z hz, mul_pow, Real.sq_sqrt hr]

/-- Certificate checker soundness: if the exact rational-square inequalities hold,
    the bounds enclose the Physlib-defined H(z) in real arithmetic. -/
theorem modelH_bounds (H₀ Ωm z lower upper : ℝ)
    (hz : 0 ≤ z) (hH : 0 ≤ H₀)
    (hr : 0 ≤ (1 - Ωm) + Ωm * (1 + z) ^ 3)
    (hLower : 0 ≤ lower) (hUpper : 0 ≤ upper)
    (hSqLower : lower ^ 2 ≤ H₀ ^ 2 * ((1 - Ωm) + Ωm * (1 + z) ^ 3))
    (hSqUpper : H₀ ^ 2 * ((1 - Ωm) + Ωm * (1 + z) ^ 3) ≤ upper ^ 2) :
    lower ≤ modelH H₀ Ωm z ∧ modelH H₀ Ωm z ≤ upper := by
  have hs := modelH_sq H₀ Ωm z hz hr
  have hn := modelH_nonneg H₀ Ωm z hH hz
  rw [← hs] at hSqLower hSqUpper
  exact ⟨(sq_le_sq₀ hLower hn).mp hSqLower, (sq_le_sq₀ hn hUpper).mp hSqUpper⟩

#print axioms modelH_eq_unconditional
#print axioms modelH_eq
#print axioms modelH_nonneg
#print axioms modelH_pos
#print axioms modelH_sq
#print axioms modelH_bounds

end HzCertificate
