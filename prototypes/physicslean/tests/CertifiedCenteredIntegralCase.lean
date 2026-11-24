import CertifiedCenteredIntegral
import HzAudit

/-! Standalone cells for a decreasing affine function and for a quadratic.
The quadratic integral is independently checked with its antiderivative. -/
namespace CertifiedCenteredIntegralCase
open Set MeasureTheory FloatLib.Numerics CertifiedCenteredIntegral

theorem affine_exact : (∫ x : ℝ in (-2)..4, 5 - 3*x) = 12 := by
  have h := cell_sound (f := fun x : ℝ => 5-3*x) (d := fun _ => -3)
    (a := -2) (b := 4) (F := ⟨2,2⟩) (D := ⟨-3,-3⟩) (by norm_num)
    (by norm_num [RationalInterval.Contains]) (by
      intro x _
      have h := (hasDerivAt_const x (5:ℝ)).sub ((hasDerivAt_id x).const_mul 3)
      simpa [Pi.sub_def] using! h.hasDerivWithinAt (s := Icc (-2) 4)) (by
      intro x _
      norm_num [RationalInterval.Contains])
  have ends : cell (-2) 4 ⟨2,2⟩ ⟨-3,-3⟩ = (⟨12,12⟩ : RationalInterval) := by decide +kernel
  rw [ends] at h
  simpa [RationalInterval.Contains, le_antisymm_iff] using h.2.antisymm h.1

theorem quadratic_bounds : (⟨6,10⟩ : RationalInterval).Contains (∫ x : ℝ in 1..3, x^2) := by
  have h := cell_sound (f := fun x : ℝ => x^2) (d := fun x => 2*x)
    (a := 1) (b := 3) (F := ⟨4,4⟩) (D := ⟨2,6⟩) (by norm_num)
    (by norm_num [RationalInterval.Contains]) (by
      intro x _
      simpa [Pi.pow_def] using! ((hasDerivAt_id x).pow 2).hasDerivWithinAt (s := Icc 1 3)) (by
      intro x hx
      norm_num at hx
      change (2:ℝ) ≤ 2*x ∧ 2*x ≤ 6
      constructor <;> linarith [hx.1,hx.2])
  have ends : cell 1 3 ⟨4,4⟩ ⟨2,6⟩ = (⟨6,10⟩ : RationalInterval) := by decide +kernel
  rwa [ends] at h

theorem quadratic_exact : (∫ x : ℝ in 1..3, x^2) = 26/3 := by
  rw [integral_pow]; norm_num

#print axioms affine_exact
#print axioms quadratic_bounds
#print axioms quadratic_exact
end CertifiedCenteredIntegralCase

run_cmd PhysicsLeanAudit.currentWitnesses
