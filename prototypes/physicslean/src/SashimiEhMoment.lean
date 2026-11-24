import CertifiedWeightedIntegral
import SashimiEhDerivativeBounds

/-! Integrate the positive-wavenumber EH98 spectrum weight on a rational cell
by extracting its exact quadratic moment. The residual range reuses the same
transfer function and coefficient boxes as the original expression. -/
namespace SashimiEhMoment

open FloatLib.Numerics SashimiEhSweep SashimiEisensteinIntegral
open CertifiedWeightedIntegral Set MeasureTheory

/-- Only the exponent in the range expression changes. EH transfer coefficients
depend on the cosmological parameters other than the primordial tilt. -/
private def residual (c : Cosmology) : Cosmology := {c with ns := c.ns - 2}

def integral? (c : Cosmology) (boxes : List RationalInterval) (a b : ℚ) (bits degree : Nat) :
    Option RationalInterval :=
  (weight? (residual c) boxes a b bits degree).map (·.scaleNonnegative (squareMoment a b))

private theorem weight_eq {c : Cosmology} {k : ℝ} (hk : 0 < k) :
    weightReal c k = k^2 * weightReal (residual c) k := by
  unfold weightReal
  rw [tilt_eq c hk, tilt_eq (residual c) hk]
  have h : (2 + ((residual c).ns : ℝ)) = (c.ns : ℝ) := by simp [residual]
  rw [h, Real.rpow_add hk, Real.rpow_two]
  change _ = k^2 * (_ * (SashimiEisensteinHu.transferReal (coefficients c) k)^2)
  ring

theorem integral_sound {c boxes cb cd a b bits degree A}
    (hb : boxes? c cb cd = some boxes) (h : integral? c boxes a b bits degree = some A) :
    A.Contains (∫ k in (a:ℝ)..(b:ℝ), weightReal c k) := by
  obtain ⟨P,hP,rfl⟩ := Option.map_eq_some_iff.mp h
  have hb' : boxes? (residual c) cb cd = some boxes := hb
  obtain ⟨ha,hab,hr⟩ := weight_sound hb' hP
  have hi := (cellIntegral_sound (SashimiEhDerivativeBounds.measurable_weight _) hab hr).1
  have he := square_weight_sound hab (hi.continuousOn_mul (continuous_id.pow 2).continuousOn) hr
  have ha' : (0:ℝ) < a := by exact_mod_cast ha
  have hab' : (a:ℝ) ≤ b := by exact_mod_cast hab
  have eqIntegral : (∫ k in (a:ℝ)..(b:ℝ), weightReal c k) =
      ∫ k in (a:ℝ)..(b:ℝ), k^2 * weightReal (residual c) k := by
    apply intervalIntegral.integral_congr
    intro k hk
    rw [uIcc_of_le hab'] at hk
    exact weight_eq (ha'.trans_le hk.1)
  rw [eqIntegral]
  exact he

#print axioms integral_sound
end SashimiEhMoment
