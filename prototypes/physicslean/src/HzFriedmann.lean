import HzCertificate
import Physlib.Cosmology.FLRW.Distances

/-!
# From a Friedmann solution to the certified Hubble function

The fluid supplies all the matter and obeys the dust density scaling law between
the two times. Radiation is the empty fluid. Flatness and the Friedmann closure
relation determine the cosmological constant's density parameter at the
reference time. The numerical value of H₀ and the resulting H use the same
inverse-time unit.
-/

namespace HzFriedmann

open Cosmology Cosmology.FLRW Cosmology.FLRW.FriedmannEquation

/-- The expanding solution at the emission time equals the exact input to the
interval certificates. Rewriting by this equality transports any certificate
for the corresponding H₀, Ωm and redshift to the solution's Hubble parameter. -/
theorem hubbleConstant_eq_modelH {U : FLRW} {F : CosmicFluid} {Λ G : ℝ}
    {c : SpeedOfLight} (sol : U.FriedmannSolution F Λ G c)
    (flat : U.geometry = .Flat) {t t₀ : Time} (ht : t ∈ U.I) (ht₀ : t₀ ∈ U.I)
    (hH₀ : 0 < hubbleConstant U.a t₀) (hH : 0 ≤ hubbleConstant U.a t)
    (dust : F.ρ t = F.ρ t₀ * (U.a t₀ / U.a t) ^ 3) :
    hubbleConstant U.a t =
      HzCertificate.modelH (hubbleConstant U.a t₀) (U.densityParameter F G t₀)
        (U.redshift t₀ t) := by
  have hK : U.curvatureDensityParameter c t₀ = 0 := by
    simp [FLRW.curvatureDensityParameter, flat]
  have hsum := sol.densityParameter_add_lambda_add_curvature ht₀ hH₀.ne'
  rw [hK, add_zero] at hsum
  have hΛ : U.lambdaDensityParameter Λ c t₀ = 1 - U.densityParameter F G t₀ := by
    exact eq_sub_of_add_eq' hsum
  have hRad : U.densityParameter CosmicFluid.empty G t₀ = 0 := by
    simp [FLRW.densityParameter, CosmicFluid.empty]
  have h := sol.hubbleConstant_eq_reducedHubble (Fm := F) (Fr := CosmicFluid.empty)
    (by simp [CosmicFluid.empty]) ht ht₀ hH₀ hH dust (by simp [CosmicFluid.empty])
  simpa only [HzCertificate.modelH, U.one_add_redshift, one_div, inv_div,
    hΛ, hRad, hK] using h

#print axioms hubbleConstant_eq_modelH

end HzFriedmann
