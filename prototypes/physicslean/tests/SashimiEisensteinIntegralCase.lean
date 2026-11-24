import HzAudit
import SashimiEisensteinIntegral

/-! A two-cell impulse check of the continuous, σ₈-normalized EH98b variance on
a narrow band; the full mass-grid experiment uses separate coefficient and panel certificates. -/
open SashimiEisensteinIntegral FloatLib.Numerics SashimiSpectrum

set_option maxRecDepth 150000 in
set_option maxHeartbeats 7000000 in
theorem checked_local_EH_variance :
    (⟨0.4, 2⟩ : RationalInterval).Contains
    (normalize8 matched.sigma8 (raw matched 0.1 0.11 8)
      (raw matched 0.1 0.11 (radius matched 1e12))) :=
  variance_contains_of_check (a := 0.1) (ends := [0.11]) (depth := 1) (bits := 80) (degree := 18)
    (tol := 2) (by decide +kernel)

#print axioms checked_local_EH_variance

run_cmd PhysicsLeanAudit.currentWitnesses
