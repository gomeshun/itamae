import HzAudit

/-! Import via the audited public root. These two downstream claims use only
the matched physical inputs, so extra continuity or transfer-envelope premises
in the production API would break this check. The remaining production
declarations are covered by HzAudit's exhaustive transitive scan. -/
namespace SashimiEisensteinAudit
open SashimiEisensteinHu SashimiEhVariance

theorem derivative {M a b : ℝ} (hm : 0 < M) (ha : 0 < a) (hab : a ≤ b) :
    HasDerivAt (fun m => normalizedS m a b) (normalizedDerivative M a b) M :=
  SashimiEhContinuity.matched_derivative hm ha hab

theorem ultraviolet {k : ℝ} (hk : 10^10 ≤ k) :
    |transferReal fiducial k| ≤ (1/10^9)/k :=
  high_envelopeBig fiducial_bounds hk

#check SashimiEhFullGrid.check_sound
#check SashimiEhFullCalculus.hasDerivAt_fullRaw
#check SashimiEhFullDerivative.hasDerivAt_fullS
#check SashimiEhYang.density_sound
#check SashimiEhDerivativeBounds.fromMesh_sound
#check SashimiHeath.check_sound

run_cmd PhysicsLeanAudit.currentWitnesses
end SashimiEisensteinAudit
