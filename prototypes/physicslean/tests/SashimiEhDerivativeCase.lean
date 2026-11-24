import HzAudit
import SashimiEhDerivativeBounds

/-! Replay the entire single-cell transfer, denominator, radius and signed
numerator for one physical halo mass on [1,1.001] h/Mpc. -/
namespace SashimiEhDerivativeCase
open FloatLib.Numerics SashimiEhDerivativeBounds SashimiEhVariance

private def out : RationalInterval :=
  ⟨-351281551169/508000000000000000000, -617749321667/1111000000000000000000⟩

set_option maxRecDepth 100000 in
set_option maxHeartbeats 16000000 in
theorem replay : fromMesh? (10^12) 1 1 [1001/1000] 48 8 = some out := by decide +kernel

theorem example_contains : out.Contains (normalizedDerivative (10^12) 1 (1001/1000)) := by
  simpa using (fromMesh_sound replay).1

theorem example_negative : normalizedDerivative (10^12) 1 (1001/1000) < 0 :=
  example_contains.2.trans_lt (by norm_num [out])

theorem actual_derivative_negative :
    deriv (fun M => normalizedS M 1 (1001/1000)) (10^12) < 0 := by
  have hc := checked_deriv (by norm_num) (by norm_num) (by decide +kernel) replay
  have h := hc.2.trans_lt (show (out.hi : ℝ) < 0 by norm_num [out])
  simpa using h

#print axioms actual_derivative_negative
end SashimiEhDerivativeCase

run_cmd PhysicsLeanAudit.currentWitnesses
