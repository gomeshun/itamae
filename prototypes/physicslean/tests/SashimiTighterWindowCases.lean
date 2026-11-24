import HzAudit
import SashimiTighterWindow

/-! Kernel-recomputed improvements on fixed dimensionless wavenumber cells.
    Both comparisons use the same endpoints, precision, and polynomial degree. -/
namespace SashimiTighterWindowCases

open FloatLib.Numerics SashimiTighterWindow

private def width (I : RationalInterval) : ℚ := I.hi - I.lo

set_option maxRecDepth 100000
set_option maxHeartbeats 16000000

theorem taylor_cell :
    10 * width (windowSq 1 (26/25) 112 16) ≤ width (SashimiTableWindow.windowSq 1 (26/25) 112 16) ∧
    10 * width (weight 1 (26/25) 112 16) ≤ width (SashimiEhDerivativeBounds.weight 1 (26/25) 112 16) := by
  decide +kernel

theorem oscillatory_cell :
    3 * width (windowSq (14/5) (293/100) 112 16) ≤ width (SashimiTableWindow.windowSq (14/5) (293/100) 112 16) ∧
    4 * width (weight (14/5) (293/100) 112 16) ≤ width (SashimiEhDerivativeBounds.weight (14/5) (293/100) 112 16) := by
  decide +kernel

#print axioms taylor_cell
#print axioms oscillatory_cell

end SashimiTighterWindowCases

run_cmd PhysicsLeanAudit.currentWitnesses
