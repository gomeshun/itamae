import SashimiEhMoment
import HzAudit

/-! Replay one finite EH98 integral with the fiducial parameters and compare
its enclosure width with that of the range-times-width rule on the same cell. -/
namespace SashimiEhMomentCase
open FloatLib.Numerics SashimiEisensteinIntegral SashimiEhSweep SashimiEhMoment

set_option maxRecDepth 150000

private def a : ℚ := 1/100
private def b : ℚ := 11/1000
private def coeff := (boxes? matched 80 16).getD []
private def moment := (integral? matched coeff a b 80 16).getD (.point 0)
private def rectangle :=
  ((weight? matched coeff a b 80 16).map (·.scaleNonnegative (b-a))).getD (.point 0)

theorem fiducial_cell : moment.Contains (∫ k in (a:ℝ)..(b:ℝ), weightReal matched k) := by
  exact integral_sound (c := matched) (cb := 80) (cd := 16) (boxes := coeff)
    (bits := 80) (degree := 16) (by decide +kernel) (by decide +kernel)

example : 0 < moment.lo ∧ moment.hi < (1/1000000000 : ℚ) := by decide +kernel
example : moment.hi - moment.lo < (2/3) * (rectangle.hi - rectangle.lo) := by decide +kernel
example : integral? matched coeff 0 1 80 16 = none := by decide +kernel
example : integral? matched coeff 1 (1/2) 80 16 = none := by decide +kernel

#print axioms fiducial_cell
end SashimiEhMomentCase

run_cmd PhysicsLeanAudit.currentWitnesses
