import SashimiGrowth
import SashimiBackground

/-! Kernel reductions of the executable enclosures. Bounds below are exact
rationals in the public units of the native ITAMAE backend. -/
namespace SashimiBackgroundCases
open FloatLib.Numerics

def enclosedIn (o : Option RationalInterval) (a b : ℚ) : Bool :=
  match o with
  | some I => decide (a ≤ I.lo ∧ I.hi ≤ b)
  | none => false

set_option maxHeartbeats 1000000 in
theorem growth_sample : enclosedIn
    (SashimiGrowth.enclose .growth (315/1000) (685/1000) (674/1000) 2 (1/10000) 32 12)
    (4175/10000) (4178/10000) = true := by decide +kernel

set_option maxHeartbeats 1000000 in
theorem density_sample : enclosedIn
    (SashimiBackground.enclose .virialOverdensity (SashimiBackground.flat (315/1000)) (674/10)
      (430091/10^14) 2 (1/1000) 32 12) (171323/1000) (171325/1000) = true := by decide +kernel

theorem analytic_rate_z2 : enclosedIn
    (SashimiGrowth.enclose .growthDerivative (315/1000) (685/1000) (674/1000) 2 (1/10000) 32 12)
    (-1336/10000) (-1334/10000) = true := by decide +kernel

theorem legacy_f_rate_z2 : enclosedIn
    (SashimiGrowth.enclose .legacyFDerivative (315/1000) (685/1000) (674/1000) 2 (1/10000) 32 12)
    (-1349/10000) (-1347/10000) = true := by decide +kernel

theorem cpt_threshold_z2 : enclosedIn
    (SashimiGrowth.enclose .collapseThreshold (315/1000) (685/1000) (674/1000) 2 (1/10000) 32 12)
    (40365/10000) (40378/10000) = true := by decide +kernel

theorem wmap7_g_reference :
    SashimiBackground.e2 (SashimiBackground.wLegacy (27/100) (469/10000) 0) 0 = 10469/10000 := by
  rw [SashimiBackground.wLegacy_g_zero]
  norm_num

theorem rejects_zero_matter :
    SashimiGrowth.enclose .growth 0 1 1 0 (1/1000) = none := by decide +kernel

theorem rejects_tight_budget :
    SashimiBackground.enclose .hubble (SashimiBackground.flat (3/10)) 70 1 1 0 12 12 = none := by
  decide +kernel

#print axioms growth_sample
#print axioms density_sample

end SashimiBackgroundCases
