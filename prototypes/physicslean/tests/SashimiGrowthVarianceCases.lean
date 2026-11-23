import SashimiGrowthVariance

/-! z=2 regression on data rows 1201–1203 of the family SASHIMI-F
Planck2018_CAMB_extrap.dat, SHA-256 91269da31616bf5f8d506ad138f157c8eceb284d4cd328cd99e730b882a8dd61.
These are the cells intersecting [0.69,0.706] Mpc⁻¹. -/
namespace SashimiGrowthVarianceCases
open FloatLib.Numerics FloatLib.Numerics.Interval SashimiTable
open SashimiGrowthVariance

private def rows : List Knot :=
  [⟨1.00926,66.0392⟩, ⟨1.03277,62.6021⟩, ⟨1.05684,59.3382⟩].map (fromHUnits 0.674)

private def raw (q : SashimiMassVariance.Quantity) (tol : ℚ) :=
  SashimiMassVariance.evaluate? q .logLog .atKnots (SashimiSpectrum.fdmPowerExpr (.const 1))
    (.fixed 1) 39700000000 500000000000 1 rows 0.69 0.706 tol 2 60 12

private def integral : RationalInterval := ⟨144405740023806901/2305843009213693952,
  313389492235947019/4611686018427387904⟩
private def slope : RationalInterval := ⟨-134193/4611686018427387904, -9055/1152921504606846976⟩

set_option maxRecDepth 16000 in
theorem accepted_raw : raw .topHat (2/100) = some integral ∧ raw .topHatDerivative (1/10^13) = some slope := by
  decide +kernel

private def enclosed (q : Quantity) (r : RationalInterval) (tol lo hi : ℚ) : Bool :=
  (atRedshift? q 0.315 0.685 2 r tol (1/10000) 60 12).any (fun i => decide (lo ≤ i.lo ∧ i.hi ≤ hi))

theorem cpt_z2 :
    enclosed .variance integral 0.01 0.0109226 0.0118522 ∧
    enclosed .sigma integral 0.01 0.1045113 0.1088677 ∧
    enclosed .derivative slope (1/10^13) (-5.076e-15) (-1.369e-15) ∧
    enclosed .legacyFDerivative slope (1/10^13) (-1.216e-14) (-3.279e-15) := by
  decide +kernel

example : atRedshift? .variance 0.315 0.685 2 integral 0 (1/10000) 60 12 = none := by decide +kernel
example : combine? .sigma (.point 1) slope 1 60 = none := by decide +kernel

#print axioms accepted_raw
#print axioms cpt_z2
end SashimiGrowthVarianceCases
