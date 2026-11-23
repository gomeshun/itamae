import SashimiMassVariance

/-! A kernel replay of the physical W sharp-k mass geometry on the two WMAP7
CAMB knots bracketing its moving cutoff. Power is deliberately the supplied
unmodified CDM field in this small regression; WDM and FDM runtime results
on full supplied files and their transfer parameters are in results/mass-variance.json. -/
open SashimiMassVariance SashimiTable FloatLib.Numerics FloatLib.Numerics.Interval
private def rows : List Knot := [⟨342032817/500000000, 6027787/34300⟩,
  ⟨349999783/500000000, 57120377/343000⟩]
private def output (c : Calibration) := evaluate? .sharp .linear .continuous (fun _ => .const 1) c
  36700000000 100000000000000 1 rows 0.696 0.699 0.001 0 40 12
set_option maxRecDepth 200000 in
set_option maxHeartbeats 4000000 in
theorem checked_mass_cell :
    ((output .wFamilyLegacy).any fun J => decide (0.006 < J.lo ∧ J.hi < 0.008)) = true := by decide +kernel
#print axioms checked_mass_cell

/- The public main cutoff lies below this same cell, so it contributes zero. -/
set_option maxRecDepth 200000 in
set_option maxHeartbeats 4000000 in
theorem checked_w_main_mass_cell :
    ((output .wMain).any fun J => decide
      (-1/10^10 < J.lo ∧ J.lo ≤ 0 ∧ 0 ≤ J.hi ∧ J.hi < 1/10^10)) = true := by decide +kernel
#print axioms checked_w_main_mass_cell
