import SashimiMassVariance

/-! A kernel replay of the physical W sharp-k mass geometry on the two WMAP7
CAMB knots bracketing its moving cutoff. Power is deliberately the supplied
unmodified CDM field in this small regression; WDM and FDM runtime results
on full supplied files and their transfer parameters are in results/mass-variance.json. -/
open SashimiMassVariance SashimiTable FloatLib.Numerics FloatLib.Numerics.Interval
private def rows : List Knot := [⟨342032817/500000000, 6027787/34300⟩,
  ⟨349999783/500000000, 57120377/343000⟩]
private def output := evaluate? .sharp .linear .continuous (fun _ => .const 1) .sashimiW
  36700000000 100000000000000 1 rows 0.696 0.699 0.001 0 40 12
set_option maxRecDepth 200000 in
set_option maxHeartbeats 4000000 in
theorem checked_mass_cell :
    (output.any fun J => decide (0.006 < J.lo ∧ J.hi < 0.008)) = true := by decide +kernel
#print axioms checked_mass_cell
