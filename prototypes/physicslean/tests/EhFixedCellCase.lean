import SashimiEhFixedReplay

/-! One actual 128-per-octave EH mesh cell, independently checked from the
cosmological parameters through the centered integral and its width. -/
namespace EhFixedCellCase
open FloatLib.Numerics SashimiEisensteinIntegral SashimiEhFixedCells

private abbrev left : ℚ := 262144/390625
private abbrev right : ℚ := 264192/390625
private def reported : RationalInterval := ⟨13343/10^11, 13356/10^11⟩

private def check? : Option Unit := do
  let params ← boxes? matched 96 16
  let q ← buildCell? matched params left right 96 16
  guard (reported.lo ≤ q.integral.lo ∧ q.integral.hi ≤ reported.hi ∧
    10*(q.integral.hi - q.integral.lo) < (q.pointwise.hi - q.pointwise.lo)*(right-left))

set_option maxRecDepth 100000 in
set_option maxHeartbeats 24000000 in
private theorem checked : check? = some () := by decide +kernel

/-- Both the physical integral enclosure and a strict factor of ten gain
over range-times-width follow from the very same fixed-cell computation. -/
theorem centered_cell :
    reported.Contains (∫ k in (left:ℝ)..(right:ℝ), SashimiEhSweep.weightReal matched k) ∧
    ∃ (params : List RationalInterval) (q : Cell),
      boxes? matched 96 16 = some params ∧ buildCell? matched params left right 96 16 = some q ∧
      10*(q.integral.hi - q.integral.lo) < (q.pointwise.hi - q.pointwise.lo)*(right-left) := by
  have h := checked
  obtain ⟨params,hp,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨q,hq,hg⟩ := Option.bind_eq_some_iff.mp h
  have hbounds : reported.lo ≤ q.integral.lo ∧ q.integral.hi ≤ reported.hi ∧
      10*(q.integral.hi - q.integral.lo) < (q.pointwise.hi - q.pointwise.lo)*(right-left) := by
    by_contra hn; simp [guard,hn] at hg
  obtain ⟨hedge,hv⟩ := buildCell_sound hp hq
  have hmass := hv.mass
  rw [hedge] at hmass
  exact ⟨⟨(by exact_mod_cast hbounds.1 : (reported.lo:ℝ) ≤ q.integral.lo).trans hmass.1,
    hmass.2.trans (by exact_mod_cast hbounds.2.1)⟩,params,q,hp,hq,hbounds.2.2⟩

/-- Exercise the public evaluator through both window integrals, the shared anchor
and the derivative with respect to physical mass on the same single-cell band. -/
private def bounded : Option Unit := do
  let v ← evaluate? (10^12) 100 left [right] 112 16
  guard (88 ≤ v.variance.lo ∧ v.variance.hi ≤ 93 ∧
    -9/10^12 ≤ v.derivative.lo ∧ v.derivative.hi ≤ -8/10^12)

set_option maxRecDepth 100000 in
set_option maxHeartbeats 40000000 in
private theorem normalized_checked : bounded = some () := by decide +kernel

theorem normalized_cell :
    (⟨88,93⟩ : RationalInterval).Contains
      (SashimiEhVariance.normalizedS ((10^12 : ℚ) : ℝ) left right) ∧
    (⟨-9/10^12,-8/10^12⟩ : RationalInterval).Contains
      (deriv (fun t => SashimiEhVariance.normalizedS t left right) ((10^12 : ℚ) : ℝ)) := by
  have h := normalized_checked
  obtain ⟨v,hv,h⟩ := Option.bind_eq_some_iff.mp h
  have hb : 88 ≤ v.variance.lo ∧ v.variance.hi ≤ 93 ∧
      -9/10^12 ≤ v.derivative.lo ∧ v.derivative.hi ≤ -8/10^12 := by
    by_contra hn; simp [guard,hn] at h
  obtain ⟨hV,hD,_⟩ := evaluate_sound (a := left) (ends := [right]) (by decide +kernel) (by decide +kernel) hv
  change v.variance.Contains (SashimiEhVariance.normalizedS _ left right) at hV
  change v.derivative.Contains (deriv (fun t => SashimiEhVariance.normalizedS t left right) _) at hD
  exact ⟨⟨(by exact_mod_cast hb.1 : (88:ℝ) ≤ v.variance.lo).trans hV.1,
      hV.2.trans (by exact_mod_cast hb.2.1)⟩,
    ⟨(by exact_mod_cast hb.2.2.1 : (((-9/10^12:ℚ):ℝ)) ≤ v.derivative.lo).trans hD.1,
      hD.2.trans (by exact_mod_cast hb.2.2.2)⟩⟩

set_option maxRecDepth 100000 in
set_option maxHeartbeats 60000000 in
/-- A falsely reported zero integral fails the fresh evaluator on four
consecutive cells of the 112-bit production grid. -/
theorem false_block_rejected :
    SashimiEhFixedReplay.block? (.point 8) 112 16 (32768/390625)
      [33024/390625,33280/390625,33536/390625,33792/390625] ≠ some (.point 0) := by
  decide +kernel

#print axioms false_block_rejected
#print axioms normalized_cell
#print axioms centered_cell
end EhFixedCellCase
