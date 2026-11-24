import CertifiedIntegralRatio

/-! The same unknown positive weight appears in both integrals. A single
unit-width cell identifies their ratio exactly even if the weight spans [1,4]. -/
namespace CertifiedIntegralRatioCase
open CertifiedIntegralRatio FloatLib.Numerics MeasureTheory Set

private def input : Cell := ⟨⟨1,4⟩, .point 1, .point 1⟩

example : contribution 1 input = .point 0 := by decide +kernel
example : ¬ 0 ≤ (contribution (3/2) input).lo := by decide +kernel
example : ((⟨1,4⟩ : RationalInterval).divPositive ⟨1,4⟩) = ⟨1/4,4⟩ := by decide +kernel

theorem ratioOfSharedWeight (w : ℝ → ℝ) (hm : Measurable w)
    (hb : ∀ k ∈ Icc (1:ℝ) 2, (⟨1,4⟩:RationalInterval).Contains (w k)) :
    (RationalInterval.point 1).Contains ((∫ k in (1:ℝ)..2, w k)/(∫ k in (1:ℝ)..2, w k)) := by
  have hd := (SashimiEhSweep.cellIntegral_sound (l:=1) (u:=2) hm (by decide) hb).2.1
  have hd' : 0 < ∫ k in (1:ℝ)..2, w k := lt_of_lt_of_le (by norm_num [RationalInterval.scaleNonnegative]) hd
  have hs : (RationalInterval.point 0).Contains ((∫ k in (1:ℝ)..2, w k) -
      (1:ℚ) * (∫ k in (1:ℝ)..2, w k)) := by simp [RationalInterval.Contains,RationalInterval.point]
  exact ratio_bounds (tn:=0) (td:=0) hd' (by decide) (by decide)
    ⟨le_rfl,by simp⟩ ⟨le_rfl,by simp⟩ hs hs (by decide +kernel) (by decide +kernel)

#print axioms ratioOfSharedWeight
end CertifiedIntegralRatioCase
