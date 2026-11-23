import SashimiSpectrumSlopeCore

namespace SashimiSpectrum
open Real Set Filter MeasureTheory
open scoped Topology
noncomputable section

private theorem logMassArgument_ne {B : ℝ} (hB : 0 < B) (p : ℝ × ℝ) :
    logMassArgument B p.1 p.2 ≠ 0 := (mul_pos (mul_pos (exp_pos _) hB) (exp_pos _)).ne'

private theorem logMassArgument_continuous (B : ℝ) :
    Continuous (fun p : ℝ × ℝ => logMassArgument B p.1 p.2) :=
  ((Real.continuous_exp.comp continuous_snd).mul_const B).mul
    (Real.continuous_exp.comp (continuous_fst.div_const 3))

private theorem smoothLog_continuous (P : ℝ → ℝ) (hP : Continuous P) (B : ℝ) (hB : 0 < B) :
    Continuous (Function.uncurry (smoothLogIntegrand P B)) ∧
    Continuous (Function.uncurry (smoothLogPrime P B)) := by
  let x := fun p : ℝ × ℝ => logMassArgument B p.1 p.2
  have hx : Continuous x := logMassArgument_continuous B
  have hn (p : ℝ × ℝ) : x p ≠ 0 := logMassArgument_ne hB p
  have hd : Continuous (fun p : ℝ × ℝ => delta2 P (exp p.2)) := by
    unfold delta2
    fun_prop
  have hc := Real.continuous_cos.comp hx
  have hs := Real.continuous_sin.comp hx
  have hw : Continuous (fun p => topHat (x p)) := by
    simp only [topHat, hn, ↓reduceIte]
    apply Continuous.div _ (hx.pow 3) (fun p => pow_ne_zero 3 (hn p))
    exact (hs.sub (hx.mul hc)).const_mul 3
  have hw' : Continuous (fun p => topHatPrime (x p)) := by
    have h₁ := (hc.const_mul (9:ℝ)).div (hx.pow 3) (fun p => pow_ne_zero 3 (hn p))
    have h₂ := (hs.const_mul (9:ℝ)).div (hx.pow 4) (fun p => pow_ne_zero 4 (hn p))
    have h₃ := (hs.const_mul (3:ℝ)).div (hx.pow 2) (fun p => pow_ne_zero 2 (hn p))
    exact (h₁.sub h₂).add h₃
  exact ⟨hd.mul (hw.pow 2), (((((hd.mul_const 2).mul hw).mul hw').mul hx).div_const 3)⟩

private theorem differentiate_band (f d : ℝ → ℝ → ℝ)
    (hf : Continuous (Function.uncurry f)) (hd : Continuous (Function.uncurry d))
    (hder : ∀ u s, HasDerivAt (fun v => f v s) (d u s) u) (a b u : ℝ) :
    HasDerivAt (fun v => ∫ s in a..b, f v s) (∫ s in a..b, d u s) u := by
  let around : Set ℝ := Icc (u-1) (u+1)
  have hnear : around ∈ 𝓝 u := Icc_mem_nhds (by linarith) (by linarith)
  have hcomp : IsCompact (around ×ˢ uIcc a b) := isCompact_Icc.prod isCompact_uIcc
  obtain ⟨C, hC⟩ := hcomp.exists_bound_of_continuousOn hd.continuousOn
  have section_cont {g : ℝ → ℝ → ℝ} (h : Continuous (Function.uncurry g)) (v : ℝ) :
      Continuous (g v) := h.comp (continuous_const.prodMk continuous_id)
  have cd := section_cont hd u
  have cf := section_cont hf
  refine (intervalIntegral.hasDerivAt_integral_of_dominated_loc_of_deriv_le
    (F := f) (F' := d) (μ := volume) (a := a) (b := b) (x₀ := u) (s := around)
    (bound := fun _ => C) hnear ?_ ((cf u).intervalIntegrable a b)
    cd.aestronglyMeasurable ?_ intervalIntegrable_const ?_).2
  · exact Filter.Eventually.of_forall fun v => (cf v).aestronglyMeasurable
  · exact Filter.Eventually.of_forall fun t ht v hv => hC (v,t) ⟨hv, uIoc_subset_uIcc ht⟩
  · exact Filter.Eventually.of_forall fun t _ v _ => hder v t

/-- Differentiate the smoothed variance under the integral on a finite table
band. Its continuous derivative has a finite bound on a compact log-mass neighborhood. -/
theorem smoothLogVariance_hasDerivAt (P : ℝ → ℝ) (hP : Continuous P)
    (B a b u : ℝ) (hB : 0 < B) :
    HasDerivAt (fun v => ∫ s in a..b, smoothLogIntegrand P B v s)
      (∫ s in a..b, smoothLogPrime P B u s) u := by
  obtain ⟨hf, hd⟩ := smoothLog_continuous P hP B hB
  exact differentiate_band _ _ hf hd (fun v s => smoothLogIntegrand_hasDerivAt P B v s hB) a b u

end
end SashimiSpectrum
