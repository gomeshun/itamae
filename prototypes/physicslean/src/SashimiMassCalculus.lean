import SashimiMassVariance
import SashimiSpectrumSlope

/-! Derivative of the finite, cellwise top-hat integral at the true physical mass. -/
namespace SashimiMassVariance
open SashimiSpectrum SashimiTable Real Set Filter
open scoped Topology

private theorem change_log {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (G : ℝ → ℝ) :
    (∫ t in log a..log b, G (exp t) * exp t) = ∫ k in a..b, G k := by
  simpa only [Function.comp_apply, exp_log ha, exp_log hb] using
    (intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg (a := log a) (b := log b)
      (g := G) continuous_exp.continuousOn (fun s _ => hasDerivAt_exp s) (fun s _ => (exp_pos s).le))

private theorem arg_eval {rho m : ℝ} (hr : 0 < rho) (hm : 0 < m) (s : ℝ) :
    logMassArgument (massRadius rho 1) (log m) s = exp s * massRadius rho m := by
  simpa only [div_one, exp_log hm] using logMassArgument_eq rho 1 (log m) s hr

/-- A single positive-k cell of a continuous real power interpolant. -/
theorem topHatCell_hasDerivAt (P : ℝ → ℝ) (hP : Continuous P)
    {rho m a b : ℝ} (hr : 0 < rho) (hm : 0 < m) (ha : 0 < a) (hb : 0 < b) (D : ℝ) :
    HasDerivAt (fun z => ∫ k in a..b, kernel .topHat P rho z D k)
      (∫ k in a..b, kernel .topHatDerivative P rho m D k) m := by
  have hB : 0 < massRadius rho 1 := by unfold massRadius; positivity
  have hd := (((smoothLogVariance_hasDerivAt P hP (massRadius rho 1) (log a) (log b) (log m) hB).comp m
    (hasDerivAt_log hm.ne')).const_mul (D^2))
  apply (hd.congr_of_eventuallyEq (f₁ := fun z => ∫ k in a..b, kernel .topHat P rho z D k) ?_).congr_deriv
  · change D^2*((∫ s in log a..log b, smoothLogPrime P (massRadius rho 1) (log m) s) * m⁻¹) = _
    rw [← change_log ha hb (kernel .topHatDerivative P rho m D)]
    rw [← intervalIntegral.integral_mul_const, ← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr
    intro s _
    simp only [smoothLogPrime, arg_eval hr hm, kernel, delta2]
    ring
  · filter_upwards [Ioi_mem_nhds hm] with z hz
    change (∫ k in a..b, kernel .topHat P rho z D k) = D^2 * ∫ s in log a..log b, smoothLogIntegrand P (massRadius rho 1) (log z) s
    rw [← change_log ha hb (kernel .topHat P rho z D), ← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_congr
    intro s _
    simp only [kernel, smoothLogIntegrand, delta2, arg_eval hr (show 0 < z from hz)]
    ring

/-- A continuous restriction on the actual cell suffices. Clamping k builds
a global continuous extension needed for the log-variable proof. -/
theorem topHatCell_of_continuousOn (P : ℝ → ℝ) {rho m a b : ℝ}
    (hr : 0 < rho) (hm : 0 < m) (ha : 0 < a) (hab : a ≤ b)
    (hP : ContinuousOn P (Icc a b)) (D : ℝ) :
    HasDerivAt (fun z => ∫ k in a..b, kernel .topHat P rho z D k)
      (∫ k in a..b, kernel .topHatDerivative P rho m D k) m := by
  let ext := fun k => P (max a (min b k))
  have he : Continuous ext := by
    apply hP.comp_continuous (continuous_const.max (continuous_const.min continuous_id))
    intro x
    exact ⟨le_max_left _ _, max_le hab (min_le_left _ _)⟩
  have hconst (q : Quantity) (z : ℝ) :
      (∫ k in a..b, kernel q ext rho z D k) = ∫ k in a..b, kernel q P rho z D k := by
    apply intervalIntegral.integral_congr
    intro k hk
    rw [uIcc_of_le hab] at hk
    have hx : ext k = P k := by simp [ext, min_eq_right hk.2, max_eq_right hk.1]
    cases q <;> simp only [kernel, hx]
  simpa only [hconst] using topHatCell_hasDerivAt ext he hr hm ha (ha.trans_le hab) D

/-- Differentiate all accepted adjacent cells at the same positive mass.
Continuity is required only on the portion of each input cell used. -/
theorem topHatBand_hasDerivAt (P : Cell → ℝ → ℝ) {rho m : ℝ} (a b : ℚ)
    (hr : 0 < rho) (hm : 0 < m) (ha : 0 < a) (D K : ℝ) (cs : List Cell)
    (hP : ∀ c ∈ cs, ContinuousOn (P c) (Icc ((max a c.left.k : ℚ):ℝ) ((min b c.right.k : ℚ):ℝ))) :
    HasDerivAt (fun z => massBand .topHat (fun c => kernel .topHat (P c) rho z D) K a b cs)
      (massBand .topHatDerivative (fun c => kernel .topHatDerivative (P c) rho m D) K a b cs) m := by
  induction cs with
  | nil => simpa [massBand] using (hasDerivAt_const m (0:ℝ))
  | cons c cs ih =>
    have hcs := ih (fun c hc => hP c (by simp [hc]))
    dsimp only [massBand]
    split
    · rename_i hlr
      have hl : 0 < ((max a c.left.k : ℚ):ℝ) := by exact_mod_cast ha.trans_le (le_max_left _ _)
      have hlr' : ((max a c.left.k : ℚ):ℝ) ≤ ((min b c.right.k : ℚ):ℝ) := by exact_mod_cast hlr.le
      exact (topHatCell_of_continuousOn (P c) hr hm hl hlr' (hP c (by simp)) D).add hcs
    · simpa only [zero_add] using hcs

/-- The exact physical cutoff moves like M^(-1/3), independently of the
choice of positive calibration. -/
theorem cutoff_hasDerivAt {rho c m : ℝ} (hr : 0 < rho) (hc : 0 < c) (hm : 0 < m) :
    HasDerivAt (fun z => c / massRadius rho z) (-(c/massRadius rho m)/(3*m)) m := by
  have hd : HasDerivAt (logCutoff rho c) (-(1/m)/3) m := by
    unfold logCutoff
    convert (((hasDerivAt_log hm.ne').const_add (log (3/(4*π*rho)))).div_const 3).const_sub (log c) using 1
    ring
  apply (hd.exp.congr_of_eventuallyEq (f₁ := fun z => c / massRadius rho z) ?_).congr_deriv
  · rw [exp_logCutoff _ _ _ hr hc hm]; ring
  · filter_upwards [Ioi_mem_nhds hm] with z hz
    exact (exp_logCutoff rho c z hr hc hz).symm

noncomputable def sharpCellSlope (P : ℝ → ℝ) (m D kc l r : ℝ) :=
  if l < kc ∧ kc < r then -(D^2)*kc^3*P kc / (6*m*π^2) else 0

/-- Variation of each cell. Cutoffs strictly outside that cell contribute
zero derivative; a strictly interior cutoff contributes its boundary power. -/
theorem sharpCell_hasDerivAt (P : ℝ → ℝ) (hP : Continuous P) {rho c m l r : ℝ} (D : ℝ)
    (hr : 0 < rho) (hc : 0 < c) (hm : 0 < m) (hlr : l < r)
    (hneL : c/massRadius rho m ≠ l) (hneR : c/massRadius rho m ≠ r) :
    HasDerivAt (fun z => ∫ k in l..upper .sharp (c/massRadius rho z) l r, kernel .sharp P rho z D k)
      (sharpCellSlope P m D (c/massRadius rho m) l r) m := by
  let K := fun z => c/massRadius rho z
  let f := fun k => D^2 * k^2 * P k / (2*π^2)
  have hf : Continuous f := by dsimp [f]; fun_prop
  have hk := cutoff_hasDerivAt hr hc hm
  change HasDerivAt (fun z => ∫ k in l..max l (min r (K z)), f k) _ m
  rcases lt_or_gt_of_ne hneL with below | above
  · have he : (fun z => ∫ k in l..max l (min r (K z)), f k) =ᶠ[𝓝 m] fun _ => 0 := by
      filter_upwards [hk.continuousAt.eventually_lt_const below] with z hz
      simp [max_eq_left ((min_le_right r (K z)).trans hz.le)]
    have hh : ¬ (l < K m ∧ K m < r) := fun h => (below.not_gt h.1)
    simpa only [sharpCellSlope, K, ite_eq_right hh] using (hasDerivAt_const m (0:ℝ)).congr_of_eventuallyEq he
  · rcases lt_or_gt_of_ne hneR with inside | aboveR
    · have he : (fun z => ∫ k in l..max l (min r (K z)), f k) =ᶠ[𝓝 m] fun z => ∫ k in l..K z, f k := by
        filter_upwards [hk.continuousAt.eventually_lt_const inside, hk.continuousAt.eventually_const_lt above] with z h1 h2
        rw [min_eq_right h1.le,max_eq_right h2.le]
      have hh : l < K m ∧ K m < r := ⟨above,inside⟩
      apply (((hf.integral_hasStrictDerivAt l (K m)).hasDerivAt.comp m hk).congr_of_eventuallyEq he).congr_deriv
      rw [sharpCellSlope, ite_eq_left (And.intro above inside)]
      unfold f K; ring
    · have he : (fun z => ∫ k in l..max l (min r (K z)), f k) =ᶠ[𝓝 m] fun _ => ∫ k in l..r, f k := by
        filter_upwards [hk.continuousAt.eventually_const_lt aboveR] with z hz
        rw [min_eq_left hz.le,max_eq_right hlr.le]
      have hh : ¬ (l < K m ∧ K m < r) := fun h => aboveR.not_gt h.2
      simpa only [sharpCellSlope, K, ite_eq_right hh] using (hasDerivAt_const m (∫ k in l..r, f k)).congr_of_eventuallyEq he

/-- Sum the actual boundary contributions across the cells. A cutoff at a
tabulated kink is excluded. The derivative computation checks a strict
bracket for the active cell; inactive cells have zero boundary contribution. -/
noncomputable def sharpBandSlope (P : Cell → ℝ → ℝ) (m D kc : ℝ) (a b : ℚ) : List Cell → ℝ
  | [] => 0
  | c::cs =>
    let l := max a c.left.k
    let r := min b c.right.k
    (if l < r then sharpCellSlope (P c) m D kc l r else 0) + sharpBandSlope P m D kc a b cs

theorem sharpBand_hasDerivAt (P : Cell → ℝ → ℝ) {rho c m : ℝ} (a b : ℚ) (D : ℝ)
    (hr : 0 < rho) (hc : 0 < c) (hm : 0 < m) (cs : List Cell)
    (hP : ∀ C ∈ cs, ContinuousOn (P C) (Icc ((max a C.left.k : ℚ):ℝ) ((min b C.right.k : ℚ):ℝ)))
    (hfree : ∀ C ∈ cs, max a C.left.k < min b C.right.k →
      c/massRadius rho m ≠ ((max a C.left.k : ℚ):ℝ) ∧ c/massRadius rho m ≠ ((min b C.right.k : ℚ):ℝ)) :
    HasDerivAt (fun z => massBand .sharp (fun C => kernel .sharp (P C) rho z D) (c/massRadius rho z) a b cs)
      (sharpBandSlope P m D (c/massRadius rho m) a b cs) m := by
  induction cs with
  | nil => simpa [massBand, sharpBandSlope] using hasDerivAt_const m (0:ℝ)
  | cons C cs ih =>
    have hcs := ih (fun j hj => hP j (by simp [hj])) (fun j hj => hfree j (by simp [hj]))
    dsimp only [massBand, sharpBandSlope]
    split
    · rename_i hlr
      let l : ℝ := ((max a C.left.k : ℚ):ℝ)
      let r : ℝ := ((min b C.right.k : ℚ):ℝ)
      have hlr' : l < r := by dsimp [l,r]; exact_mod_cast hlr
      let ext := fun k => P C (max l (min r k))
      have he : Continuous ext := by
        apply (hP C (by simp)).comp_continuous (continuous_const.max (continuous_const.min continuous_id))
        intro x
        exact ⟨le_max_left _ _, max_le hlr'.le (min_le_left _ _)⟩
      have heq (z : ℝ) : (∫ k in l..upper .sharp (c/massRadius rho z) l r, kernel .sharp ext rho z D k)
          = ∫ k in l..upper .sharp (c/massRadius rho z) l r, kernel .sharp (P C) rho z D k := by
        apply intervalIntegral.integral_congr
        intro k hk
        change k ∈ uIcc l (max l (min r _)) at hk
        rw [uIcc_of_le (le_max_left _ _)] at hk
        have hkr := hk.2.trans (max_le hlr'.le (min_le_left _ _))
        simp only [kernel,ext,min_eq_right hkr,max_eq_right hk.1]
      have hsl : sharpCellSlope ext m D (c/massRadius rho m) l r = sharpCellSlope (P C) m D (c/massRadius rho m) l r := by
        unfold sharpCellSlope; split
        · rename_i hin
          simp only [ext, min_eq_right hin.2.le, max_eq_right hin.1.le]
        · rfl
      obtain ⟨h1,h2⟩ := hfree C (by simp) hlr
      have h := sharpCell_hasDerivAt ext he D hr hc hm hlr' h1 h2
      simp only [heq,hsl] at h
      exact h.add hcs
    · simpa only [zero_add] using hcs

#print axioms sharpBand_hasDerivAt
#print axioms sharpCell_hasDerivAt
#print axioms topHatBand_hasDerivAt
end SashimiMassVariance
