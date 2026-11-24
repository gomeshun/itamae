import SashimiEhPointDerivative
import SashimiEhSlopeTails

/-! Differentiation of the improper EH variance integral. For masses near m
the radius stays above rmin and the mass above m/2; this provides a single
integrable bound for the derivative over the entire positive k axis. -/
namespace SashimiEhFullCalculus
open SashimiEhVariance SashimiEisensteinIntegral SashimiEhFull SashimiEhFullGrid
open SashimiSpectrum SashimiEhSweep SashimiEhDerivativeBounds SashimiEhPointDerivative
open SashimiEhSlopeTails
open MeasureTheory Set Filter
open scoped Topology

theorem hasDerivAt_fullRaw {m rmin : ℝ} (hm : 0 < m) (hr : 0 < rmin)
    (hl : rmin < massRadius massDensity m) (hKR : 1 ≤ (highEdge:ℝ)*rmin) :
    HasDerivAt (fun x => fullRaw matched (massRadius massDensity x))
      (fullSlope (massRadius massDensity m) / m) m := by
  let K : ℝ := highEdge
  let A : ℝ := 84*(highConstant:ℝ)^2/rmin^3
  let small := (Ioc (0:ℝ) K).indicator (fun _ => 36*K^3)
  let large := (Ioi K).indicator (fun k : ℝ => A*k^(-2:ℝ))
  have hK : 0 < K := by norm_num [K, highEdge]
  have hsmall : Integrable small volume :=
    (integrableOn_const (s := Ioc (0:ℝ) K) (C := 36*K^3) (measure_Ioc_lt_top.ne)).integrable_indicator measurableSet_Ioc
  have hlarge : Integrable large volume :=
    (integrable_indicator_iff measurableSet_Ioi).mpr
      ((integrableOn_Ioi_rpow_of_lt (by norm_num : (-2:ℝ) < -1) hK).const_mul A)
  let bound := fun k => (small k + large k) / (m/2)
  have hbound : IntegrableOn bound (Ioi (0:ℝ)) := ((hsmall.add hlarge).div_const (m/2)).integrableOn
  let radiusAt : ℝ → ℝ := massRadius massDensity
  have hnear : {x | m/2 < x ∧ rmin < radiusAt x} ∈ 𝓝 m := by
    have hc : ContinuousAt radiusAt m := (radius_hasDerivAt (by norm_num [massDensity]) hm).continuousAt
    exact Filter.inter_mem (Ioi_mem_nhds (by linarith : m/2 < m)) (hc (Ioi_mem_nhds hl))
  let F := fun x => spectrumWeight matched (radiusAt x)
  let G := fun x k => slope (radiusAt x) k / x
  have h_deriv := hasDerivAt_integral_of_dominated_loc_of_deriv_le (F := F) (F' := G)
    (μ := volume.restrict (Ioi (0:ℝ))) (bound := bound) hnear
    (Eventually.of_forall fun x => (measurableWeight _ _).aestronglyMeasurable) (fullRaw_integrable hr hl.le hKR)
    (((measurableSlope _).div_const m).aestronglyMeasurable) (by
      filter_upwards [ae_restrict_mem (μ := volume) (measurableSet_Ioi (a := (0:ℝ)))] with k hk
      intro x hx
      obtain ⟨hmx, hrx⟩ := hx
      have hslope : |slope (radiusAt x) k| ≤ small k + large k := by
        by_cases h : k ≤ K
        · have hk' : k ∈ Ioc 0 K := ⟨hk,h⟩
          have hnot : k ∉ Ioi K := not_lt.mpr h
          simpa [small, large, hk', hnot] using midBound (hr.trans hrx) hk h (by norm_num [K, highEdge])
        · have hfar : k ∈ Ioi K := lt_of_not_ge h
          have hnot : k ∉ Ioc 0 K := fun p => h p.2
          simpa [small, large, hfar, hnot, A] using highBound hr hrx.le le_rfl hKR hfar
      have hm2 : 0 < m/2 := by linarith
      have hxpos : 0 < x := hm2.trans hmx
      rw [Real.norm_eq_abs, show G x k = slope (radiusAt x) k / x from rfl, abs_div, abs_of_pos hxpos]
      exact div_le_div₀ ((abs_nonneg _).trans hslope) hslope hm2 hmx.le)
    hbound (by
      filter_upwards [ae_restrict_mem (μ := volume) (measurableSet_Ioi (a := (0:ℝ)))] with k hk
      intro x hx
      have hdx := point_deriv (show 0 < x by have hh := hx.1; linarith) hk
      apply hdx.congr_of_eventuallyEq
      filter_upwards with y
      change spectrumWeight matched (radiusAt y) k = _
      exact SashimiEhFull.weight_eq matched (radiusAt y) hk)
  simpa only [F, G, fullRaw, fullSlope, radiusAt, integral_div] using h_deriv.2

#print axioms hasDerivAt_fullRaw

end SashimiEhFullCalculus
