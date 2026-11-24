import SashimiEhDerivativeBounds

/-! Differentiation in physical halo mass at fixed positive wavenumber.
The exact factor from the Lagrangian radius is R/(3 M). -/
namespace SashimiEhPointDerivative
open SashimiSpectrum SashimiEhDerivativeBounds SashimiEhVariance SashimiEisensteinIntegral SashimiEhSweep

theorem radius_hasDerivAt {rho m : ℝ} (hr : 0 < rho) (hm : 0 < m) :
    HasDerivAt (massRadius rho) (massRadius rho m / (3*m)) m := by
  let x := 3*m/(4*Real.pi*rho)
  have hx : 0 < x := by dsimp [x]; positivity
  have hb := (((hasDerivAt_id m).const_mul 3).div_const (4*Real.pi*rho))
  have hp := (Real.hasDerivAt_rpow_const (p := (1/3:ℝ)) (Or.inl hx.ne')).comp m hb
  apply hp.congr_deriv
  rw [Real.rpow_sub_one hx.ne']
  unfold massRadius
  dsimp [x]
  field_simp

/-- P is the k-dependent weight and is held fixed as the mass varies. -/
theorem kernel_hasDerivAt (P : ℝ → ℝ) {rho m k : ℝ}
    (hr : 0 < rho) (hm : 0 < m) (hk : 0 < k) :
    HasDerivAt (fun t => P k * topHat (k * massRadius rho t)^2)
      (P k * shape (k * massRadius rho m) / m) m := by
  have hR : 0 < massRadius rho m := by unfold massRadius; positivity
  have pos := mul_pos hk hR
  have hd := (((topHat_hasDerivAt _ pos).comp m ((radius_hasDerivAt hr hm).const_mul k)).pow 2).const_mul (P k)
  apply hd.congr_deriv
  rw [shape_eq pos]
  simp only [Function.comp_apply]
  ring

theorem point_deriv {m k : ℝ} (hm : 0 < m) (hk : 0 < k) :
    HasDerivAt (fun t => weightReal matched k * topHat (k * massRadius massDensity t)^2)
      (weightReal matched k * shape (k * massRadius massDensity m) / m) m :=
  kernel_hasDerivAt _ (by norm_num [massDensity]) hm hk

/-- Positive k keeps the tilting logarithm and transfer function continuous. -/
theorem weight_continuous : ContinuousOn (weightReal matched) (Set.Ioi 0) := by
  have hc : coefficients matched = SashimiEisensteinHu.fiducial := by
    norm_num [coefficients, matched, SashimiEisensteinHu.fiducial]
  have b := SashimiEisensteinHu.fiducial_bounds
  have ht := SashimiEhContinuity.transfer_continuous (lt_of_lt_of_le (by norm_num) b.h.1)
    (lt_of_lt_of_le (by norm_num) b.keq.1) (lt_of_lt_of_le (by norm_num) b.s.1)
    (lt_of_lt_of_le (by norm_num) b.silk.1) (lt_of_lt_of_le (by norm_num) b.ac.1)
    (lt_of_lt_of_le (by norm_num) b.bc.1) b.bb.1 b.node.1
  unfold weightReal
  rw [hc]
  apply ContinuousOn.mul _ (ht.pow 2)
  intro k hk
  exact (Real.continuous_exp.continuousAt.comp
    (continuousAt_const.mul (Real.continuousAt_log (show 0 < k from hk).ne'))).continuousWithinAt

#print axioms point_deriv
end SashimiEhPointDerivative
