import SashimiEhSweep
import SashimiEisensteinTail

/-! Positive-k continuity of the baryonic EH transfer on the cosmology used
in the matched calculation, including its oscillatory baryon term. -/
namespace SashimiEhContinuity
open SashimiEisensteinHu SashimiEhVariance Real Set

private theorem contPow {f : ℝ → ℝ} {x : ℝ} (hf : ContinuousAt f x) (h : 0 < f x) (p : ℝ) :
    ContinuousAt (fun t => powPos (f t) p) x :=
  Real.continuous_exp.continuousAt.comp ((hf.log h.ne').mul_const p)

private theorem smooth_at {q : ℝ → ℝ} {k alpha beta : ℝ} (hq : ContinuousAt q k)
    (hpos : 0 < q k) (ha : 0 < alpha) (hb : 0 < beta) :
    ContinuousAt (fun t => smooth (q t) alpha beta) k := by
  let e := fun t => exp 1 + 1.8*beta*q t
  let u := fun t => log (e t)
  let d := fun t => 1+69.9*powPos (q t) 1.08
  have he : ContinuousAt e k := continuousAt_const.add (continuousAt_const.mul hq)
  have e1 : 1 < e k := by
    dsimp [e]; have he' := Real.one_lt_exp_iff (x := 1) |>.2 (by norm_num)
    linarith [mul_pos (mul_pos (by norm_num : (0:ℝ) < 1.8) hb) hpos]
  have hu : ContinuousAt u k := he.log (zero_lt_one.trans e1).ne'
  have hu' : 0 < u k := log_pos e1
  have hd : ContinuousAt d k := continuousAt_const.add (continuousAt_const.mul (contPow hq hpos 1.08))
  have hd' : 0 < d k := by dsimp [d,powPos]; positivity
  have hc : ContinuousAt (fun t => 14.2/alpha + 386/d t) k :=
    continuousAt_const.add (continuousAt_const.div hd hd'.ne')
  exact hu.div (hu.add ((hc.mul hq).mul hq)) (by dsimp; positivity)

theorem transfer_continuous {a : Coefficients ℝ} (hh : 0<a.h) (hk : 0<a.kEquality)
    (hs : 0<a.soundHorizon) (hSilk : 0<a.kSilk) (ha : 0<a.alphaCDM) (hb : 0<a.betaCDM)
    (hbb : 0≤a.betaBaryon) (hN : 0≤a.betaNode) : ContinuousOn (transferReal a) (Ioi 0) := by
  apply continuousOn_of_forall_continuousAt
  intro k hx
  have h0 : 0<k := hx
  let K := fun t : ℝ => t*a.h
  let U := fun t : ℝ => acousticK a t
  let Q := fun t : ℝ => reducedK a t
  have cK : ContinuousAt K k := continuousAt_id.mul_const _
  have cU : ContinuousAt U k := cK.mul_const _
  have cQ : ContinuousAt Q k := (cK.div_const _).div_const _
  have hK : 0 < K k := mul_pos h0 hh
  have hQ : 0 < Q k := by dsimp [Q,reducedK]; positivity
  have cf : ContinuousAt (fun t => blend (U t)) k :=
    continuousAt_const.div (continuousAt_const.add ((cU.div_const _).pow 4)) (by change (1 + (U k/5.4)^4 : ℝ) ≠ 0; positivity)
  have cc : ContinuousAt (cold a) k :=
    (cf.mul (smooth_at cQ hQ (by norm_num) hb)).add
      ((continuousAt_const.sub cf).mul (smooth_at cQ hQ ha hb))
  let base := fun t => 1+(a.betaNode / K t / a.soundHorizon)^3
  have cb : ContinuousAt base k :=
    continuousAt_const.add (((continuousAt_const.div cK hK.ne').div_const _).pow 3)
  have hbase : 0 < base k := by dsimp [base]; positivity
  have cp := contPow cb hbase (1/3)
  have ph : ContinuousAt (phase a) k :=
    cK.mul (continuousAt_const.div cp (Real.exp_pos _).ne')
  have phpos : 0 < phase a k := by unfold phase powPos; positivity
  have damp : ContinuousAt (damping a) k := Real.continuous_exp.continuousAt.comp
    (continuousAt_const.sub (contPow (cK.div_const _) (div_pos hK hSilk) 1.4))
  have fb : ContinuousAt (fun t => 1+(a.betaBaryon / K t / a.soundHorizon)^3) k :=
    continuousAt_const.add (((continuousAt_const.div cK hK.ne').div_const _).pow 3)
  have bc : ContinuousAt (baryon a) k :=
    (((smooth_at cQ hQ (by norm_num) (by norm_num)).div
      (continuousAt_const.add ((cU.div_const _).pow 2)) (by dsimp; positivity)).add
        ((continuousAt_const.div fb (by dsimp; positivity)).mul damp)).mul
          ((Real.continuous_sin.continuousAt.comp ph).div ph phpos.ne')
  rw [show transferReal a = fun t => a.baryonFraction*baryon a t + a.cdmFraction*cold a t from
    funext (transferReal_eq a)]
  exact (continuousAt_const.mul bc).add (continuousAt_const.mul cc)

theorem realPower_continuous : ContinuousOn realPower (Ioi 0) := by
  have b := fiducial_bounds
  have ht := transfer_continuous (lt_of_lt_of_le (by norm_num) b.h.1)
    (lt_of_lt_of_le (by norm_num) b.keq.1) (lt_of_lt_of_le (by norm_num) b.s.1)
    (lt_of_lt_of_le (by norm_num) b.silk.1) (lt_of_lt_of_le (by norm_num) b.ac.1)
    (lt_of_lt_of_le (by norm_num) b.bc.1) b.bb.1 b.node.1
  apply ContinuousOn.mul _ (ht.pow 2)
  apply continuousOn_of_forall_continuousAt
  intro k hk
  exact continuousAt_id.rpow_const (Or.inl (show 0 < k from hk).ne')

/-- No external continuity premise is needed for the matched model. -/
theorem matched_derivative {M a b : ℝ} (hm : 0 < M) (ha : 0 < a) (hab : a ≤ b) :
    HasDerivAt (fun t => normalizedS t a b) (normalizedDerivative M a b) M :=
  massDerivative_spec hm ha hab (realPower_continuous.mono fun _ hk => ha.trans_le hk.1)

#print axioms matched_derivative
end SashimiEhContinuity
