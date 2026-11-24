import SashimiEisensteinHu
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Sinc

/-! Analytic EH98 bounds outside a quadrature band. k is in h/Mpc. Every
coefficient inequality used for the fiducial cosmology is kernel checked. -/
namespace SashimiEisensteinHu

structure TailBounds (a : Coefficients ℝ) : Prop where
  h : 0.65 ≤ a.h ∧ a.h ≤ 1
  fb : 0 ≤ a.baryonFraction ∧ a.baryonFraction ≤ 0.2
  fc : 0 ≤ a.cdmFraction ∧ a.cdmFraction ≤ 1
  keq : 0.01 ≤ a.kEquality ∧ a.kEquality ≤ 0.011
  s : 100 ≤ a.soundHorizon ∧ a.soundHorizon ≤ 200
  silk : 0.08 ≤ a.kSilk ∧ a.kSilk ≤ 0.1
  ac : 0.7 ≤ a.alphaCDM ∧ a.alphaCDM ≤ 1
  bc : 1 ≤ a.betaCDM ∧ a.betaCDM ≤ 2
  ab : 0 ≤ a.alphaBaryon ∧ a.alphaBaryon ≤ 1
  bb : 0 ≤ a.betaBaryon ∧ a.betaBaryon ≤ 8
  node : 0 ≤ a.betaNode ∧ a.betaNode ≤ 4

noncomputable def fiducial := prepareReal 0.674 0.315 0.049 2.7255
private def within? (i : Nat) (l u : ℚ) : Option Unit := do
  let bs ← coefficientBoxes? 0.674 0.315 0.049 2.7255 48 16
  let b ← bs[i]?
  guard (l ≤ b.lo ∧ b.hi ≤ u)

private theorem checked {i : Nat} {l u : ℚ} (h : within? i l u = some ()) :
    (l:ℝ) ≤ (fiducial.toList[i]?).getD 0 ∧ (fiducial.toList[i]?).getD 0 ≤ (u:ℝ) := by
  obtain ⟨bs,hc,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨b,hb,h⟩ := Option.bind_eq_some_iff.mp h
  have hw : l ≤ b.lo ∧ b.hi ≤ u := by by_contra hn; simp [guard, hn] at h
  have hi := coefficientBoxes_sound hc hb
  exact ⟨le_trans (by exact_mod_cast hw.1) hi.1, le_trans hi.2 (by exact_mod_cast hw.2)⟩

theorem fiducial_bounds : TailBounds fiducial := by
  constructor
  · exact checked (i := 0) (l := 0.65) (u := 1) (by decide +kernel)
  · exact checked (i := 1) (l := 0) (u := 0.2) (by decide +kernel)
  · exact checked (i := 2) (l := 0) (u := 1) (by decide +kernel)
  · exact checked (i := 3) (l := 0.01) (u := 0.011) (by decide +kernel)
  · exact checked (i := 4) (l := 100) (u := 200) (by decide +kernel)
  · exact checked (i := 5) (l := 0.08) (u := 0.1) (by decide +kernel)
  · exact checked (i := 6) (l := 0.7) (u := 1) (by decide +kernel)
  · exact checked (i := 7) (l := 1) (u := 2) (by decide +kernel)
  · exact checked (i := 8) (l := 0) (u := 1) (by decide +kernel)
  · exact checked (i := 9) (l := 0) (u := 8) (by decide +kernel)
  · exact checked (i := 10) (l := 0) (u := 4) (by decide +kernel)

private theorem reduced_lower {a : Coefficients ℝ} (ha : TailBounds a) {k : ℝ} (hk : 0 < k) :
    4*k ≤ reducedK a k := by
  have he := ha.keq
  have hh := ha.h
  unfold reducedK
  apply (le_div_iff₀ (by linarith : 0 < a.kEquality)).mpr
  apply (le_div_iff₀ (by norm_num : (0:ℝ) < 13.41)).mpr
  nlinarith [mul_nonneg hk.le (sub_nonneg.mpr he.2), mul_nonneg hk.le (sub_nonneg.mpr hh.1)]

private theorem smooth_bounds {q alpha beta : ℝ} (hq : 0 < q) (ha : 0 < alpha) (hb : 0 ≤ beta) :
    0 ≤ smooth q alpha beta ∧ smooth q alpha beta ≤ 1 := by
  let l := Real.log (Real.exp 1 + 1.8*beta*q)
  have hl : 0 < l := Real.log_pos (by nlinarith [Real.exp_one_gt_two, mul_nonneg hb hq.le])
  dsimp [l] at hl
  have hp : 0 < powPos q 1.08 := Real.exp_pos _
  unfold smooth
  have hterm : 0 ≤ (14.2/alpha+386/(1+69.9*powPos q 1.08))*q*q := by positivity
  constructor
  · positivity
  · exact (div_le_one (by linarith)).mpr (by linarith)

private theorem blend_bounds (u : ℝ) : 0 ≤ blend u ∧ blend u ≤ 1 := by
  unfold blend
  have h : 0 ≤ (u/5.4)^4 := by positivity
  exact ⟨by positivity,(div_le_one (by linarith)).mpr (by linarith)⟩

private theorem transfer_parts {a : Coefficients ℝ} (ha : TailBounds a) {k : ℝ} (hk : 0 < k) :
    |cold a k| ≤ 1 ∧ |baryon a k| ≤ 2 := by
  have hq := reduced_lower ha hk
  have pos : 0 < reducedK a k := by linarith
  have hb := ha.bc
  have hac := ha.ac
  have c1 := smooth_bounds pos zero_lt_one (show 0 ≤ a.betaCDM by linarith)
  have c2 := smooth_bounds pos (show 0 < a.alphaCDM by linarith) (show 0 ≤ a.betaCDM by linarith)
  have c3 := smooth_bounds pos zero_lt_one zero_le_one
  have hf := blend_bounds (acousticK a k)
  have f0 := hf.1
  have f1 := sub_nonneg.mpr hf.2
  have z1 := c1.1
  have z2 := c2.1
  have z3 := c3.1
  constructor
  · unfold cold
    rw [abs_of_nonneg (by positivity)]
    nlinarith [mul_nonneg hf.1 (sub_nonneg.mpr c1.2),
      mul_nonneg (sub_nonneg.mpr hf.2) (sub_nonneg.mpr c2.2)]
  · have hh : 0 < a.h := by linarith [ha.h.1]
    have hs : 0 < a.soundHorizon := by linarith [ha.s.1]
    have hbb := ha.bb.1
    have hab := ha.ab.1
    have hab1 := ha.ab.2
    have hr : 0 ≤ a.betaBaryon/(k*a.h)/a.soundHorizon := by positivity
    have hd : 0 < damping a k := Real.exp_pos _
    have hd1 : damping a k ≤ 1 := Real.exp_le_one_iff.mpr (by unfold powPos; linarith [Real.exp_pos (Real.log (k*a.h/a.kSilk)*(1.4:ℝ))])
    have hw : 0 < 1+(acousticK a k/5.2)^2 := by positivity
    have hv : 0 < 1+(a.betaBaryon/(k*a.h)/a.soundHorizon)^3 := by positivity
    have ht : smooth (reducedK a k) 1 1 / (1+(acousticK a k/5.2)^2) ≤ 1 :=
      (div_le_one hw).mpr (by nlinarith [sq_nonneg (acousticK a k/5.2)])
    have hu : a.alphaBaryon / (1+(a.betaBaryon/(k*a.h)/a.soundHorizon)^3) ≤ 1 :=
      (div_le_one hv).mpr (by nlinarith [pow_nonneg hr 3])
    have hphase : 0 < phase a k := by unfold phase powPos; positivity
    have sinc : |Real.sin (phase a k) / phase a k| ≤ 1 := by
      rw [abs_div, abs_of_pos hphase, div_le_one hphase]
      simpa only [abs_of_pos hphase] using Real.abs_sin_le_abs (x := phase a k)
    unfold baryon
    rw [abs_mul, abs_of_nonneg (by positivity)]
    have H := mul_le_mul (a := a.alphaBaryon / (1+(a.betaBaryon/(k*a.h)/a.soundHorizon)^3))
      hu hd1 hd.le zero_le_one
    apply (mul_le_mul_of_nonneg_left sinc (by positivity)).trans
    nlinarith

theorem low_envelope {a : Coefficients ℝ} (ha : TailBounds a) {k : ℝ} (hk : 0 < k) :
    |transferReal a k| ≤ 3 := by
  have h := transfer_parts ha hk
  rw [transferReal_eq]
  refine (abs_add_le ..).trans ?_
  rw [abs_mul, abs_mul, abs_of_nonneg ha.fb.1, abs_of_nonneg ha.fc.1]
  nlinarith [mul_le_mul ha.fb.2 h.2 (abs_nonneg _) (by norm_num),
    mul_le_mul ha.fc.2 h.1 (abs_nonneg _) zero_le_one]

/-- Logarithmic growth of the no-acoustic part is bounded at the left
endpoint of the tail. This keeps the constant proportional to log(K)/K. -/
private theorem smooth_tail {k K L q alpha beta : ℝ} (hK : 16 ≤ K) (hk : K ≤ k)
    (hL : 1 ≤ L) (hexp : 20*K ≤ Real.exp L) (hq : 4*k ≤ q)
    (ha : 0 < alpha ∧ alpha ≤ 1) (hb : 0 ≤ beta ∧ beta ≤ 2) :
    smooth q alpha beta ≤ L / (224*K*k) := by
  have k0 : 0 < k := by linarith
  have K0 : 0 < K := by linarith
  have q0 : 0 < q := by linarith
  have L0 : 0 ≤ L := by linarith
  let l := Real.log (Real.exp 1 + 1.8*beta*q)
  have y0 : 1 < Real.exp 1 + 1.8*beta*q := by
    nlinarith [Real.exp_one_gt_two, mul_nonneg hb.1 q0.le]
  have l0 : 0 < l := Real.log_pos y0
  have hlog : 4*K*l ≤ L*q := by
    have hm : 1 ≤ q/(4*K) := (one_le_div (by positivity)).mpr (by linarith)
    have bound : Real.exp 1 + 1.8*beta*q ≤ 5*q := by
      nlinarith [Real.exp_one_lt_three, mul_le_mul_of_nonneg_right hb.2 q0.le]
    have hre : 5*q = 20*K*(q/(4*K)) := by field_simp; ring
    have H : l ≤ L * (q/(4*K)) := by
      calc l ≤ Real.log (5*q) := Real.log_le_log (by linarith) bound
           _ = Real.log (20*K) + Real.log (q/(4*K)) := by rw [hre, Real.log_mul (by positivity) (by positivity)]
           _ ≤ L + (q/(4*K)-1) := add_le_add ((Real.log_le_iff_le_exp (by positivity)).mpr hexp)
                (Real.log_le_sub_one_of_pos (by positivity))
           _ ≤ _ := by nlinarith [mul_nonneg (sub_nonneg.mpr hL) (sub_nonneg.mpr hm)]
    calc _ = l * (4*K) := mul_comm _ _
         _ ≤ L * (q/(4*K)) * (4*K) := mul_le_mul_of_nonneg_right H (by positivity)
         _ = _ := by field_simp
  have c14 : 14 ≤ 14.2/alpha := (le_div_iff₀ ha.1).mpr (by linarith [ha.2])
  have p0 : 0 < powPos q 1.08 := Real.exp_pos _
  let d := l + (14.2/alpha+386/(1+69.9*powPos q 1.08))*q*q
  have hd : 14*q*q ≤ d := by
    have H : 14 ≤ 14.2/alpha + 386/(1+69.9*powPos q 1.08) := by
      exact c14.trans (le_add_of_nonneg_right (by positivity))
    dsimp [d]
    linarith [mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right H q0.le) q0.le]
  have d0 : 0 < d := lt_of_lt_of_le (by positivity) hd
  change l / d ≤ _
  apply (div_le_div_iff₀ d0 (by positivity)).mpr
  calc l*(224*K*k) = (4*K*l)*(56*k) := by ring
       _ ≤ (L*q)*(14*q) := mul_le_mul hlog (by linarith) (by positivity) (by positivity)
       _ = L*(14*q*q) := by ring
       _ ≤ _ := mul_le_mul_of_nonneg_left hd L0

private theorem damping_tail {a : Coefficients ℝ} (ha : TailBounds a) {k K : ℝ}
    (hK : 16 ≤ K) (hk : K ≤ k) : damping a k ≤ 1/(9*K*k) := by
  have k0 : 0 < k := by linarith
  have K0 : 0 < K := by linarith
  have ks0 : 0 < a.kSilk := by linarith [ha.silk.1]
  let r := k*a.h/a.kSilk
  have hr : 6*k ≤ r := (le_div_iff₀ ks0).mpr (by
    nlinarith [mul_le_mul_of_nonneg_left ha.silk.2 k0.le,
               mul_le_mul_of_nonneg_left ha.h.1 k0.le])
  let y := powPos r 1.4
  have hlog : 0 ≤ Real.log r := Real.log_nonneg (by linarith)
  have hy : 6*k ≤ y := by
    calc _ ≤ r := hr
         _ = Real.exp (Real.log r) := (Real.exp_log (by linarith)).symm
         _ ≤ _ := Real.exp_le_exp.mpr (by nlinarith)
  have hh : 3*k ≤ Real.exp (y/2) := by linarith [Real.add_one_le_exp (y/2)]
  have he : 9*K*k ≤ Real.exp y := by
    calc 9*K*k ≤ (3*k)^2 := by nlinarith [mul_le_mul_of_nonneg_right hk k0.le]
         _ ≤ (Real.exp (y/2))^2 := pow_le_pow_left₀ (by positivity) hh 2
         _ = _ := by rw [sq, ← Real.exp_add]; ring_nf
  change Real.exp (0-y) ≤ _
  rw [zero_sub, Real.exp_neg, ← one_div]
  exact one_div_le_one_div_of_le (by positivity) he

theorem high_envelope {a : Coefficients ℝ} (ha : TailBounds a) {k K L : ℝ}
    (hK : 16 ≤ K) (hk : K ≤ k) (hL : 1 ≤ L) (hexp : 20*K ≤ Real.exp L) :
    |transferReal a k| ≤ (L/128 + 1/32)/(K*k) := by
  let q := reducedK a k
  let B := L/(224*K*k)
  let E := 1/(9*K*k)
  have K0 : 0 < K := by linarith
  have k0 : 0 < k := by linarith
  have L0 : 0 < L := by linarith
  have qlow := reduced_lower ha k0
  have q0 : 0 < q := by dsimp [q]; linarith
  have c : 0 < a.alphaCDM ∧ a.alphaCDM ≤ 1 := ⟨by linarith [ha.ac.1], ha.ac.2⟩
  have b : 0 ≤ a.betaCDM ∧ a.betaCDM ≤ 2 := ⟨by linarith [ha.bc.1], ha.bc.2⟩
  have hs1 := smooth_tail hK hk hL hexp qlow ⟨zero_lt_one, le_rfl⟩ b
  have hs2 := smooth_tail hK hk hL hexp qlow c b
  have hs3 := smooth_tail hK hk hL hexp qlow ⟨zero_lt_one, le_rfl⟩ ⟨zero_le_one, by norm_num⟩
  have hs1' := (smooth_bounds q0 zero_lt_one b.1).1
  have hs2' := (smooth_bounds q0 c.1 b.1).1
  have hs3' := (smooth_bounds q0 zero_lt_one zero_le_one).1
  have hf := blend_bounds (acousticK a k)
  have Hcold : |cold a k| ≤ B := by
    have f0 := hf.1
    have f1 := sub_nonneg.mpr hf.2
    unfold cold
    rw [abs_of_nonneg (add_nonneg (mul_nonneg f0 hs1') (mul_nonneg f1 hs2'))]
    dsimp [B]
    nlinarith [mul_le_mul_of_nonneg_left hs1 f0, mul_le_mul_of_nonneg_left hs2 f1]
  have Hbar : |baryon a k| ≤ B+E := by
    have hh : 0 < a.h := by linarith [ha.h.1]
    have hs : 0 < a.soundHorizon := by linarith [ha.s.1]
    have hab := ha.ab.1
    have hbb := ha.bb.1
    have hr : 0 ≤ a.betaBaryon/(k*a.h)/a.soundHorizon := by positivity
    have hd : 0 ≤ damping a k := (Real.exp_pos _).le
    have hdt := damping_tail ha hK hk
    have hw : 0 < 1+(acousticK a k/5.2)^2 := by positivity
    have hv : 0 < 1+(a.betaBaryon/(k*a.h)/a.soundHorizon)^3 := by positivity
    have ht : smooth q 1 1 / (1+(acousticK a k/5.2)^2) ≤ B := by
      refine le_trans ((div_le_self hs3' ?_)) hs3
      nlinarith [sq_nonneg (acousticK a k/5.2)]
    have hu : a.alphaBaryon / (1+(a.betaBaryon/(k*a.h)/a.soundHorizon)^3) ≤ 1 :=
      (div_le_one hv).mpr (by nlinarith [ha.ab.2, pow_nonneg hr 3])
    have hphase : 0 < phase a k := by unfold phase powPos; positivity
    have sinc : |Real.sin (phase a k) / phase a k| ≤ 1 := by
      rw [abs_div, abs_of_pos hphase, div_le_one hphase]
      simpa only [abs_of_pos hphase] using Real.abs_sin_le_abs (x := phase a k)
    change |(smooth q 1 1 / _ + _ * damping a k) * _| ≤ _
    rw [abs_mul, abs_of_nonneg (by positivity)]
    apply (mul_le_mul_of_nonneg_left sinc (by positivity)).trans
    rw [mul_one]
    exact add_le_add ht ((mul_le_mul hu hdt hd zero_le_one).trans (by simp [E]))
  rw [transferReal_eq]
  refine (abs_add_le ..).trans ?_
  rw [abs_mul, abs_mul, abs_of_nonneg ha.fb.1, abs_of_nonneg ha.fc.1]
  calc a.baryonFraction * |baryon a k| + a.cdmFraction * |cold a k| ≤ 0.2*(B+E) + 1*B :=
          add_le_add (mul_le_mul ha.fb.2 Hbar (abs_nonneg _) (by norm_num))
            (mul_le_mul ha.fc.2 Hcold (abs_nonneg _) zero_le_one)
       _ ≤ _ := by dsimp [B,E]; field_simp; nlinarith

theorem high_envelope16 {a : Coefficients ℝ} (ha : TailBounds a) {k : ℝ} (hk : 16 ≤ k) :
    |transferReal a k| ≤ (1/200)/k := by
  have he : (20:ℝ)*16 ≤ Real.exp 6 := by
    rw [show (6:ℝ) = (6:ℕ)*1 by norm_num, Real.exp_nat_mul]
    have H : (2.7:ℝ) ≤ Real.exp 1 := by linarith [Real.exp_one_gt_d9]
    exact le_trans (by norm_num : (20:ℝ)*16 ≤ 2.7^6) (pow_le_pow_left₀ (by norm_num) H 6)
  refine (high_envelope ha (K := 16) (L := 6) (by norm_num) hk (by norm_num) he).trans ?_
  have : 0 < k := by linarith
  field_simp; linarith

theorem high_envelopeBig {a : Coefficients ℝ} (ha : TailBounds a) {k : ℝ} (hk : 10^10 ≤ k) :
    |transferReal a k| ≤ (1/10^9)/k := by
  have he : (20:ℝ)*10^10 ≤ Real.exp 40 := by
    rw [show (40:ℝ) = (40:ℕ)*1 by norm_num, Real.exp_nat_mul]
    exact le_trans (by norm_num : (20:ℝ)*10^10 ≤ 2^40)
      (pow_le_pow_left₀ (by norm_num) Real.exp_one_gt_two.le 40)
  refine (high_envelope ha (K := 10^10) (L := 40) (by norm_num) hk (by norm_num) he).trans ?_
  have : 0 < k := by linarith
  field_simp; linarith

#print axioms low_envelope
#print axioms high_envelopeBig
end SashimiEisensteinHu
