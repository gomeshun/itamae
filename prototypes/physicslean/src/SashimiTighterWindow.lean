import SashimiEhDerivativeBounds

/-! Whole-cell refinements of the top-hat window and the signed logarithmic-mass
weight 2 W(sinc−W). All polynomials below enclose the real functions; their
coefficients and all subsequent interval operations are rational. -/
namespace SashimiTighterWindow
set_option maxHeartbeats 1200000
open Real Set FloatLib.Numerics FloatLib.Numerics.Interval SashimiSpectrum

private theorem le_of_derivatives {f g df dg : ℝ → ℝ}
    (hf : ∀ t, HasDerivAt f (df t) t) (hg : ∀ t, HasDerivAt g (dg t) t)
    (h0 : f 0 ≤ g 0) (hs : ∀ t, 0 ≤ t → df t ≤ dg t) {x : ℝ} (hx : 0≤x) : f x ≤ g x := by
  have hd (t : ℝ) := (hg t).sub (hf t)
  have hm : MonotoneOn (fun t => g t - f t) (Ici 0) :=
    monotoneOn_of_hasDerivWithinAt_nonneg (convex_Ici _)
      (fun t _ => (hd t).continuousAt.continuousWithinAt)
      (fun t _ => (hd t).hasDerivWithinAt)
      (fun t ht => sub_nonneg.mpr (hs t (interior_subset ht)))
  have hh := hm (show (0:ℝ) ∈ Ici (0:ℝ) from (by simp)) hx hx
  dsimp at hh
  linarith

private noncomputable def s3 (x : ℝ) := x - x^3/6
private noncomputable def c4 (x : ℝ) := 1 - x^2/2 + x^4/24
private theorem dc4 (x : ℝ) : HasDerivAt c4 (-s3 x) x := by
  convert! (((hasDerivAt_const x (1:ℝ)).sub (((hasDerivAt_id' x).pow 2).div_const 2)).add (((hasDerivAt_id' x).pow 4).div_const 24)) using 1
  simp [s3]; ring
private noncomputable def s5 (x : ℝ) := x - x^3/6 + x^5/120
private theorem ds5 (x : ℝ) : HasDerivAt s5 (c4 x) x := by
  convert! (((hasDerivAt_id' x).sub (((hasDerivAt_id' x).pow 3).div_const 6)).add (((hasDerivAt_id' x).pow 5).div_const 120)) using 1
  simp [c4]; ring
private noncomputable def c6 (x : ℝ) := 1 - x^2/2 + x^4/24 - x^6/720
private theorem dc6 (x : ℝ) : HasDerivAt c6 (-s5 x) x := by
  convert! ((dc4 x).sub (((hasDerivAt_id' x).pow 6).div_const 720)) using 1
  simp [s5,s3]; ring
private noncomputable def s7 (x : ℝ) := x - x^3/6 + x^5/120 - x^7/5040
private theorem ds7 (x : ℝ) : HasDerivAt s7 (c6 x) x := by
  convert! ((ds5 x).sub (((hasDerivAt_id' x).pow 7).div_const 5040)) using 1
  simp [c6,c4]; ring
private noncomputable def c8 (x : ℝ) := 1 - x^2/2 + x^4/24 - x^6/720 + x^8/40320
private theorem dc8 (x : ℝ) : HasDerivAt c8 (-s7 x) x := by
  convert! ((dc6 x).add (((hasDerivAt_id' x).pow 8).div_const 40320)) using 1
  simp [s7,s5]; ring
private noncomputable def s9 (x : ℝ) := x - x^3/6 + x^5/120 - x^7/5040 + x^9/362880
private theorem ds9 (x : ℝ) : HasDerivAt s9 (c8 x) x := by
  convert! ((ds7 x).add (((hasDerivAt_id' x).pow 9).div_const 362880)) using 1
  simp [c8,c6]; ring
private noncomputable def c10 (x : ℝ) := 1 - x^2/2 + x^4/24 - x^6/720 + x^8/40320 - x^10/3628800
private theorem dc10 (x : ℝ) : HasDerivAt c10 (-s9 x) x := by
  convert! ((dc8 x).sub (((hasDerivAt_id' x).pow 10).div_const 3628800)) using 1
  simp [s9,s7]; ring
private noncomputable def s11 (x : ℝ) := x - x^3/6 + x^5/120 - x^7/5040 + x^9/362880 - x^11/39916800
private theorem ds11 (x : ℝ) : HasDerivAt s11 (c10 x) x := by
  convert! ((ds9 x).sub (((hasDerivAt_id' x).pow 11).div_const 39916800)) using 1
  simp [c10,c8]; ring


private theorem trig_bounds {x : ℝ} (hx : 0≤x) :
    s11 x ≤ sin x ∧ sin x ≤ s9 x ∧ c10 x ≤ cos x ∧ cos x ≤ c8 x := by
  have h3 (t : ℝ) (ht : 0≤t) : s3 t ≤ sin t := sin_ge_sub_cube ht
  have h4 (x : ℝ) (hx : 0≤x) : cos x ≤ c4 x :=
    le_of_derivatives Real.hasDerivAt_cos dc4 (by norm_num [c4]) (fun t ht => neg_le_neg (h3 t ht)) hx
  have h5 (x : ℝ) (hx : 0≤x) : sin x ≤ s5 x :=
    le_of_derivatives Real.hasDerivAt_sin ds5 (by norm_num [s5]) (fun t ht => h4 t ht) hx
  have h6 (x : ℝ) (hx : 0≤x) : c6 x ≤ cos x :=
    le_of_derivatives dc6 Real.hasDerivAt_cos (by norm_num [c6]) (fun t ht => neg_le_neg (h5 t ht)) hx
  have h7 (x : ℝ) (hx : 0≤x) : s7 x ≤ sin x :=
    le_of_derivatives ds7 Real.hasDerivAt_sin (by norm_num [s7]) (fun t ht => h6 t ht) hx
  have h8 (x : ℝ) (hx : 0≤x) : cos x ≤ c8 x :=
    le_of_derivatives Real.hasDerivAt_cos dc8 (by norm_num [c8]) (fun t ht => neg_le_neg (h7 t ht)) hx
  have h9 (x : ℝ) (hx : 0≤x) : sin x ≤ s9 x :=
    le_of_derivatives Real.hasDerivAt_sin ds9 (by norm_num [s9]) (fun t ht => h8 t ht) hx
  have h10 (x : ℝ) (hx : 0≤x) : c10 x ≤ cos x :=
    le_of_derivatives dc10 Real.hasDerivAt_cos (by norm_num [c10]) (fun t ht => neg_le_neg (h9 t ht)) hx
  have h11 (x : ℝ) (hx : 0≤x) : s11 x ≤ sin x :=
    le_of_derivatives ds11 Real.hasDerivAt_sin (by norm_num [s11]) (fun t ht => h10 t ht) hx
  exact ⟨h11 _ hx,h9 _ hx,h10 _ hx,h8 _ hx⟩

/-- Horner interval evaluation. The highest-degree coefficients below are
intervals because the Taylor upper and lower polynomials differ. -/
private def horner (t : RationalInterval) : List RationalInterval → RationalInterval
  | [] => .point 0
  | c::cs => c.add (t.mul (horner t cs))

private noncomputable def poly (x : ℝ) (get : RationalInterval → ℚ) : List RationalInterval → ℝ
  | [] => 0
  | c::cs => (get c:ℝ) + x*poly x get cs

private theorem horner_sound {I : RationalInterval} {x : ℝ} (hx : I.Contains x)
    (cs : List RationalInterval) (get : RationalInterval → ℚ)
    (hg : ∀ c ∈ cs, c.Contains (get c)) : (horner I cs).Contains (poly x get cs) := by
  induction cs with
  | nil => simp [horner,poly,RationalInterval.point,RationalInterval.Contains]
  | cons c cs ih =>
    exact RationalInterval.contains_add (hg _ (by simp))
      (RationalInterval.contains_mul hx (ih (fun q hq => hg _ (by simp [hq]))))

private def wCoeffs : List RationalInterval :=
  [.point 1,.point (-1/10),.point (1/280),.point (-1/15120),⟨-1/13305600,1/1209600⟩]
private def dCoeffs : List RationalInterval :=
  [.point 0,.point (-1/15),.point (1/210),.point (-1/7560),⟨1/518400,113/39916800⟩,⟨-1/39916800,0⟩]

private theorem W_bounds {x : ℝ} (hx : 0<x) :
    poly (x^2) (·.lo) wCoeffs ≤ topHat x ∧ topHat x ≤ poly (x^2) (·.hi) wCoeffs := by
  obtain ⟨hS,hS',hC,hC'⟩ := trig_bounds hx.le
  rw [topHat, ite_eq_right hx.ne']
  constructor
  · apply (le_div_iff₀ (by positivity)).2
    norm_num [poly,wCoeffs,RationalInterval.point,s11,c8] at hS hC' ⊢
    linarith only [hS, mul_le_mul_of_nonneg_left hC' hx.le]
  · apply (div_le_iff₀ (by positivity)).2
    norm_num [poly,wCoeffs,RationalInterval.point,s9,c10] at hS' hC ⊢
    linarith only [hS', mul_le_mul_of_nonneg_left hC hx.le]

private theorem D_bounds {x : ℝ} (hx : 0<x) :
    poly (x^2) (·.lo) dCoeffs ≤ sin x/x - topHat x ∧ sin x/x - topHat x ≤ poly (x^2) (·.hi) dCoeffs := by
  obtain ⟨hl,hu,_,_⟩ := trig_bounds hx.le
  have hL := div_le_div_of_nonneg_right hl hx.le
  have hU := div_le_div_of_nonneg_right hu hx.le
  have w := W_bounds hx
  constructor
  · apply le_trans _ (sub_le_sub hL w.2)
    norm_num [poly,dCoeffs,wCoeffs,RationalInterval.point,s11]
    apply le_of_eq; field_simp [hx.ne']; ring
  · apply le_trans (sub_le_sub hU w.1) _
    norm_num [poly,dCoeffs,wCoeffs,RationalInterval.point,s9]
    apply le_of_eq; field_simp [hx.ne']; ring

private def taylor (cs : List RationalInterval) (l u : ℚ) :=
  horner ⟨l^2,u^2⟩ cs

private theorem taylor_sound {l u : ℚ} {x y : ℝ} (hl : 0≤l)
    (hx : x ∈ Icc (l:ℝ) u) (cs : List RationalInterval)
    (hg : ∀ c∈cs, c.lo ≤ c.hi)
    (h : poly (x^2) (·.lo) cs ≤ y ∧ y ≤ poly (x^2) (·.hi) cs) :
    (taylor cs l u).Contains y := by
  have hl' : (0:ℝ) ≤ l := by exact_mod_cast hl
  have hI : (⟨l^2,u^2⟩ : RationalInterval).Contains (x^2) := by
    constructor <;> push_cast
    · exact pow_le_pow_left₀ hl' hx.1 _
    · exact pow_le_pow_left₀ (hl'.trans hx.1) hx.2 _
  refine ⟨(horner_sound hI cs (·.lo) ?_).1.trans h.1,h.2.trans (horner_sound hI cs (·.hi) ?_).2⟩
  · intro c hc; exact ⟨le_rfl, by exact_mod_cast hg c hc⟩
  · intro c hc; exact ⟨by exact_mod_cast hg c hc,le_rfl⟩

private theorem w_order : ∀ c ∈ wCoeffs, c.lo ≤ c.hi := by decide +kernel
private theorem d_order : ∀ c ∈ dCoeffs, c.lo ≤ c.hi := by decide +kernel

private def meet (A B : RationalInterval) : RationalInterval := ⟨max A.lo B.lo,min A.hi B.hi⟩
private theorem meet_sound {A B : RationalInterval} {x : ℝ} (hA : A.Contains x) (hB : B.Contains x) :
    (meet A B).Contains x := by simpa [meet,RationalInterval.Contains] using And.intro (max_le hA.1 hB.1) (le_min hA.2 hB.2)

private def smallSq (l u : ℚ) : RationalInterval :=
  let W := taylor wCoeffs l u
  W.mul W

private def smallWeight (l u : ℚ) : RationalInterval :=
  ((taylor wCoeffs l u).mul (taylor dCoeffs l u)).scaleNonnegative 2

private theorem smallSq_sound {l u : ℚ} {x : ℝ} (hl : 0≤l)
    (hx : x ∈ Icc (l:ℝ) u) : (smallSq l u).Contains (topHat x ^ 2) := by
  have W : (taylor wCoeffs l u).Contains (topHat x) := by
    apply taylor_sound hl hx _ w_order
    rcases ((by exact_mod_cast hl : (0:ℝ)≤l).trans hx.1).eq_or_lt with e|p
    · simp [← e,poly,wCoeffs,RationalInterval.point,topHat]
    · exact W_bounds p
  simpa [smallSq, sq] using RationalInterval.contains_mul W W

private theorem smallWeight_sound {l u : ℚ} {x : ℝ} (hl : 0<l) (hx : x ∈ Icc (l:ℝ) u) :
    (smallWeight l u).Contains (SashimiEhDerivativeBounds.shape x) := by
  have hp := (by exact_mod_cast hl : (0:ℝ)<l).trans_le hx.1
  have hw := taylor_sound hl.le hx _ w_order (W_bounds hp)
  have hd := taylor_sound hl.le hx _ d_order (D_bounds hp)
  exact_mod_cast RationalInterval.contains_scaleNonnegative (RationalInterval.contains_mul hw hd) (by decide : (0:ℚ)≤2)

/-- Value and ordinary x derivative over the same cell. Algebraically,
W = 3 (sinc − cos)/x² and W' = 3 (sinc − W)/x. -/
private structure Jet where
  val : RationalInterval
  slope : RationalInterval

private def jet (signed : Bool) (l u : ℚ) (si cs : RationalInterval) : Jet := Id.run do
  let X : RationalInterval := ⟨l,u⟩
  let S := si.divPositive X
  let W := ((S.sub cs).scaleNonnegative 3).divPositive ⟨l^2,u^2⟩
  let D := S.sub W
  let dW := (D.scaleNonnegative 3).divPositive X
  if signed then
    let dS := (cs.sub S).divPositive X
    return ⟨(W.mul D).scaleNonnegative 2,
      ((dW.mul D).add (W.mul (dS.sub dW))).scaleNonnegative 2⟩
  return ⟨W.mul W,(W.mul dW).scaleNonnegative 2⟩

private noncomputable def value (signed : Bool) (x : ℝ) :=
  if signed then SashimiEhDerivativeBounds.shape x else topHat x^2

private theorem jet_sound {sg : Bool} {l u : ℚ} {x : ℝ} {si cs : RationalInterval}
    (hl : 0<l) (hx : x ∈ Icc (l:ℝ) u) (hsi : si.Contains (sin x)) (hcs : cs.Contains (cos x)) :
    (jet sg l u si cs).val.Contains (value sg x) ∧
    (jet sg l u si cs).slope.Contains (deriv (value sg) x) ∧ DifferentiableAt ℝ (value sg) x := by
  let X : RationalInterval := ⟨l,u⟩
  have hX : X.Contains x := hx
  have hl' : (0:ℝ)<l := by exact_mod_cast hl
  have hp := hl'.trans_le hx.1
  have hS := RationalInterval.contains_divPositive hsi hX hl
  have hX2 : (⟨l^2,u^2⟩:RationalInterval).Contains (x^2) := by
    constructor <;> push_cast
    · exact pow_le_pow_left₀ hl'.le hx.1 _
    · exact pow_le_pow_left₀ hp.le hx.2 _
  have W := RationalInterval.contains_divPositive
    (RationalInterval.contains_scaleNonnegative (RationalInterval.contains_sub hS hcs) (by decide : (0:ℚ)≤3))
    hX2 (by dsimp; positivity)
  have eqW : ((3:ℚ):ℝ)*(sin x/x-cos x)/x^2 = topHat x := by
    simp only [topHat,ite_eq_right hp.ne',Rat.cast_ofNat]; field_simp
  rw [eqW] at W
  have D := RationalInterval.contains_sub hS W
  have dW := RationalInterval.contains_divPositive
    (RationalInterval.contains_scaleNonnegative D (by decide : (0:ℚ)≤3)) hX hl
  have eqP : ((3:ℚ):ℝ)*(sin x/x-topHat x)/x=topHatPrime x := by
    simp only [topHat,topHatPrime,ite_eq_right hp.ne',Rat.cast_ofNat]; field_simp; ring
  rw [eqP] at dW
  have hW := topHat_hasDerivAt x hp
  cases sg with
  | false =>
    have hv := hW.pow 2
    refine ⟨?_, ?_, hv.differentiableAt⟩
    · rw [show value false x = topHat x*topHat x by simp [value,sq]]
      exact RationalInterval.contains_mul W W
    · change (jet false l u si cs).slope.Contains (deriv (topHat^2) x)
      rw [hv.deriv]
      change (jet false l u si cs).slope.Contains (2*topHat x^1*topHatPrime x)
      rw [pow_one,mul_assoc]
      have hd := RationalInterval.contains_scaleNonnegative (RationalInterval.contains_mul W dW) (by decide : (0:ℚ)≤2)
      simp only [Rat.cast_ofNat] at hd
      exact hd
  | true =>
    have hds : HasDerivAt (fun y => sin y/y) ((cos x-sin x/x)/x) x := by
      apply ((Real.hasDerivAt_sin x).div (hasDerivAt_id' x) hp.ne').congr_deriv
      field_simp
    have hd := ((hW.mul (hds.sub hW)).const_mul 2)
    have hd' : HasDerivAt (value true)
        (2*(topHatPrime x*(sin x/x-topHat x)+topHat x*((cos x-sin x/x)/x-topHatPrime x))) x := hd
    have dS := RationalInterval.contains_divPositive (RationalInterval.contains_sub hcs hS) hX hl
    refine ⟨?_,?_,?_⟩
    · exact_mod_cast RationalInterval.contains_scaleNonnegative (RationalInterval.contains_mul W D) (by decide : (0:ℚ)≤2)
    · rw [hd'.deriv]
      exact_mod_cast RationalInterval.contains_scaleNonnegative (RationalInterval.contains_add
        (RationalInterval.contains_mul dW D) (RationalInterval.contains_mul W (RationalInterval.contains_sub dS dW)))
        (by decide : (0:ℚ)≤2)
    · exact hd.differentiableAt

private def jet? (sg : Bool) (l u : ℚ) (bits degree : Nat) : Option Jet := do
  let si ← SashimiNumerics.exprRange (Expr.unary .sin (.var 0)) bits degree l u
  let cs ← SashimiNumerics.exprRange (Expr.unary .cos (.var 0)) bits degree l u
  return jet sg l u si cs

private theorem jet?_sound {sg l u bits degree j} (h : jet? sg l u bits degree = some j)
    (hl : 0<l) {x : ℝ} (hx : x ∈ Icc (l:ℝ) u) : j.val.Contains (value sg x) ∧
    j.slope.Contains (deriv (value sg) x) ∧ DifferentiableAt ℝ (value sg) x := by
  unfold jet? at h
  obtain ⟨si,hsi,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨cs,hcs,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  have hs : si.Contains (sin x) := SashimiNumerics.exprRange_sound (.unary .sin (.var 0)) bits degree _ _ _ hsi x hx
  have hc : cs.Contains (cos x) := SashimiNumerics.exprRange_sound (.unary .cos (.var 0)) bits degree _ _ _ hcs x hx
  exact jet_sound (sg:=sg) hl hx hs hc

/-- A sign-definite derivative makes the endpoint values sharp bounds. On a
cell crossing a stationary point, intersect with the mean-value estimate from
the midpoint. The two trigonometric data ranges are shared in each jet. -/
private def middle? (sg : Bool) (l u : ℚ) (bits degree : Nat) : Option RationalInterval := do
  let j ← jet? sg l u bits degree
  if 0≤j.slope.lo ∨ j.slope.hi≤0 then
    let a ← jet? sg l l bits degree
    let b ← jet? sg u u bits degree
    return meet j.val ⟨min a.val.lo b.val.lo,max a.val.hi b.val.hi⟩
  else
    let mid := (l+u)/2
    let m ← jet? sg mid mid bits degree
    let e := max (-j.slope.lo) j.slope.hi * ((u-l)/2)
    return meet j.val ⟨m.val.lo-e,m.val.hi+e⟩

private theorem middle_sound {sg l u bits degree out} (h : middle? sg l u bits degree = some out)
    (hl : 0<l) {x : ℝ} (hx : x ∈ Icc (l:ℝ) u) : out.Contains (value sg x) := by
  unfold middle? at h
  obtain ⟨j,hj,h⟩ := Option.bind_eq_some_iff.mp h
  have H := fun (t:ℝ) (ht:t ∈ Icc (l:ℝ) u) => jet?_sound hj hl ht
  have hLU := hx.1.trans hx.2
  split at h
  · rename_i hsign
    obtain ⟨a,ha,h⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨b,hb,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    apply meet_sound (H _ hx).1
    have hu : 0<u := hl.trans_le (by exact_mod_cast hLU)
    have A := (jet?_sound ha hl (show (l:ℝ) ∈ Icc (l:ℝ) l from ⟨le_rfl,le_rfl⟩)).1
    have B := (jet?_sound hb hu (show (u:ℝ) ∈ Icc (u:ℝ) u from ⟨le_rfl,le_rfl⟩)).1
    have hd (t : ℝ) (ht : t∈Icc (l:ℝ) u) := (H t ht).2.2
    have hc : ContinuousOn (value sg) (Icc (l:ℝ) u) := fun t ht => (hd t ht).continuousAt.continuousWithinAt
    have hdiff : DifferentiableOn ℝ (value sg) (interior (Icc (l:ℝ) u)) :=
      fun t ht => (hd t (interior_subset ht)).differentiableWithinAt
    constructor <;> push_cast
    all_goals rcases hsign with hs|hs
    all_goals first
      | have hm := monotoneOn_of_deriv_nonneg (convex_Icc (l:ℝ) u) hc hdiff
          (fun t ht => (by exact_mod_cast hs : (0:ℝ)≤j.slope.lo).trans (H t (interior_subset ht)).2.1.1)
      | have hm := antitoneOn_of_deriv_nonpos (convex_Icc (l:ℝ) u) hc hdiff
          (fun t ht => (H t (interior_subset ht)).2.1.2.trans (by exact_mod_cast hs))
    · exact (min_le_left _ _).trans (A.1.trans (hm ⟨le_rfl,hLU⟩ hx hx.1))
    · exact (min_le_right _ _).trans (B.1.trans (hm hx ⟨hLU,le_rfl⟩ hx.2))
    · exact (hm hx ⟨hLU,le_rfl⟩ hx.2).trans (B.2.trans (le_max_right _ _))
    · exact (hm ⟨le_rfl,hLU⟩ hx hx.1).trans (A.2.trans (le_max_left _ _))
  · obtain ⟨m,hm,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    apply meet_sound (H _ hx).1
    let mid : ℚ := (l+u)/2
    have hmR : (mid:ℝ) = ((l:ℝ)+(u:ℝ))/2 := by simp [mid]
    have hmI : (mid:ℝ) ∈ Icc (l:ℝ) u := by rw [hmR]; constructor <;> linarith
    have hp : 0<mid := hl.trans_le (by exact_mod_cast hmI.1)
    have V := (jet?_sound hm hp (show (mid:ℝ) ∈ Icc (mid:ℝ) mid from ⟨le_rfl,le_rfl⟩)).1
    let C : ℝ := max (-(j.slope.lo:ℝ)) (j.slope.hi:ℝ)
    have bound (t : ℝ) (ht : t∈Icc (l:ℝ) u) : ‖deriv (value sg) t‖ ≤ C := by
      rw [Real.norm_eq_abs,abs_le]
      exact ⟨(neg_le_neg (le_max_left _ _)).trans (by simpa using (H t ht).2.1.1),
        (H t ht).2.1.2.trans (le_max_right _ _)⟩
    have hC : 0≤C := (norm_nonneg _).trans (bound x hx)
    have hv := (convex_Icc (l:ℝ) u).norm_image_sub_le_of_norm_hasDerivWithin_le
      (fun t ht => (H t ht).2.2.hasDerivAt.hasDerivWithinAt) bound hmI hx
    simp only [Real.norm_eq_abs] at hv
    have dx : |x-(mid:ℝ)| ≤ ((u:ℝ)-(l:ℝ))/2 := by rw [hmR,abs_le]; constructor <;> linarith [hx.1,hx.2]
    have hv' := abs_le.mp (hv.trans (mul_le_mul_of_nonneg_left dx hC))
    constructor <;> push_cast <;> dsimp [C] at hv' <;> linarith [V.1,V.2]

private def tighten (sg : Bool) (old : RationalInterval) (l u : ℚ) (bits degree : Nat) :=
  if u≤2 then meet old (if sg then smallWeight l u else smallSq l u)
  else if 0<l ∧ u≤128 then
    match middle? sg l u bits degree with
    | some I => meet old I
    | none => old
  else old

private theorem tighten_refines (sg old l u bits degree) :
    old.lo ≤ (tighten sg old l u bits degree).lo ∧ (tighten sg old l u bits degree).hi ≤ old.hi := by
  unfold tighten
  split
  · exact ⟨le_max_left ..,min_le_left ..⟩
  · split
    · split
      · exact ⟨le_max_left ..,min_le_left ..⟩
      · exact ⟨le_rfl,le_rfl⟩
    · exact ⟨le_rfl,le_rfl⟩

private theorem tighten_sound {sg : Bool} {old : RationalInterval} {l u : ℚ} {bits degree : Nat} {x : ℝ} (hx : x∈Icc (l:ℝ) u)
    (hold : old.Contains (value sg x))
    (hsmall : (if sg then smallWeight l u else smallSq l u).Contains (value sg x)) :
    (tighten sg old l u bits degree).Contains (value sg x) := by
  unfold tighten
  split
  · exact meet_sound hold hsmall
  · split
    · rename_i hdom
      split
      · rename_i I hI; exact meet_sound hold (middle_sound hI hdom.1 hx)
      · exact hold
    · exact hold

/-- Refined entire-cell range of the squared spherical top-hat. -/
def windowSq (l u : ℚ) (bits degree : Nat) :=
  tighten false (SashimiTableWindow.windowSq l u bits degree) l u bits degree

theorem windowSq_sound {x : ℝ} {l u : ℚ} (hl : 0 ≤ l) (hx : (l:ℝ) ≤ x)
    (hu : x ≤ (u:ℝ)) (bits degree : Nat) : (windowSq l u bits degree).Contains ((topHat x)^2) :=
  tighten_sound (sg:=false) (old:=SashimiTableWindow.windowSq l u bits degree)
    (bits:=bits) (degree:=degree) (l:=l) (u:=u) (x:=x) ⟨hx,hu⟩ (SashimiTableWindow.windowSq_sound hl hx hu bits degree)
    (smallSq_sound hl ⟨hx,hu⟩)

theorem windowSq_refines (l u : ℚ) (bits degree : Nat) :
    (SashimiTableWindow.windowSq l u bits degree).lo ≤ (windowSq l u bits degree).lo ∧
    (windowSq l u bits degree).hi ≤ (SashimiTableWindow.windowSq l u bits degree).hi :=
  tighten_refines false (SashimiTableWindow.windowSq l u bits degree) l u bits degree

/-- Refined signed logarithmic-mass weight 2W(sinc−W). -/
def weight (l u : ℚ) (bits degree : Nat) :=
  tighten true (SashimiEhDerivativeBounds.weight l u bits degree) l u bits degree

theorem weight_sound {l u : ℚ} {x : ℝ} (hl : 0 < l) (hx : x ∈ Icc (l:ℝ) u) (bits degree : Nat) :
    (weight l u bits degree).Contains (SashimiEhDerivativeBounds.shape x) :=
  tighten_sound (sg:=true) (old:=SashimiEhDerivativeBounds.weight l u bits degree)
    (bits:=bits) (degree:=degree) (l:=l) (u:=u) (x:=x) hx (SashimiEhDerivativeBounds.weight_sound hl hx bits degree)
    (smallWeight_sound hl hx)

theorem weight_refines (l u : ℚ) (bits degree : Nat) :
    (SashimiEhDerivativeBounds.weight l u bits degree).lo ≤ (weight l u bits degree).lo ∧
    (weight l u bits degree).hi ≤ (SashimiEhDerivativeBounds.weight l u bits degree).hi :=
  tighten_refines true (SashimiEhDerivativeBounds.weight l u bits degree) l u bits degree

#print axioms windowSq_sound
#print axioms weight_sound

end SashimiTighterWindow
