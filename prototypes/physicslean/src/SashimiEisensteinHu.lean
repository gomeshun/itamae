import SashimiSpectrumExpr

/-! The with-baryon transfer fit in Colossus 1.4.0
`colossus.cosmology.power_spectrum.modelEisenstein98`, Eisenstein & Hu (1998),
ApJ 496, 605, equations (2)–(24). k is in h/Mpc, the sound horizon is in Mpc,
and the physical wavenumber inside the equations is k*h in inverse Mpc.
The amplitude T enters the linear power as k^ns * T² prior to normalization. -/
namespace SashimiEisensteinHu
open FloatLib.Numerics FloatLib.Numerics.Interval

/-- Equality and drag scales, CDM/baryon mixture and acoustic-shape coefficients.
Compute once per cosmology; no field depends on the queried wavenumber. -/
structure Coefficients (α : Type) where
  h : α
  baryonFraction : α
  cdmFraction : α
  kEquality : α
  soundHorizon : α
  kSilk : α
  alphaCDM : α
  betaCDM : α
  alphaBaryon : α
  betaBaryon : α
  betaNode : α

private structure Elementary (α : Type) where
  ofRat : ℚ → α
  pow : α → α → α
  sqrt : α → α
  log : α → α
  exp : α → α
  sin : α → α

private def prepare [Add α] [Sub α] [Mul α] [Div α] (op : Elementary α)
    (h om ob tcmb : α) : Coefficients α := Id.run do
  let c := op.ofRat
  let sq := fun x : α => x*x
  let p := fun x y => op.pow x (c y)
  let fb := ob/om
  let fc := (om-ob)/om
  let wm := om * sq h
  let wb := ob * sq h
  let theta2 := sq (tcmb/c 2.7)
  let theta4 := sq theta2
  let zeq := c 25000 * wm / theta4
  let keq := c 0.0746 * wm / theta2
  let b1d := c 0.313 * p wm (-0.419) * (c 1 + c 0.607 * p wm 0.674)
  let b2d := c 0.238 * p wm 0.223
  let zd := c 1291 * p wm 0.251 / (c 1 + c 0.659 * p wm 0.828) * (c 1 + b1d * op.pow wb b2d)
  let rd := c 31.5 * wb / theta4 / (zd/c 1000)
  let req := c 31.5 * wb / theta4 / (zeq/c 1000)
  let s := c (2/3) / keq * op.sqrt (c 6 / req) *
    op.log ((op.sqrt (c 1 + rd) + op.sqrt (rd + req)) / (c 1 + op.sqrt req))
  let ksilk := c 1.6 * p wb 0.52 * p wm 0.73 * (c 1 + p (c 10.4 * wm) (-0.95))
  let a1 := p (c 46.9 * wm) 0.670 * (c 1 + p (c 32.1 * wm) (-0.532))
  let a2 := p (c 12 * wm) 0.424 * (c 1 + p (c 45 * wm) (-0.582))
  let ac := op.pow a1 (c 0 - fb) * op.pow a2 (c 0 - fb*fb*fb)
  let b1 := c 0.944 / (c 1 + p (c 458 * wm) (-0.708))
  let b2 := p (c 0.395 * wm) (-0.0266)
  let bc := c 1 / (c 1 + b1 * (op.pow fc b2 - c 1))
  let y := (c 1 + zeq) / (c 1 + zd)
  let t := op.sqrt (c 1 + y)
  let gy := y * ((c 0 - c 6*t) + (c 2 + c 3*y) * op.log ((t+c 1)/(t-c 1)))
  let ab := c 2.07 * keq * s * p (c 1 + rd) (-3/4) * gy
  let bb := c 0.5 + fb + (c 3 - c 2*fb) * op.sqrt (sq (c 17.2 * wm) + c 1)
  let bnode := c 8.41 * p wm 0.435
  return ⟨h, fb, fc, keq, s, ksilk, ac, bc, ab, bb, bnode⟩

private def transfer [Add α] [Sub α] [Mul α] [Div α] (op : Elementary α)
    (a : Coefficients α) (k : α) (silk : Option α := none) : α := Id.run do
  let c := op.ofRat
  let sq := fun x : α => x*x
  let cube := fun x : α => x*x*x
  let physicalK := k*a.h
  let ks := physicalK*a.soundHorizon
  let q := physicalK / c 13.41 / a.kEquality
  let f := c 1 / (c 1 + sq (sq (ks / c 5.4)))
  let tail := c 386 / (c 1 + c 69.9 * op.pow q (c 1.08))
  let t0 := fun alpha beta =>
    let l := op.log (op.exp (c 1) + c 1.8*beta*q)
    l / (l + (c 14.2 / alpha + tail) * q*q)
  let tc := f * t0 (c 1) a.betaCDM + (c 1 - f) * t0 a.alphaCDM a.betaCDM
  let st := a.soundHorizon / op.pow (c 1 + cube (a.betaNode / physicalK / a.soundHorizon)) (c (1/3))
  let x := physicalK*st
  let oscillations := op.sin x / x
  let damping := match silk with
    | none => op.exp (c 0 - op.pow (physicalK / a.kSilk) (c 1.4))
    | some d => d
  let tb := (t0 (c 1) (c 1) / (c 1 + sq (ks/c 5.2)) +
    a.alphaBaryon / (c 1 + cube (a.betaBaryon / physicalK / a.soundHorizon)) * damping) * oscillations
  return a.baryonFraction*tb + a.cdmFraction*tc

private def floatOps : Elementary Float :=
  ⟨fun q => Float.ofInt q.num / Float.ofNat q.den, Float.pow, Float.sqrt, Float.log, Float.exp, Float.sin⟩

/-- Native binary64 evaluator of these equations; interval certification uses the Expr backend. -/
def prepareFloat (h om ob tcmb : Float) : Coefficients Float := prepare floatOps h om ob tcmb
def transferFloat (a : Coefficients Float) (k : Float) : Float := transfer floatOps a k

local instance : Sub Expr := ⟨.binary .sub⟩
local instance : Div Expr := ⟨.binary .div⟩
private def exprOps : Elementary Expr :=
  ⟨.const, SashimiSpectrum.positivePowExpr, .unary .sqrt, .unary .log, .unary .exp, .unary .sin⟩

def prepareExpr (h om ob tcmb : ℚ) : Coefficients Expr :=
  prepare exprOps (.const h) (.const om) (.const ob) (.const tcmb)

def transferExpr (a : Coefficients Expr) (k : Expr) : Expr := transfer exprOps a k

/-- Positive-base real powers use exp(log(x)*p), the FloatLib expression for x^p. -/
private noncomputable def realOps : Elementary ℝ :=
  ⟨fun q => q, fun x y => Real.exp (Real.log x * y), Real.sqrt, Real.log, Real.exp, Real.sin⟩

noncomputable def prepareReal (h om ob tcmb : ℝ) : Coefficients ℝ := prepare realOps h om ob tcmb
noncomputable def transferReal (a : Coefficients ℝ) (k : ℝ) : ℝ := transfer realOps a k

/-- These real formulae expose the intermediate terms for calculus and asymptotic bounds. -/
noncomputable def powPos (x y : ℝ) := Real.exp (Real.log x * y)
noncomputable def reducedK (a : Coefficients ℝ) (k : ℝ) := k * a.h / 13.41 / a.kEquality
noncomputable def acousticK (a : Coefficients ℝ) (k : ℝ) := k*a.h*a.soundHorizon
noncomputable def blend (u : ℝ) := 1/(1+(u/5.4)^4)
noncomputable def smooth (q alpha beta : ℝ) :=
  let l := Real.log (Real.exp 1 + 1.8*beta*q)
  l / (l + (14.2/alpha + 386/(1+69.9*powPos q 1.08))*q*q)
noncomputable def phase (a : Coefficients ℝ) (k : ℝ) :=
  k*a.h * (a.soundHorizon / powPos (1+(a.betaNode/(k*a.h)/a.soundHorizon)^3) (1/3))
noncomputable def damping (a : Coefficients ℝ) (k : ℝ) :=
  Real.exp (0 - powPos (k*a.h/a.kSilk) 1.4)
noncomputable def baryon (a : Coefficients ℝ) (k : ℝ) :=
  (smooth (reducedK a k) 1 1 / (1+(acousticK a k / 5.2)^2) +
    a.alphaBaryon / (1+(a.betaBaryon/(k*a.h)/a.soundHorizon)^3) * damping a k) *
      (Real.sin (phase a k) / phase a k)
noncomputable def cold (a : Coefficients ℝ) (k : ℝ) :=
  let f := blend (acousticK a k)
  f*smooth (reducedK a k) 1 a.betaCDM + (1-f)*smooth (reducedK a k) a.alphaCDM a.betaCDM
theorem transferReal_eq (a : Coefficients ℝ) (k : ℝ) :
    transferReal a k = a.baryonFraction*baryon a k + a.cdmFraction*cold a k := by
  unfold transferReal transfer baryon cold damping phase blend acousticK reducedK smooth powPos
  dsimp [realOps]
  push_cast
  norm_num [pow_succ]
  ring_nf
  simp

/-- Physical domain of the Colossus positive-baryon fitting formula. -/
def valid (h om ob tcmb : ℚ) : Prop := 0 < h ∧ 0 < tcmb ∧ 1/10^20 ≤ ob ∧ ob < om ∧ om ≤ 1
instance (h om ob tcmb : ℚ) : Decidable (valid h om ob tcmb) := by unfold valid; infer_instance

def Coefficients.toList (a : Coefficients α) : List α :=
  [a.h,a.baryonFraction,a.cdmFraction,a.kEquality,a.soundHorizon,a.kSilk,
   a.alphaCDM,a.betaCDM,a.alphaBaryon,a.betaBaryon,a.betaNode]

private def parameterVariables : Coefficients Expr :=
  ⟨.var 1,.var 2,.var 3,.var 4,.var 5,.var 6,.var 7,.var 8,.var 9,.var 10,.var 11⟩

/-- k in slot 0 followed by coefficients in Coefficients.toList order. -/
def variableTransfer : Expr := transferExpr parameterVariables (.var 0)

private def encloseList (bits degree : Nat) : List Expr → Option (List RationalInterval)
  | [] => some []
  | x :: xs => do
    let I ← x.gridEnclosure? bits degree
    let Is ← encloseList bits degree xs
    return I :: Is

def coefficientBoxes? (h om ob tcmb : ℚ) (bits degree : Nat) : Option (List RationalInterval) := do
  guard (valid h om ob tcmb)
  encloseList bits degree (prepareExpr h om ob tcmb).toList

private theorem encloseList_sound {bits degree es boxes}
    (hv : encloseList bits degree es = some boxes) {i : Nat} {I : RationalInterval} (hi : boxes[i]? = some I) :
    I.Contains (((es.map (·.eval fun _ => 0))[i]?).getD 0) := by
  induction es generalizing boxes i with
  | nil => simp [encloseList] at hv; subst boxes; simp at hi
  | cons e es ih =>
    obtain ⟨J,hJ,h⟩ := Option.bind_eq_some_iff.mp hv
    obtain ⟨cs,hcs,heq⟩ := Option.bind_eq_some_iff.mp h
    cases heq
    cases i with
    | zero => simp at hi; subst I; exact Expr.contains_gridEnclosure? hJ
    | succ i => simpa using ih hcs (by simpa using hi)

def env (k : ℝ) (a : Coefficients ℝ) : Nat → ℝ
  | 0 => k
  | n+1 => (a.toList[n]?).getD 0

theorem eval_variableTransfer (k : ℝ) (a : Coefficients ℝ) :
    variableTransfer.eval (env k a) = transferReal a k := by rfl

set_option maxRecDepth 8000 in
theorem eval_coefficients (h om ob tcmb : ℚ) :
    ((prepareExpr h om ob tcmb).toList.map (·.eval fun _ => 0)) =
      (prepareReal h om ob tcmb).toList := by rfl

theorem coefficientBoxes_sound {h om ob tcmb bits degree boxes}
    (hc : coefficientBoxes? h om ob tcmb bits degree = some boxes) {i : Nat} {I : RationalInterval}
    (hi : boxes[i]? = some I) :
    I.Contains (((prepareReal h om ob tcmb).toList[i]?).getD 0) := by
  obtain ⟨_,_,hb⟩ := Option.bind_eq_some_iff.mp hc
  simpa only [eval_coefficients] using encloseList_sound hb hi

/-- Enclose T on a full positive k interval, reusing already-enclosed constants. -/
def range? (boxes : List RationalInterval) (K : RationalInterval) (bits degree : Nat) :
    Option RationalInterval := do
  guard (0 < K.lo ∧ K.lo ≤ K.hi ∧ boxes.length = 11)
  variableTransfer.gridEnclosureWith? (K :: boxes) bits degree

theorem range_sound {a : Coefficients ℝ} {K boxes bits degree I} {k : ℝ}
    (hk : K.Contains k) (hb : ∀ (i : Nat) (J : RationalInterval), boxes[i]? = some J → J.Contains ((a.toList[i]?).getD 0))
    (h : range? boxes K bits degree = some I) : I.Contains (transferReal a k) := by
  obtain ⟨_,_,he⟩ := Option.bind_eq_some_iff.mp h
  rw [← eval_variableTransfer]
  apply Expr.contains_gridEnclosureWith? (values := env k a) (h := he)
  intro i J hj
  cases i with
  | zero => simp at hj; subst J; exact hk
  | succ i => exact hb i J (by simpa using hj)

/-- At k*h/kSilk ≥ 128, replacing the Silk exponential by the *interval*
[0,2⁻¹²⁸] avoids forming exp of extremely large arguments. All other terms,
including the acoustic phase, are still evaluated over the full input. -/
def dampedTransfer : Expr := transfer exprOps parameterVariables (.var 0) (some (.var 12))

private noncomputable def silk (a : Coefficients ℝ) (k : ℝ) : ℝ :=
  Real.exp (0 - Real.exp (Real.log (k*a.h/a.kSilk) * (1.4:ℝ)))

private def silkBox : RationalInterval := ⟨0, 1/2^128⟩

private theorem silk_small {a : Coefficients ℝ} {k : ℝ} (hr : 128 ≤ k*a.h/a.kSilk) :
    silkBox.Contains (silk a k) := by
  let r := k*a.h/a.kSilk
  have hr' : 128 ≤ r := hr
  have hlog : 0 ≤ Real.log r := Real.log_nonneg (by linarith)
  have hy : 128 ≤ Real.exp (Real.log r * 1.4) := by
    calc 128 ≤ r := hr'
         _ = Real.exp (Real.log r) := (Real.exp_log (by linarith)).symm
         _ ≤ _ := Real.exp_le_exp.mpr (by linarith)
  have hl : Real.exp (-1) ≤ (1/2:ℝ) := by
    rw [Real.exp_neg, one_div]
    exact inv_anti₀ (by norm_num) (by linarith [Real.add_one_le_exp (1:ℝ)])
  have he : Real.exp (-128) ≤ (1/2:ℝ)^128 := by
    rw [show (-128 : ℝ) = (128:ℕ)*(-1:ℝ) by norm_num, Real.exp_nat_mul]
    gcongr
  change (0:ℝ) ≤ Real.exp (0 - Real.exp (Real.log r * (1.4:ℝ))) ∧
    Real.exp (0 - Real.exp (Real.log r * (1.4:ℝ))) ≤ ((1/2^128:ℚ):ℝ)
  refine ⟨(Real.exp_pos _).le, ?_⟩
  calc _ ≤ Real.exp (-128) := Real.exp_le_exp.mpr (by linarith)
       _ ≤ _ := he.trans (by push_cast; rw [div_pow, one_pow])

private def silkRatio : Expr := .binary .div (.var 0 * .var 1) (.var 6)

/-- Integral-ready transfer enclosure. The exact analytic case split only
changes how the proof bounds the Silk factor when the exponential is tiny. -/
def rangeFast? (boxes : List RationalInterval) (K : RationalInterval) (bits degree : Nat) :
    Option RationalInterval := do
  guard (0 < K.lo ∧ K.lo ≤ K.hi ∧ boxes.length = 11)
  let r ← silkRatio.gridEnclosureWith? (K::boxes) bits degree
  if 128 ≤ r.lo then
    dampedTransfer.gridEnclosureWith? ((K::boxes) ++ [silkBox]) bits degree
  else range? boxes K bits degree

theorem rangeFast_sound {a : Coefficients ℝ} {K boxes bits degree I} {k : ℝ}
    (hk : K.Contains k) (hb : ∀ (i : Nat) (J : RationalInterval), boxes[i]? = some J →
      J.Contains ((a.toList[i]?).getD 0))
    (h : rangeFast? boxes K bits degree = some I) : I.Contains (transferReal a k) := by
  unfold rangeFast? at h
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨r,hr,h⟩ := Option.bind_eq_some_iff.mp h
  split at h
  · rename_i test
    have len : boxes.length = 11 := by by_contra hl; simp [guard, hl] at hg
    have inputs : ∀ (i : Nat) (J : RationalInterval), (K::boxes)[i]? = some J → J.Contains (env k a i) := by
      intro i J hj; cases i with
      | zero => simp at hj; subst J; exact hk
      | succ n => exact hb n J (by simpa using hj)
    have hratio := Expr.contains_gridEnclosureWith? inputs hr
    have low : (128 : ℝ) ≤ r.lo := by exact_mod_cast test
    have hd := silk_small (a := a) (k := k) (low.trans hratio.1)
    let v : Nat → ℝ := fun i => if i = 12 then silk a k else env k a i
    have inputs' : ∀ (i : Nat) (J : RationalInterval), ((K::boxes) ++ [silkBox])[i]? = some J → J.Contains (v i) := by
      intro i J hj
      by_cases hi : i < (K::boxes).length
      · rw [List.getElem?_append_left hi] at hj
        have hn : i ≠ 12 := by simp [len] at hi; omega
        simpa [v, hn] using inputs i J hj
      · have hn := Nat.le_of_not_gt hi
        rw [List.getElem?_append_right hn] at hj
        simp [len] at hj
        by_cases he : i = 12
        · subst i; simp at hj; subst J; simpa [v] using hd
        · have hp : 0 < i-12 := by simp [len] at hn; omega
          have none : [silkBox][i-12]? = none := List.getElem?_eq_none (by simp; omega)
          simp [none] at hj
    exact Expr.contains_gridEnclosureWith? inputs' h
  · exact range_sound hk hb h

def point? (h om ob tcmb k : ℚ) (bits degree : Nat) : Option RationalInterval := do
  let boxes ← coefficientBoxes? h om ob tcmb bits degree
  range? boxes (.point k) bits degree

theorem point_sound {h om ob tcmb k bits degree I}
    (hp : point? h om ob tcmb k bits degree = some I) :
    I.Contains (transferReal (prepareReal h om ob tcmb) k) := by
  obtain ⟨boxes,hb,hr⟩ := Option.bind_eq_some_iff.mp hp
  exact range_sound (RationalInterval.contains_point k) (fun _ _ => coefficientBoxes_sound hb) hr

#print axioms point_sound
end SashimiEisensteinHu
