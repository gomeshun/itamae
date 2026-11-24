import SashimiEhSweep

/-! A second representation of the Eisenstein–Hu spectrum on a complete cell.
Cancelling the leading powers of wavenumber before interval evaluation preserves
the correlation between the spectral tilt and the decaying CDM transfer. -/
namespace SashimiEhTighterPower
open SashimiEisensteinHu SashimiEisensteinIntegral SashimiEhSweep
open FloatLib.Numerics FloatLib.Numerics.Interval Set Real

local instance : Div Expr := ⟨.binary .div⟩
local instance : Sub Expr := ⟨.binary .sub⟩
private abbrev c := Expr.const
private abbrev power := SashimiSpectrum.positivePowExpr
private def sqr (x : Expr) := x*x
private def cube (x : Expr) := x*x*x

private theorem scaled_ratio (k v L C : ℝ) (hk : k ≠ 0) :
    L / (L / k^2 + C*v*v) = k^2 * (L / (L + C*(k*v)*(k*v))) := by
  rw [show L / k^2 + C*v*v = (L + C*(k*v)*(k*v))/k^2 by field_simp]
  rw [div_div_eq_mul_div]
  ring

private def scaledSmooth (k v alpha beta : Expr) : Expr :=
  let q := k*v
  let L := .unary .log (.unary .exp (c 1) + c 1.8*beta*q)
  let C := c 14.2/alpha + c 386/(c 1 + c 69.9*power q (c 1.08))
  L / (L / sqr k + C*v*v)

private theorem evAdd (a b : Expr) (f : Nat → ℝ) : (a+b).eval f = a.eval f + b.eval f := rfl
private theorem evMul (a b : Expr) (f : Nat → ℝ) : (a*b).eval f = a.eval f * b.eval f := rfl
private theorem evDiv (a b : Expr) (f : Nat → ℝ) : (a/b).eval f = a.eval f / b.eval f := rfl
private theorem evSub (a b : Expr) (f : Nat → ℝ) : (a-b).eval f = a.eval f - b.eval f := rfl
private theorem evC (r : ℚ) (f : Nat → ℝ) : (c r).eval f = (r:ℝ) := rfl
private theorem evVar (i : Nat) (f : Nat → ℝ) : (Expr.var i).eval f = f i := rfl
private theorem evSin (e : Expr) (f : Nat → ℝ) : (Expr.unary .sin e).eval f = sin (e.eval f) := rfl
private theorem evLog (e : Expr) (f : Nat → ℝ) : (Expr.unary .log e).eval f = log (e.eval f) := rfl
private theorem evExp (e : Expr) (f : Nat → ℝ) : (Expr.unary .exp e).eval f = exp (e.eval f) := rfl
private theorem evPow (a b : Expr) (f : Nat → ℝ) : (power a b).eval f = powPos (a.eval f) (b.eval f) := rfl

private theorem eval_scaledSmooth (k v alpha beta : Expr) (f : Nat → ℝ)
    (hk : k.eval f ≠ 0) :
    (scaledSmooth k v alpha beta).eval f =
      (k.eval f)^2 * smooth (k.eval f * v.eval f) (alpha.eval f) (beta.eval f) := by
  simp only [scaledSmooth, sqr, smooth, evLog, evExp, evAdd, evMul, evDiv, evC, evPow]
  push_cast
  rw [← sq]
  apply scaled_ratio _ _ _ _ hk

/-- Variable 0 is the independently enclosed Silk factor, 1 is k, and the
remaining eleven variables follow the public EH coefficient order. -/
private def scaledTransfer : Expr := Id.run do
  let D := Expr.var 0
  let k := Expr.var 1
  let ph := k * .var 2
  let u := ph * .var 6
  let v := .var 2 / c 13.41 / .var 5
  let z := u / c 5.4
  let mix := c 1 / (c 1 + sqr (sqr z))
  let tc := mix * scaledSmooth k v (c 1) (.var 9) +
    (c 1-mix) * scaledSmooth k v (.var 8) (.var 9)
  let x := ph * (.var 6 / power (c 1 + cube (.var 12 / ph / .var 6)) (c (1/3)))
  let tb := (scaledSmooth k v (c 1) (c 1) / (c 1 + sqr (u / c 5.2)) +
    (sqr k * .var 10) / (c 1 + cube (.var 11 / ph / .var 6)) * D) * (.unary .sin x / x)
  return .var 3*tb + .var 4*tc

private noncomputable def withDecay (a : Coefficients ℝ) (k : ℝ) : Nat → ℝ
  | 0 => damping a k
  | n+1 => SashimiEisensteinHu.env k a n

private theorem eval_scaledTransfer (a : Coefficients ℝ) {k : ℝ} (hk : 0 < k) :
    scaledTransfer.eval (withDecay a k) = k^2 * transferReal a k := by
  have hk' : (Expr.var 1).eval (withDecay a k) ≠ 0 := hk.ne'
  have hs (x y : Expr) := eval_scaledSmooth (Expr.var 1) (Expr.var 2 / c 13.41 / Expr.var 5) x y (withDecay a k) hk'
  dsimp only [scaledTransfer, Id.run, pure, Pure.pure]
  simp only [evAdd, evSub, evMul, evDiv, evC, evPow, evSin, hs, sqr, cube, evVar,
    withDecay, SashimiEisensteinHu.env, Coefficients.toList]
  rw [transferReal_eq]
  unfold cold baryon acousticK reducedK blend phase
  norm_num only [List.getElem?_cons_zero, List.getElem?_cons_succ, Option.getD_some, Rat.cast_ofScientific,
    Rat.cast_one, Rat.cast_div, Rat.cast_ofNat, pow_succ, pow_zero, mul_one]
  ring_nf

private def exponent : Expr := power (.var 0 * .var 1 / .var 6) (c 1.4)

private def decay? (K : RationalInterval) (boxes : List RationalInterval) (bits degree : Nat) :
    Option RationalInterval := do
  let t ← exponent.gridEnclosureWith? (K::boxes) bits degree
  if 128 ≤ t.lo then return ⟨0,1/2^128⟩
  else SashimiNumerics.exprRange (.unary .exp (.const 0 - .var 0)) bits degree t.lo t.hi

private theorem tiny {t : ℝ} (ht : 128 ≤ t) :
    (⟨0,1/2^128⟩ : RationalInterval).Contains (exp (0-t)) := by
  have he : exp (-1) ≤ (1/2:ℝ) := by
    rw [exp_neg, one_div]
    exact inv_anti₀ (by norm_num) (by linarith [add_one_le_exp (1:ℝ)])
  constructor
  · exact (exp_pos _).le
  · calc exp (0-t) ≤ exp ((128:ℕ)*(-1:ℝ)) := exp_le_exp.mpr (by norm_num; linarith)
         _ = _ := exp_nat_mul _ _
         _ ≤ (1/2:ℝ)^128 := by gcongr
         _ = _ := by norm_num

private theorem decay_sound {a : Coefficients ℝ} {k : ℝ} {K : RationalInterval} {boxes : List RationalInterval} {bits degree D}
    (inputs : ∀ (i : Nat) (J : RationalInterval), (K::boxes)[i]? = some J → J.Contains (SashimiEisensteinHu.env k a i))
    (h : decay? K boxes bits degree = some D) : D.Contains (damping a k) := by
  unfold decay? at h
  obtain ⟨t, ht, h⟩ := Option.bind_eq_some_iff.mp h
  have hh := Expr.contains_gridEnclosureWith? inputs ht
  change t.Contains (powPos (k*a.h/a.kSilk) 1.4) at hh
  split at h
  · rename_i large
    cases h
    exact tiny ((show (128:ℝ) ≤ t.lo by exact_mod_cast large).trans hh.1)
  · exact SashimiNumerics.exprRange_sound _ _ _ _ _ _ h _ hh

/-- Evaluate k²T and combine it with the exact compensating tilt. -/
def scaledWeight? (csm : Cosmology) (boxes : List RationalInterval) (l u : ℚ) (bits degree : Nat) :
    Option RationalInterval := do
  guard (0 < l ∧ l ≤ u)
  let K : RationalInterval := ⟨l,u⟩
  let D ← decay? K boxes bits degree
  let T ← scaledTransfer.gridEnclosureWith? (D::K::boxes) bits degree
  let V ← SashimiNumerics.exprRange (.unary .exp (.const (csm.ns-2) * .unary .log (.var 0))) bits degree l u
  let Q := T.mul T
  return V.mul ⟨max 0 Q.lo, Q.hi⟩

theorem scaledWeight_sound {csm boxes cb cd bits degree l u P}
    (hb : boxes? csm cb cd = some boxes) (h : scaledWeight? csm boxes l u bits degree = some P) :
    0<l ∧ l≤u ∧ ∀ k ∈ Icc (l:ℝ) u, P.Contains (weightReal csm k) := by
  unfold scaledWeight? at h
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  have hp : 0 < l ∧ l ≤ u := by by_contra hn; simp [guard, hn] at hg
  obtain ⟨D,hd,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨T,ht,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨V,hv,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  refine ⟨hp.1,hp.2,?_⟩
  intro k hk
  have hk' : 0 < k := (show 0 < (l:ℝ) by exact_mod_cast hp.1).trans_le hk.1
  let a := coefficients csm
  let K : RationalInterval := ⟨l,u⟩
  have hc : ∀ (i : Nat) (J : RationalInterval), (K::boxes)[i]? = some J → J.Contains (SashimiEisensteinHu.env k a i) := by
    intro i J hJ
    cases i with
    | zero => simp at hJ; subst J; exact hk
    | succ n => exact coefficientBoxes_sound hb (by simpa using hJ)
  have hd' := decay_sound hc hd
  have hT := Expr.contains_gridEnclosureWith? (values := withDecay a k) (h := ht) (by
    intro i J hJ
    cases i with
    | zero => simp at hJ; subst J; exact hd'
    | succ n => exact hc n J (by simpa using hJ))
  rw [eval_scaledTransfer a hk'] at hT
  have hs := RationalInterval.contains_mul hT hT
  rw [← sq] at hs
  have sqbox : (⟨max 0 (T.mul T).lo,(T.mul T).hi⟩ : RationalInterval).Contains
      ((k^2*transferReal a k)^2) :=
    ⟨by simpa only [Rat.cast_max,Rat.cast_zero] using max_le (sq_nonneg _) hs.1,hs.2⟩
  have hV := SashimiNumerics.exprRange_sound _ _ _ _ _ _ hv k hk
  change V.Contains (exp (((csm.ns-2:ℚ):ℝ)*log k)) at hV
  have hw := RationalInterval.contains_mul hV sqbox
  unfold weightReal
  have he : (((2+csm.ns:ℚ):ℝ)*log k) = (((csm.ns-2:ℚ):ℝ)*log k) + (4:ℕ)*log k := by
    push_cast; ring
  rw [he, exp_add, exp_nat_mul, exp_log hk']
  convert hw using 1; ring

/-- Tighten an existing enclosure using the algebraically scaled expression.
An unavailable scaled enclosure leaves the supplied endpoints usable. -/
def refine (csm : Cosmology) (boxes : List RationalInterval) (l u : ℚ) (bits degree : Nat)
    (original : RationalInterval) : RationalInterval :=
  match scaledWeight? csm boxes l u bits degree with
  | some J => ⟨max original.lo J.lo, min original.hi J.hi⟩
  | none => original

theorem refine_nested (csm : Cosmology) (boxes : List RationalInterval) (l u : ℚ)
    (bits degree : Nat) (A : RationalInterval) :
    A.lo ≤ (refine csm boxes l u bits degree A).lo ∧
      (refine csm boxes l u bits degree A).hi ≤ A.hi := by
  unfold refine
  split
  · exact ⟨le_max_left _ _, min_le_left _ _⟩
  · exact ⟨le_rfl,le_rfl⟩

theorem refine_sound {csm boxes cb cd bits degree} {l u : ℚ} {A : RationalInterval} {k : ℝ}
    (hb : boxes? csm cb cd = some boxes) (hk : k ∈ Icc (l:ℝ) u) (ha : A.Contains (weightReal csm k)) :
    (refine csm boxes l u bits degree A).Contains (weightReal csm k) := by
  unfold refine
  split
  · rename_i J hj
    have ht := (scaledWeight_sound hb hj).2.2 k hk
    exact ⟨by simpa only [Rat.cast_max] using max_le ha.1 ht.1,
      by simpa only [Rat.cast_min] using le_min ha.2 ht.2⟩
  · exact ha

/-- Standalone single-cell evaluator, always at least as tight as the natural extension. -/
def weight? (csm : Cosmology) (boxes : List RationalInterval) (l u : ℚ) (bits degree : Nat) :
    Option RationalInterval :=
  (SashimiEhSweep.weight? csm boxes l u bits degree).map (refine csm boxes l u bits degree)

theorem weight_sound {csm boxes cb cd bits degree l u P}
    (hb : boxes? csm cb cd = some boxes) (h : weight? csm boxes l u bits degree = some P) :
    0<l ∧ l≤u ∧ ∀ k ∈ Icc (l:ℝ) u, P.Contains (weightReal csm k) := by
  unfold weight? at h
  obtain ⟨A,hA,rfl⟩ := Option.map_eq_some_iff.mp h
  obtain ⟨hl,hu,hv⟩ := SashimiEhSweep.weight_sound hb hA
  exact ⟨hl,hu,fun k hk => refine_sound hb hk (hv k hk)⟩

/-- Reuse a checked mesh by tightening each cell, with the same endpoints. -/
def refineBins (csm : Cosmology) (boxes : List RationalInterval) (bits degree : Nat) (a : ℚ) :
    List Bin → List Bin
  | [] => []
  | q::qs => ⟨q.right, refine csm boxes a q.right bits degree q.power⟩ ::
      refineBins csm boxes bits degree q.right qs

theorem refineBins_good {csm boxes cb cd bits degree a bins}
    (hb : boxes? csm cb cd = some boxes) (hg : Good (weightReal csm) a bins) :
    Good (weightReal csm) a (refineBins csm boxes bits degree a bins) := by
  induction bins generalizing a with
  | nil => trivial
  | cons q qs ih =>
    exact ⟨hg.1,hg.2.1,fun k hk => refine_sound hb hk (hg.2.2.1 k hk),ih hg.2.2.2⟩

/-- A standalone refined mesh. Its untrusted candidate cells can be replayed
in short blocks using the same shape as the natural extension's checker. -/
def build? (csm : Cosmology) (boxes : List RationalInterval) (bits degree : Nat) (a : ℚ)
    (ends : List ℚ) : Option (List Bin) :=
  (SashimiEhSweep.build? csm boxes bits degree a ends).map (refineBins csm boxes bits degree a)

theorem build_sound {csm boxes cb cd bits degree a ends bins}
    (hb : boxes? csm cb cd = some boxes) (h : build? csm boxes bits degree a ends = some bins) :
    Good (weightReal csm) a bins := by
  unfold build? at h
  obtain ⟨original, ho, rfl⟩ := Option.map_eq_some_iff.mp h
  exact refineBins_good hb (SashimiEhSweep.build_sound hb ho)

#print axioms weight_sound
#print axioms refineBins_good
end SashimiEhTighterPower
