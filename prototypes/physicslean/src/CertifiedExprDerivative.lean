import CertifiedDyadicInterval
import FloatLib.Numerics.Enclosure.Expression.ClosedProof
import FloatLib.Numerics.Enclosure.Rational.Converges
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Deriv
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.LinearCombination

/-! Forward interval differentiation in the first real variable of a FloatLib expression.
Other inputs may be known only through boxes. Each derivative describes the first
variable with all those other real inputs held fixed. -/
namespace CertifiedExprDerivative

open FloatLib.Numerics FloatLib.Numerics.Interval Set

structure Jet where
  value : RationalInterval
  slope : RationalInterval
deriving Repr

/-- Snap outward after each expression node to cap rational denominator growth.
Guard bits retain precision in small intermediate values. -/
private def trim (bits : Nat) (a : Jet) : Jet :=
  ⟨CertifiedDyadicInterval.rounded (bits+128) a.value, CertifiedDyadicInterval.rounded (bits+128) a.slope⟩

private def elem (u : UnaryOp) (I : RationalInterval) (bits degree : Nat) :=
  (Expr.unary u (.var 0)).gridEnclosureWith? [I] bits degree

private theorem elem_sound {u I bits degree O x} (hv : I.Contains x)
    (h : elem u I bits degree = some O) : O.Contains (u.eval x) := by
  have he := Expr.contains_gridEnclosureWith? (values := fun _ => x) (e := .unary u (.var 0))
    (fun i J hi => by cases i with
                      | zero => simp at hi; subst J; exact hv
                      | succ n => simp at hi) h
  exact he

private def unary? (op : UnaryOp) (a : Jet) (bits degree : Nat) : Option Jet := do
  match op with
  | .neg => return ⟨a.value.neg, a.slope.neg⟩
  | .exp =>
    let v ← elem .exp a.value bits degree
    return ⟨v, v.mul a.slope⟩
  | .log =>
    guard (0 < a.value.lo)
    let v ← elem .log a.value bits degree
    return ⟨v, a.slope.divPositive a.value⟩
  | .sin =>
    let v ← elem .sin a.value bits degree
    let c ← elem .cos a.value bits degree
    return ⟨v, c.mul a.slope⟩
  | .cos =>
    let v ← elem .cos a.value bits degree
    let c ← elem .sin a.value bits degree
    return ⟨v, c.neg.mul a.slope⟩
  | _ => none

private def binary? (op : BinaryOp) (a b : Jet) : Option Jet := do
  match op with
  | .add => return ⟨a.value.add b.value, a.slope.add b.slope⟩
  | .sub => return ⟨a.value.sub b.value, a.slope.sub b.slope⟩
  | .mul => return ⟨a.value.mul b.value, (a.slope.mul b.value).add (a.value.mul b.slope)⟩
  | .div =>
    guard (0 < b.value.lo)
    let v := a.value.divPositive b.value
    return ⟨v, (a.slope.sub (v.mul b.slope)).divPositive b.value⟩
  | _ => none

/-- A bounded forward pass for arithmetic, exp/log and sin/cos. Division
and log require a strictly positive boxed argument. Unsupported operations return none. -/
def jet? (e : Expr) (boxes : List RationalInterval) (bits degree : Nat) : Option Jet :=
  match e with
  | .const q => some ⟨.point q, .point 0⟩
  | .var i => do return ⟨← boxes[i]?, .point (if i = 0 then 1 else 0)⟩
  | .unary op e => do return trim bits (← unary? op (← jet? e boxes bits degree) bits degree)
  | .binary op e f => do return trim bits (← binary? op (← jet? e boxes bits degree) (← jet? f boxes bits degree))
  | _ => none
termination_by structural e

private def Models (j : Jet) (f : ℝ → ℝ) (x : ℝ) : Prop :=
  j.value.Contains (f x) ∧ ∃ d, j.slope.Contains d ∧ HasDerivAt f d x

private theorem model_trim {j f x} (h : Models j f x) (bits : Nat) : Models (trim bits j) f x :=
  ⟨CertifiedDyadicInterval.round_sound _ h.1, by
    obtain ⟨d,hd,hf⟩ := h.2; exact ⟨d,CertifiedDyadicInterval.round_sound _ hd,hf⟩⟩

private theorem model_unary {j k f x op bits degree} (hf : Models j f x)
    (h : unary? op j bits degree = some k) : Models k (fun t => op.eval (f t)) x := by
  obtain ⟨hv,d,hd,df⟩ := hf
  cases op <;> simp only [unary?] at h <;> try contradiction
  case neg =>
    cases h
    exact ⟨RationalInterval.contains_neg hv,-d,RationalInterval.contains_neg hd,df.neg⟩
  case exp =>
    obtain ⟨v,he,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    have hv' := elem_sound hv he
    exact ⟨hv',_,RationalInterval.contains_mul hv' hd,df.exp⟩
  case log =>
    obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
    have hl : 0 < j.value.lo := by by_contra hn; simp [guard,hn] at hg
    obtain ⟨v,he,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    have hpos : 0 < f x := lt_of_lt_of_le (by exact_mod_cast hl : (0:ℝ)<j.value.lo) hv.1
    exact ⟨elem_sound hv he,_,RationalInterval.contains_divPositive hd hv hl,df.log hpos.ne'⟩
  case sin =>
    obtain ⟨v,he,h⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨w,hc,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    exact ⟨elem_sound hv he,_,RationalInterval.contains_mul (elem_sound hv hc) hd,df.sin⟩
  case cos =>
    obtain ⟨v,he,h⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨w,hc,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    exact ⟨elem_sound hv he, _, RationalInterval.contains_mul
      (RationalInterval.contains_neg (elem_sound hv hc)) hd, df.cos⟩

private theorem model_binary {i j k f g x op} (hf : Models i f x) (hg : Models j g x)
    (h : binary? op i j = some k) : Models k (fun t => op.eval (f t) (g t)) x := by
  obtain ⟨hv,d,hd,df⟩ := hf
  obtain ⟨hw,e,he,dg⟩ := hg
  cases op <;> simp only [binary?] at h <;> try contradiction
  case add =>
    cases h; exact ⟨RationalInterval.contains_add hv hw,_,RationalInterval.contains_add hd he,df.add dg⟩
  case sub =>
    cases h; exact ⟨RationalInterval.contains_sub hv hw,_,RationalInterval.contains_sub hd he,df.sub dg⟩
  case mul =>
    cases h
    exact ⟨RationalInterval.contains_mul hv hw, _, RationalInterval.contains_add
      (RationalInterval.contains_mul hd hw) (RationalInterval.contains_mul hv he), df.mul dg⟩
  case div =>
    obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
    have hl : 0 < j.value.lo := by by_contra hn; simp [guard,hn] at hg
    cases h
    have hne : g x ≠ 0 := (lt_of_lt_of_le (by exact_mod_cast hl : (0:ℝ)<j.value.lo) hw.1).ne'
    have hquot := RationalInterval.contains_divPositive hv hw hl
    refine ⟨hquot,_,RationalInterval.contains_divPositive
      (RationalInterval.contains_sub hd (RationalInterval.contains_mul hquot he)) hw hl,?_⟩
    change HasDerivAt (fun t => f t / g t) ((d - (f x / g x)*e) / g x) x
    rw [show (d - (f x / g x)*e) / g x = (d*g x - f x*e)/(g x)^2 by field_simp]
    exact df.div dg hne

/-- Replace the differentiable variable, leaving all other arguments fixed. -/
def vary (values : Nat → ℝ) (x : ℝ) (i : Nat) : ℝ := if i = 0 then x else values i

theorem jet_sound (e : Expr) {boxes : List RationalInterval} {bits degree : Nat} {J : Jet}
    {values : Nat → ℝ} {x : ℝ}
    (hv : ∀ (i : Nat) (I : RationalInterval), boxes[i]? = some I → I.Contains (vary values x i))
    (h : jet? e boxes bits degree = some J) :
    J.value.Contains (e.eval (vary values x)) ∧
    ∃ d, J.slope.Contains d ∧ HasDerivAt (fun t => e.eval (vary values t)) d x := by
  change Models J (fun t => e.eval (vary values t)) x
  cases e with
  | const q =>
    simp [jet?] at h; subst J
    exact ⟨RationalInterval.contains_point q,0,by norm_num [RationalInterval.point,RationalInterval.Contains],
      hasDerivAt_const x _⟩
  | var i =>
    unfold jet? at h
    obtain ⟨I,hI,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    refine ⟨hv i I hI,?_⟩
    by_cases h0 : i=0
    · subst i
      exact ⟨1, by norm_num [RationalInterval.point, RationalInterval.Contains],
        by simpa [vary, Expr.eval] using hasDerivAt_id' x⟩
    · exact ⟨0, by norm_num [h0, RationalInterval.point, RationalInterval.Contains],
        by simpa [vary, h0, Expr.eval] using hasDerivAt_const x (values i)⟩
  | unary op e =>
    unfold jet? at h
    obtain ⟨a,ha,h⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨b,hb,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    exact model_trim (model_unary (j := a) (f := fun t => e.eval (vary values t)) (x := x) (jet_sound e hv ha) hb) _
  | binary op e f =>
    unfold jet? at h
    obtain ⟨a,ha,h⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨b,hb,h⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨k,hk,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    exact model_trim (model_binary (i := a) (j := b) (x := x)
      (f := fun t => e.eval (vary values t)) (g := fun t => f.eval (vary values t))
      (jet_sound e hv ha) (jet_sound f hv hb) hk) _
  | ternary => simp [jet?] at h
  | call => simp [jet?] at h
termination_by structural e

private theorem chord {l u m k : ℝ} {f : ℝ → ℝ} {D : RationalInterval}
    (h : ∀ x ∈ Icc l u, ∃ d, D.Contains d ∧ HasDerivAt f d x)
    (hm : m ∈ Icc l u) (hk : k ∈ Icc l u) :
    ∃ d, D.Contains d ∧ f k = f m + d*(k-m) := by
  have step (a b : ℝ) (ha : a ∈ Icc l u) (hb : b ∈ Icc l u) (hab : a < b) :
      ∃ d, D.Contains d ∧ f b - f a = d*(b-a) := by
    have sub (x : ℝ) (hx : x ∈ Icc a b) : x ∈ Icc l u := ⟨ha.1.trans hx.1,hx.2.trans hb.2⟩
    have diff (x : ℝ) (hx : x ∈ Icc a b) : DifferentiableAt ℝ f x := by
      obtain ⟨_,_,hd⟩ := h x (sub x hx); exact hd.differentiableAt
    obtain ⟨y,hy,hd⟩ := exists_deriv_eq_slope f hab
      (fun x hx => (diff x hx).continuousAt.continuousWithinAt)
      (fun x hx => (diff x ⟨hx.1.le,hx.2.le⟩).differentiableWithinAt)
    obtain ⟨d,hD,hd'⟩ := h y (sub y ⟨hy.1.le,hy.2.le⟩)
    rw [hd'.deriv] at hd
    exact ⟨d,hD,((eq_div_iff (sub_ne_zero.mpr hab.ne')).mp hd).symm⟩
  rcases lt_trichotomy m k with hmk | heq | hkm
  · obtain ⟨d,hd,hs⟩ := step m k hm hk hmk; exact ⟨d,hd,by linarith⟩
  · subst k; obtain ⟨d,hd,_⟩ := h m hm; exact ⟨d,hd,by ring⟩
  · obtain ⟨d,hd,hs⟩ := step k m hk hm hkm; exact ⟨d,hd,by linear_combination -hs⟩

/-- FloatLib at the exact midpoint, with an interval jet over the complete cell. -/
def exprCenter? (e : Expr) (args : List RationalInterval) (l u : ℚ) (bits degree : Nat) :
    Option (RationalInterval × RationalInterval) := do
  let c ← e.gridEnclosureWith? (.point ((l+u)/2)::args) bits degree
  let j ← jet? e (⟨l,u⟩::args) bits degree
  return (c,j.slope)

theorem exprCenter_sound {e args l u bits degree CD} {values : Nat → ℝ}
    (ha : ∀ (i : Nat) (I : RationalInterval), args[i]? = some I → I.Contains (values (i+1)))
    (h : exprCenter? e args l u bits degree = some CD) (k : ℝ) (hk : k ∈ Icc (l:ℝ) u) :
    ∃ v : ℝ, CD.1.Contains v ∧ ∃ d : ℝ, CD.2.Contains d ∧
      e.eval (vary values k) = v + d*(k - (((l+u)/2 : ℚ):ℝ)) := by
  obtain ⟨C,hC,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨J,hJ,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  let m := ((l+u)/2 : ℚ)
  have inputs {I : RationalInterval} {t : ℝ} (ht : I.Contains t) :
      ∀ (i : Nat) (B : RationalInterval), (I::args)[i]? = some B → B.Contains (vary values t i) := by
    intro i B hb
    cases i with
    | zero => simp at hb; subst B; exact ht
    | succ i => exact ha i B (by simpa using hb)
  have hc := Expr.contains_gridEnclosureWith? (values := vary values (m:ℝ))
    (inputs (RationalInterval.contains_point m)) hC
  have hdiff (t : ℝ) (ht : t ∈ Icc (l:ℝ) u) := (jet_sound e
    (inputs (I := ⟨l,u⟩) ht) hJ).2
  have hm : (m:ℝ) ∈ Icc (l:ℝ) u := by
    have hle := hk.1.trans hk.2
    change (l:ℝ) ≤ (m:ℝ) ∧ (m:ℝ) ≤ u
    dsimp [m]; push_cast; constructor <;> linarith
  exact ⟨_,hc,chord hdiff hm hk⟩

#print axioms jet_sound
#print axioms exprCenter_sound
end CertifiedExprDerivative
