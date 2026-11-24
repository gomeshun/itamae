import CertifiedExprDerivative
import SashimiEhSweep

/-! Certified centered forms on a positive k cell. Midpoint values and derivative
bounds use one fixed set of nuisance parameters for each evaluation; this
also accommodates the high-frequency Silk envelope. -/
namespace EhCenteredWeight

open FloatLib.Numerics FloatLib.Numerics.Interval CertifiedExprDerivative Set
open SashimiEisensteinHu SashimiEisensteinIntegral SashimiEhSweep

private def expression (c : Cosmology) (transfer : Expr) : Expr :=
  -- Balance the growing tilt against T before squaring. A grid-rounded T² alone
  -- would lose relative precision at very large wavenumbers.
  .binary .mul (.binary .mul (tiltExpr c) transfer) transfer

def center? (c : Cosmology) (boxes : List RationalInterval) (l u : ℚ) (bits degree : Nat) :
    Option (RationalInterval × RationalInterval) := do
  let (t, args) ← fastForm? boxes ⟨l,u⟩ bits degree
  exprCenter? (expression c t) args l u bits degree

/-- An affine inclusion, equally applicable to a point range or a cell integral.
For the fast Silk form, the existentially chosen midpoint value may depend on k. -/
theorem center_sound {c boxes cb cd l u bits degree CD} (hb : boxes? c cb cd = some boxes)
    (h : center? c boxes l u bits degree = some CD) (k : ℝ) (hk : k ∈ Icc (l:ℝ) u) :
    ∃ v : ℝ, CD.1.Contains v ∧ ∃ d : ℝ, CD.2.Contains d ∧
      weightReal c k = v + d*(k-(((l+u)/2:ℚ):ℝ)) := by
  obtain ⟨⟨e,as⟩,hform,he⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨env,henv,hs,hval⟩ := fastForm_sound (a := coefficients c) (K := ⟨l,u⟩) hk
    (fun _ _ hi => coefficientBoxes_sound hb hi) hform
  have ht : ∀ (i : Nat) (I : RationalInterval), as[i]? = some I → I.Contains (env (i+1)) := by
    intro i I hi; exact hs (i+1) I (by simpa using hi)
  obtain ⟨v,hv,d,hd,H⟩ := exprCenter_sound ht he k hk
  refine ⟨v,hv,d,hd,?_⟩
  have hh : vary env k = env := by funext i; cases i <;> simp [vary,henv]
  change Real.exp ((((2+c.ns:ℚ):ℝ))*Real.log k) * (transferReal (coefficients c) k)^2 = _
  rw [hh] at H
  change (Real.exp ((((2+c.ns:ℚ):ℝ))*Real.log (env 0)) * e.eval env) * e.eval env = _ at H
  simpa only [henv,hval,pow_two,mul_assoc] using H

/-- Turn an affine inclusion into a pointwise range on the entire k cell. -/
def weight? (c : Cosmology) (boxes : List RationalInterval) (l u : ℚ) (bits degree : Nat) :
    Option RationalInterval := do
  guard (0 < l ∧ l ≤ u)
  let (C,D) ← center? c boxes l u bits degree
  return C.add (D.mul ⟨(l-u)/2,(u-l)/2⟩)

theorem weight_sound {c boxes cb cd l u bits degree P} (hb : boxes? c cb cd = some boxes)
    (h : weight? c boxes l u bits degree = some P) :
    0 < l ∧ l ≤ u ∧ ∀ k ∈ Icc (l:ℝ) u, P.Contains (weightReal c k) := by
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  have hab : 0 < l ∧ l ≤ u := by by_contra hn; simp [guard,hn] at hg
  obtain ⟨⟨C,D⟩,hcd,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  refine ⟨hab.1,hab.2,?_⟩
  intro k hk
  obtain ⟨v,hv,d,hd,he⟩ := center_sound hb hcd k hk
  rw [he]
  apply RationalInterval.contains_add hv (RationalInterval.contains_mul hd ?_)
  change (((l-u)/2:ℚ):ℝ) ≤ _ ∧ _ ≤ (((u-l)/2:ℚ):ℝ)
  push_cast
  constructor <;> linarith [hk.1,hk.2]

#print axioms center_sound
#print axioms weight_sound
end EhCenteredWeight
