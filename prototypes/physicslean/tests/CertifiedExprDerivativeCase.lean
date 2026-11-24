import CertifiedExprDerivative

namespace CertifiedExprDerivativeCase
open FloatLib.Numerics FloatLib.Numerics.Interval CertifiedExprDerivative Set

-- t²+p with t∈[1,2] and a fixed uncertain p∈[-1,1].
private def e : Expr := .binary .add (.binary .mul (.var 0) (.var 0)) (.var 1)
private def C : RationalInterval := ⟨5/4,13/4⟩
private def D : RationalInterval := ⟨2,4⟩

private theorem result : exprCenter? e [⟨-1,1⟩] 1 2 64 12 = some (C,D) := by decide +kernel

theorem centered (t p : ℝ) (ht : t ∈ Icc 1 2) (hp : p ∈ Icc (-1) 1) :
    ∃ v : ℝ, C.Contains v ∧ ∃ d : ℝ, D.Contains d ∧ t*t+p = v+d*(t-3/2) := by
  have hm : ∀ (i : Nat) (I : RationalInterval), [⟨-1,1⟩][i]? = some I → I.Contains p := by
    intro i I hi
    cases i with
    | zero => simp at hi; subst I; simpa [RationalInterval.Contains] using hp
    | succ i => simp at hi
  obtain ⟨v,hv,d,hd,he⟩ := exprCenter_sound (values := fun _ => p) hm result t (by simpa using ht)
  exact ⟨v,hv,d,hd,by norm_num [e, Expr.eval, BinaryOp.eval, vary] at he; exact he⟩

#print axioms centered
end CertifiedExprDerivativeCase
