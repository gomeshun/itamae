import FloatLib.Numerics.Enclosure.Rational.Proof
import Mathlib.Algebra.Order.Floor.Ring
import Mathlib.Tactic.Positivity

/-! Outward rounding of a rational interval onto a fixed dyadic grid. -/
namespace CertifiedDyadicInterval

open FloatLib.Numerics

/-- Round the lower endpoint down and the upper endpoint up to multiples of
`2^(-bits)`. The input endpoints may be arbitrary rational numbers. -/
def rounded (bits : Nat) (I : RationalInterval) : RationalInterval :=
  let s : ℚ := 2^bits
  ⟨⌊I.lo*s⌋/s,⌈I.hi*s⌉/s⟩

theorem round_sound (bits : Nat) {I : RationalInterval} {x : ℝ} (h : I.Contains x) :
    (rounded bits I).Contains x := by
  let s : ℚ := 2^bits
  have hs : 0 < s := by dsimp [s]; positivity
  have hl : (⌊I.lo*s⌋ : ℚ)/s ≤ I.lo := (div_le_iff₀ hs).2 (Int.floor_le _)
  have hu : I.hi ≤ (⌈I.hi*s⌉ : ℚ)/s := (le_div_iff₀ hs).2 (Int.le_ceil _)
  exact ⟨(by exact_mod_cast hl : (((rounded bits I).lo : ℚ) : ℝ) ≤ I.lo).trans h.1,
    h.2.trans (by exact_mod_cast hu : (I.hi : ℝ) ≤ ((rounded bits I).hi : ℝ))⟩

#print axioms round_sound
end CertifiedDyadicInterval
