import HzRoundoff

/-! Native binary64 evaluation. FloatLib proves the arithmetic and square-root
laws for Lean's logical model of `Float`. Execution on the host relies on its
compiler and processor implementing those laws. -/

namespace HzNative

open FloatLib.Floats FloatLib.Floats.Formats.BinaryInterchange

/-- Evaluate H(z) on stored inputs, multiplying twice for the cube. -/
def hubble (h₀ Ωm z : Float) : Float :=
  let a := 1 + z
  h₀ * Float.sqrt ((1 - Ωm) + Ωm * ((a * a) * a))

/-- Read the result through FloatLib's binary64 bridge. The generated theorem
requires finite inputs, valid intermediate ranges and the square-root domain.
Its Physlib reference uses the exact binary64 values of the stored inputs. -/
roundoff mixed theorem hubble_roundoff using HzRoundoff.real_reference :=
  fun (h₀ Ωm z : Float) => (ExecFloat.Binary.ofFloat (hubble h₀ Ωm z) : Binary64).toReal

/-- Exact output words of the logical native program, including the rounding
of the source-level decimal literals. -/
theorem reference_words :
    (hubble 70 0.3 1).toBits = 0x405ecfda9cbb1a15 ∧
    (hubble 67.4 0.315 2).toBits = 0x40698a578a7d1554 := by
  decide +kernel

#print axioms hubble_roundoff
#print axioms reference_words

end HzNative
