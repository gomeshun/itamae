import CertifiedIntegralRatio
import SashimiTableRequest
import SashimiTighterWindow
open FloatLib.Numerics SashimiEhSweep CertifiedIntegralRatio

private def rational (s : String) : IO ℚ := do
  let .ok v := SashimiTableRequest.rational s | throw (IO.userError "invalid rational")
  return v

private def units (I : RationalInterval) (bits : Nat) : String :=
  let scale : ℚ := 2^bits
  s!"{(I.lo*scale).num} {(I.hi*scale).num}"

/-- Emit exact dyadic cell bounds for the independent rational sign checker. -/
def main (args : List String) : IO Unit := do
  let [b,d] := args | throw (IO.userError "expected bits and degree")
  let some bits := b.toNat? | throw (IO.userError "expected bits")
  let some degree := d.toNat? | throw (IO.userError "expected degree")
  let input ← IO.getStdin
  let mut fixed := false
  let tighter := fun (r : RationalInterval) (l u : ℚ) =>
    rounded bits (SashimiTighterWindow.windowSq (l*r.lo) (u*r.hi) bits degree)
  let mut ranges : Array (ℚ × ℚ) := #[]
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    match line.trimAscii.toString.splitOn " " with
    | ["p", l, u, p, q] =>
      let lo ← rational l
      let hi ← rational u
      let P : RationalInterval := ⟨← rational p, ← rational q⟩
      let E := rounded bits (P.scaleNonnegative (hi-lo))
      let A := window (.point 8) lo hi bits degree
      IO.println s!"b {units E bits} {units A bits}"
      ranges := ranges.push (lo,hi)
    | ["f", l, u, al, au] =>
      fixed := true
      let lo ← rational l
      let hi ← rational u
      let E := rounded bits ⟨← rational al, ← rational au⟩
      let A := tighter (.point 8) lo hi
      IO.println s!"b {units E bits} {units A bits}"
      ranges := ranges.push (lo,hi)
    | ["r", id, l, u] =>
      let r : RationalInterval := ⟨← rational l, ← rational u⟩
      for (a,b) in ranges do
        let X := if fixed then tighter r a b else window r a b bits degree
        IO.println s!"w {id} {units X bits}"
    | _ => throw (IO.userError "expected spectral or radius row")
