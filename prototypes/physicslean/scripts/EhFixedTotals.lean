import SashimiEhFixedCells
import SashimiTableRequest

open FloatLib.Numerics SashimiEisensteinIntegral SashimiEhFixedCells

private def rational (s : String) : IO ℚ := do
  let .ok x := SashimiTableRequest.rational s | throw (IO.userError "invalid rational input")
  return x

/-- Integrate both windows on the same immutable power cells for a group of masses. -/
def main (args : List String) : IO Unit := do
  let [p,d] := args | throw (IO.userError "expected bits and degree")
  let some bits := p.toNat? | throw (IO.userError "invalid precision")
  let some degree := d.toNat? | throw (IO.userError "invalid degree")
  let input ← IO.getStdin
  let a : ℚ := 1/100000000
  let mut edge := a
  let mut cells : Array Cell := #[]
  let mut masses := #[]
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    match line.trimAscii.toString.splitOn " " with
    | ["cell", l, u, p, q, al, au] =>
      let left ← rational l
      let right ← rational u
      let W := (⟨← rational p, ← rational q⟩ : RationalInterval)
      let A := (⟨← rational al, ← rational au⟩ : RationalInterval)
      if left != edge || left ≥ right || W.lo > W.hi || A.lo > A.hi then
        throw (IO.userError "unordered power mesh")
      cells := cells.push ⟨right,W,A⟩
      edge := right
    | ["mass", m] => masses := masses.push (← rational m)
    | _ => throw (IO.userError "expected cell or mass")
  if edge != 10000000000 then throw (IO.userError "incomplete power mesh")
  let cs := cells.toList
  let A := variance (.point 8) bits degree a cs
  IO.println s!"anchor {A.lo} {A.hi}"
  for m in masses do
    let some v := mass? m 1e20 a cs A bits degree | throw (IO.userError s!"unresolved mass {m}")
    IO.println s!"mass {m} {v.radius.lo} {v.radius.hi} {v.numerator.lo} {v.numerator.hi} {v.variance.lo} {v.variance.hi}"
    IO.println s!"slope {m} {v.derivative.lo} {v.derivative.hi}"
