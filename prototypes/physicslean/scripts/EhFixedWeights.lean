import SashimiEhFixedWeight
import SashimiTableRequest

open FloatLib.Numerics SashimiEisensteinIntegral SashimiEhFixedCells

/-- Worker for a consecutive block of prescribed cells. Input gives both exact edges. -/
def main (args : List String) : IO Unit := do
  let [p,d] := args | throw (IO.userError "expected bits and degree")
  let some bits := p.toNat? | throw (IO.userError "invalid precision")
  let some degree := d.toNat? | throw (IO.userError "invalid degree")
  let some C := boxes? matched bits degree | throw (IO.userError "unresolved coefficients")
  let input ← IO.getStdin
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    let [l,u] := line.trimAscii.toString.splitOn " " | throw (IO.userError "expected cell edges")
    let .ok lo := SashimiTableRequest.rational l | throw (IO.userError "invalid lower edge")
    let .ok hi := SashimiTableRequest.rational u | throw (IO.userError "invalid upper edge")
    let some q := buildCell? matched C lo hi bits degree | throw (IO.userError s!"unresolved on {l}..{u}")
    IO.println s!"cell {lo} {q.right} {q.pointwise.lo} {q.pointwise.hi} {q.integral.lo} {q.integral.hi}"
