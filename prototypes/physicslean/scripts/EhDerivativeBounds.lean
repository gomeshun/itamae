import SashimiEhDerivativeBounds
import SashimiTableRequest

open FloatLib.Numerics SashimiEisensteinIntegral SashimiEhSweep

private def q (s : String) : IO ℚ := do
  let .ok x := SashimiTableRequest.rational s | throw (IO.userError "invalid fraction")
  return x

/-- Read the shared certified transfer ranges and then an exact mass grid.
Each row is either "cell k_right lo hi" or "mass M". -/
def main (args : List String) : IO Unit := do
  let [aa, prec, deg, tol] := args.take 4 | throw (IO.userError "expected lower wavenumber, bits, degree, width")
  let blocks ← match args.drop 4 with
    | [] => pure 0
    | [s] => match s.toNat? with
      | some (n+1) => pure (n+1)
      | _ => throw (IO.userError "invalid block size")
    | _ => throw (IO.userError "unexpected arguments")
  let a ← q aa
  let tol ← q tol
  let some bits := prec.toNat? | throw (IO.userError "invalid bits")
  let some degree := deg.toNat? | throw (IO.userError "invalid degree")
  let src ← IO.getStdin
  let mut bins := #[]
  let mut masses := #[]
  let mut left := a
  repeat
    let line ← src.getLine
    if line.isEmpty then break
    match line.trimAscii.toString.splitOn " " with
    | ["cell",kr,lower,upper] =>
      let r ← q kr
      let lP ← q lower
      let uP ← q upper
      if r ≤ left || lP > uP then throw (IO.userError "invalid shared cell or interval")
      bins := bins.push (⟨r,⟨lP,uP⟩⟩ : Bin)
      left := r
    | ["mass",mm] => masses := masses.push (← q mm)
    | _ => throw (IO.userError "invalid input")
  let cells := bins.toList
  let A := total (.point 8) bits degree a cells
  IO.println s!"anchor {A.lo} {A.hi}"
  for m in masses do
    if blocks > 0 then
      let some R := radiusBox? matched m bits degree | throw (IO.userError "invalid radius")
      IO.println s!"radius {m} {R.lo} {R.hi}"
      let mut lo := a
      for i in [:((bins.size + blocks - 1) / blocks)] do
        let part := bins.extract (i * blocks) ((i + 1) * blocks)
        let v := SashimiEhDerivativeBounds.sum R bits degree lo part.toList
        IO.println s!"part {m} {i} {v.lo} {v.hi}"
        if let some q := part.back? then lo := q.right
    if let some I := SashimiEhDerivativeBounds.derivative? m tol a cells A bits degree then
      IO.println s!"enclosed {m} {I.lo} {I.hi}"
    else
      IO.println s!"unresolved {m}"
