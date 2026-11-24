import SashimiEhSweep
import SashimiTableRequest
open SashimiEhSweep SashimiEisensteinIntegral FloatLib.Numerics

private def showI (I : RationalInterval) := s!"{I.lo} {I.hi}"

/-- Stdin supplies one bounded physical mass per line, in exact solar masses. -/
def main (args : List String) : IO Unit := do
  let [sub,prec,deg] := args | throw (IO.userError "expected subdivisions, bits, degree")
  let some steps := sub.toNat? | throw (IO.userError "invalid subdivisions")
  let some bits := prec.toNat? | throw (IO.userError "invalid bits")
  let some degree := deg.toNat? | throw (IO.userError "invalid degree")
  if steps = 0 then throw (IO.userError "empty mesh")
  let a : ℚ := 1/100000000
  let b : ℚ := 10000000000
  let some C := boxes? matched bits degree | throw (IO.userError "coefficient box unresolved")
  let mut l := a
  let mut bins : Array Bin := #[]
  while l < b do
    let r := min (2*l) b
    for j in [0:steps] do
      let lo := l + (r-l)*(j:ℚ)/steps
      let hi := l + (r-l)*((j+1:ℕ):ℚ)/steps
      let some I := weight? matched C lo hi bits degree | throw (IO.userError s!"unresolved weighted range on {lo}..{hi}")
      bins := bins.push ⟨hi,I⟩
      IO.println s!"bin {lo} {hi} {showI I}"
    l := r
  let bs := bins.toList
  let A := total (.point 8) bits degree a bs
  IO.println s!"anchor {showI A}"
  let input ← IO.getStdin
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    let .ok m := SashimiTableRequest.rational line.trimAscii.toString | throw (IO.userError "invalid rational mass")
    let some R := radiusBox? matched m bits degree | throw (IO.userError "unresolved mass radius")
    if R.lo < 0 then throw (IO.userError "unresolved positive radius")
    let I := total R bits degree a bs
    let some Z := SashimiNormalization.apply? matched.sigma8 1e20 I A | throw (IO.userError "normalization unresolved")
    IO.println s!"mass {m} {showI R} {showI I} {showI Z}"
