import SashimiTableRequest

/-! Rational stdin/stdout transport for the pure table request interpreter. -/
open SashimiTable

private def checked {α} : Except String α → IO α
  | .ok x => pure x
  | .error message => throw (IO.userError message)

def main (args : List String) : IO Unit := do
  let input ← IO.getStdin
  let mut rows : Array Knot := #[]
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    let [k, power] := line.trimAscii.toString.splitOn " " | throw (IO.userError "malformed row")
    rows := rows.push ⟨← checked (SashimiTableRequest.rational k),
      ← checked (SashimiTableRequest.rational power)⟩
  match ← checked (SashimiTableRequest.run args rows.toList) with
  | some I => IO.println s!"enclosed {I.lo} {I.hi}"
  | none => IO.println "unresolved"
