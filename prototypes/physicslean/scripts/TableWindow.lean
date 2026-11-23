import SashimiTableWindow

/-! Runtime for a whole-band spherical top-hat moment. Each rational input is
a numerator and positive denominator; Python only parses tabulated decimals. -/
open FloatLib.Numerics SashimiTable SashimiTableWindow

private def q (n d : String) : IO ℚ := do
  let some n := n.toInt? | throw (IO.userError "invalid numerator")
  let some d := d.toNat? | throw (IO.userError "invalid denominator")
  if d = 0 then throw (IO.userError "zero denominator")
  return mkRat n d

def main (args : List String) : IO Unit := do
  let chunkSize := (args[17]?).bind String.toNat? |>.getD 0
  let [model, mode, mn, md, rln, rld, run, rud, an, ad, bn, bd, tn, td, depth, bits, degree] := args.take 17
    | throw (IO.userError "expected table_window.py request")
  let p ← q mn md
  let radius := RationalInterval.mk (← q rln rld) (← q run rud)
  let a ← q an ad
  let b ← q bn bd
  let tol ← q tn td
  let some depth := depth.toNat? | throw (IO.userError "invalid depth")
  let some bits := bits.toNat? | throw (IO.userError "invalid precision")
  let some degree := degree.toNat? | throw (IO.userError "invalid degree")
  let mode ← match mode with
    | "log-log" => pure Interpolation.logLog
    | "linear" => pure Interpolation.linear
    | _ => throw (IO.userError "invalid interpolation")
  let ratio : ℚ → Option RationalInterval ← match model with
    | "cdm" => pure fun _ => some (.point 1)
    | "fdm" => pure fun k => SashimiSpectrum.fdm? p k bits degree
    | _ => throw (IO.userError "unknown window transfer")
  let stdin ← IO.getStdin
  let mut rows : Array Knot := #[]
  repeat
    let line ← stdin.getLine
    if line.isEmpty then break
    let [kn, kd, pn, pd] := line.trimAscii.toString.splitOn " "
      | throw (IO.userError "invalid knot")
    rows := rows.push ⟨← q kn kd, ← q pn pd⟩
  match SashimiTableWindow.integrate? 0 mode (weighted? ratio) rows.toList radius a b tol depth bits degree with
  | none => IO.println "unresolved"
  | some I =>
    if chunkSize > 0 then
      let all := (cells rows.toList).toArray
      let mut start := 0
      while start < all.size do
        let portion := (all.extract start (start+chunkSize)).toList
        let some part := sumCells? 0 mode (weighted? ratio) radius a b depth bits degree portion
          | throw (IO.userError "inconsistent chunk computation")
        IO.println s!"part {start} {part.lo} {part.hi}"
        start := start + chunkSize
    IO.println s!"enclosed {I.lo} {I.hi}"
