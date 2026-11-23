import SashimiMassVariance
import SashimiGrowthVariance

open FloatLib.Numerics FloatLib.Numerics.Interval SashimiTable SashimiMassVariance

private def rat (s : String) : IO ℚ := do
  match s.splitOn "/" with
  | [n,d] =>
    let some n := n.toInt? | throw (IO.userError "invalid numerator")
    let some d := d.toNat? | throw (IO.userError "invalid denominator")
    if d = 0 then throw (IO.userError "zero denominator")
    return mkRat n d
  | [n] =>
    let some n := n.toInt? | throw (IO.userError "invalid rational")
    return n
  | _ => throw (IO.userError "invalid rational")

private def natural (s : String) : IO Nat := do
  let some n := s.toNat? | throw (IO.userError "invalid natural number")
  return n

def main (args : List String) : IO Unit := do
  let [op,coords,loc,model,param,omega,hh,scale,rho,mass,growth,a,b,tol,depth,bits,degree] := args.take 17
    | throw (IO.userError "expected arguments from mass_variance.py")
  let mode ← match coords with
    | "linear" => pure Interpolation.linear
    | "log-log" => pure Interpolation.logLog
    | _ => throw (IO.userError "unknown interpolant")
  let placement ← match loc with
    | "at-knots" => pure Placement.atKnots
    | "continuous" => pure Placement.continuous
    | _ => throw (IO.userError "unknown placement")
  let p ← rat param
  let omega ← rat omega
  let h ← rat hh
  let ax : Expr := if omega > 0 && h > 0 then
    .binary .div (SashimiSpectrum.wdmAlphaExpr (.const p) (.const omega) (.const h)) (.const h)
    else .const p
  let ratio : Expr → Expr ← match model with
    | "cdm" => pure (fun _ => .const 1)
    | "fdm" => pure (SashimiSpectrum.fdmPowerExpr (.const p))
    | "wdm-q5" => pure (SashimiSpectrum.wdmPowerExpr .publishedQ5 ax)
    | "wdm-q10" => pure (SashimiSpectrum.wdmPowerExpr .standardT2Q10 ax)
    | _ => throw (IO.userError "unknown transfer")
  let cal ← match scale with
    | "sashimi-w" | "w-family-legacy" => pure Calibration.wFamilyLegacy
    | "w-main" => pure Calibration.wMain
    | _ => do
      let c ← rat scale
      if c ≤ 0 then throw (IO.userError "positive mass calibration required")
      pure (Calibration.fixed c)
  let rho ← rat rho
  let mass ← rat mass
  let growth ← rat growth
  let a ← rat a
  let b ← rat b
  let tol ← rat tol
  let depth ← natural depth
  let bits ← natural bits
  let degree ← natural degree
  let (D, rawTol, finish) : ℚ × ℚ × (RationalInterval → Option RationalInterval) ←
    match args.drop 17 with
    | [] => pure (growth, tol, some)
    | [source, output, m, l, z, growthTol, rawTol] => do
      let q ← match output with
        | "variance" => pure SashimiGrowthVariance.Quantity.variance
        | "sigma" => pure SashimiGrowthVariance.Quantity.sigma
        | "derivative" => pure SashimiGrowthVariance.Quantity.derivative
        | "legacy-f-linear" => pure SashimiGrowthVariance.Quantity.legacyFDerivative
        | _ => throw (IO.userError "unknown growth observable")
      let rawTol ← rat rawTol
      if source = "fixed" then
        pure (1, rawTol, fun raw => SashimiGrowthVariance.combine? q (.point growth) raw tol bits)
      else if source = "cpt" then
        let m ← rat m
        let l ← rat l
        let z ← rat z
        let gt ← rat growthTol
        pure (1, rawTol, fun raw => SashimiGrowthVariance.atRedshift? q m l z raw tol gt bits degree)
      else throw (IO.userError "unknown growth model")
    | _ => throw (IO.userError "invalid growth parameters")
  let input ← IO.getStdin
  let mut rows : Array Knot := #[]
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    let [k,power] := line.trimAscii.toString.splitOn " " | throw (IO.userError "invalid table row")
    rows := rows.push ⟨← rat k, ← rat power⟩
  let out ← if op = "sharp-derivative" then
    pure (sharpDerivative? mode placement ratio cal rho mass D rows.toList a b rawTol bits degree)
  else do
    let q ← match op with
      | "sharp" => pure Quantity.sharp
      | "top-hat" => pure Quantity.topHat
      | "top-hat-derivative" => pure Quantity.topHatDerivative
      | _ => throw (IO.userError "unknown mass operation")
    pure (evaluate? q mode placement ratio cal rho mass D rows.toList a b rawTol depth bits degree)
  match out.bind finish with
  | some I => IO.println s!"enclosed {I.lo} {I.hi}"
  | none => IO.println "unresolved"
