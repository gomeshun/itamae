import SashimiEpsKernel
import SashimiHalo
import SashimiTableRequest

/-! Feed any selected variance and growth calculation through the *same* guarded
SASHIMI-C downstream formulas. Inputs to this driver are exact rational values
or boxes. Error in an upstream numerical quadrature is not included in a point
input. The halo history here specializes the public assembly API to zi = 0. -/
open FloatLib.Numerics FloatLib.Numerics.Interval

private abbrev I := RationalInterval

private def solve (tag : String) (args : List ℚ) : Option (List I) := do
  match tag, args with
  | "F", [b1lo,b1hi,b0lo,b0hi,s1,s0] =>
    let f ← EpsKernel.enclose? ⟨b1lo,b1hi⟩ ⟨b0lo,b0hi⟩ (.point s1) (.point s0) 96 40
    return [f]
  | "host", [z,dDlo,dDhi,sq,s0] =>
    guard (0 ≤ z ∧ dDlo ≤ dDhi ∧ s0 < sq)
    let ds := Expr.const (sq-s0)
    let slope := Expr.var 0
    let m := SashimiHalo.massHistory (.const (10^12)) (.const z) (.const 1) slope ds
    let matter := (0.315 : ℚ)*(1+z)^3
    let Δ := SashimiHalo.virialOverdensity (.const (matter/(0.685+matter)))
    let vir := SashimiHalo.virialMassFit m (SashimiHalo.correa200 m z) Δ
    let rate := SashimiHalo.massRate vir (.const z) (.const 1) slope ds
    let mbox ← m.gridEnclosureWith? [⟨dDlo,dDhi⟩] 96 32
    let ratebox ← rate.gridEnclosureWith? [⟨dDlo,dDhi⟩] 96 32
    return [mbox, ratebox]
  | _, _ => none

def main : IO Unit := do
  let input ← IO.getStdin
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    let tag :: args := line.trimAscii.toString.splitOn " " | throw (IO.userError "missing calculation")
    let .ok xs := args.mapM SashimiTableRequest.rational | throw (IO.userError "invalid rational input")
    let some boxes := solve tag xs | throw (IO.userError s!"unresolved {tag} for {args}")
    IO.println (String.intercalate " " (boxes.flatMap fun b => [toString b.lo, toString b.hi]))
