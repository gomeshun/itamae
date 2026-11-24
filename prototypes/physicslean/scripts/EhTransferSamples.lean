import SashimiEisensteinHu

/-! Evaluate the matched baryonic transfer on the exact binary64 bit patterns
used by the comparison library. Corresponding real proofs live separately. -/
open SashimiEisensteinHu

def main : IO Unit := do
  let input ← IO.getStdin
  let coefficients := prepareFloat 0.674 0.315 0.049 2.7255
  repeat
    let line ← input.getLine
    if line.isEmpty then break
    let some bits := line.trimAscii.toString.toNat? | throw (IO.userError "invalid IEEE bits")
    if bits ≥ 2^64 then throw (IO.userError "invalid binary64")
    let k := Float.ofBits (UInt64.ofNat bits)
    unless k.isFinite && k > 0 do throw (IO.userError "wavenumber must be positive and finite")
    IO.println (transferFloat coefficients k).toBits.toNat
