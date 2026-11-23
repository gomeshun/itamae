import HzExample
import HzExceptional
import HzFriedmann
import HzInterval
import HzNative
import HzRoundoffExamples
import Lean.Elab.Command
import Lean.Util.CollectAxioms

/-! Build-time check of the transitive axioms of the prototype's exported declarations. -/

open Lean in
run_cmd
  let env ← getEnv
  let mut checked : Nat := 0
  for (name, _) in env.constants do
    let some idx := env.getModuleIdxFor? name | continue
    unless env.header.moduleNames[idx.toNat]!.toString.startsWith "Hz" do continue
    for ax in ← collectAxioms name do
      unless #[``propext, ``Classical.choice, ``Quot.sound].contains ax do
        throwError "{name} depends on the unexpected axiom {ax}"
    checked := checked + 1
  if checked == 0 then throwError "no PhysicsLean declarations were audited"
  logInfo m!"PhysicsLean axiom audit passed on {checked} declarations"
