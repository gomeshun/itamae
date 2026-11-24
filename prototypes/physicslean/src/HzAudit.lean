import Sashimi
import HzExample
import HzExceptional
import HzFriedmann
import HzInterval
import HzNative
import HzRoundoffExamples
import Lean.Elab.Command
import Lean.Util.CollectAxioms

/-! Build-time check of the transitive axioms of the prototype's exported declarations. -/

namespace PhysicsLeanAudit
open Lean Elab Command

private def checkAllowed (name : Name) : CommandElabM Unit := do
  for ax in ← collectAxioms name do
    unless #[``propext, ``Classical.choice, ``Quot.sound].contains ax do
      throwError "{name} depends on the unexpected axiom {ax}"

/-- Apply the library's whitelist to theorems and axioms in the current proof
file, including private intermediate witnesses. A printout alone cannot fail
a build, and the final general theorem need not depend on every local witness. -/
def currentWitnesses : CommandElabM Unit := do
  let env ← getEnv
  let mut checked : Nat := 0
  for (name, info) in env.constants do
    if (env.getModuleIdxFor? name).isSome then continue
    match info with
    | .thmInfo _ | .axiomInfo _ =>
      checkAllowed name
      checked := checked + 1
    | _ => pure ()
  if checked == 0 then throwError "no local proof declarations were audited"
  logInfo m!"Concrete axiom audit passed on {checked} local proof declarations"

run_cmd
  let env ← getEnv
  let mut checked : Nat := 0
  for (name, _) in env.constants do
    let some idx := env.getModuleIdxFor? name | continue
    let mod := env.header.moduleNames[idx.toNat]!.toString
    unless ["Hz", "Sashimi", "Certified", "EhCentered"].any (fun p => mod.startsWith p) do continue
    checkAllowed name
    checked := checked + 1
  if checked == 0 then throwError "no PhysicsLean declarations were audited"
  logInfo m!"PhysicsLean axiom audit passed on {checked} declarations"

end PhysicsLeanAudit
