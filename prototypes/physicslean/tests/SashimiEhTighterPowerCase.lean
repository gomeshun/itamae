import HzAudit
import SashimiEhTighterPower

/-! Replay a real EH cell straddling k = 1 h/Mpc. The width inequality
compares two kernel-recomputed enclosures of the same continuum function. -/
namespace SashimiEhTighterPowerCase
open SashimiEisensteinIntegral SashimiEhSweep FloatLib.Numerics Set

private def left : ℚ := 376832/390625
private def right : ℚ := 393216/390625
private def original : RationalInterval :=
  ⟨345041510457782592745924945645853065853491860508688766620599913992271401160657707161536634849933 /
      17498005798264095394980017816940970922825355447145699491406164851279623993595007385788105416184430592,
   956788806059204154704230474753778318646467040158641629936260278648611616614912189642752786528775 /
      34996011596528190789960035633881941845650710894291398982812329702559247987190014771576210832368861184⟩
private def tighter : RationalInterval :=
  ⟨12220572868330742076805101588035549399165814080484657859184575229897293613309665621431395915475 /
      546812681195752981093125556779405341338292357723303109106442651602488249799843980805878294255763456,
   26381033485524482124970251021098104205137527231036540225289909236319274962128868365105294833745 /
      1093625362391505962186251113558810682676584715446606218212885303204976499599687961611756588511526912⟩

private def outcomes : Option (RationalInterval × RationalInterval) := do
  let boxes ← boxes? matched 112 16
  let before ← SashimiEhSweep.weight? matched boxes left right 112 16
  let after ← SashimiEhTighterPower.weight? matched boxes left right 112 16
  return (before, after)

set_option maxRecDepth 100000 in
set_option maxHeartbeats 16000000 in
private theorem checked : outcomes = some (original, tighter) := by decide +kernel

theorem contains {k : ℝ} (hk : k ∈ Icc (left : ℝ) right) : tighter.Contains (weightReal matched k) := by
  have ho := checked
  unfold outcomes at ho
  obtain ⟨boxes,hb,ho⟩ := Option.bind_eq_some_iff.mp ho
  obtain ⟨before,_,ho⟩ := Option.bind_eq_some_iff.mp ho
  obtain ⟨after,ha,ho⟩ := Option.bind_eq_some_iff.mp ho
  have he : (before,after) = (original,tighter) := Option.some.inj ho
  have eqAfter : after = tighter := congrArg Prod.snd he
  rw [← eqAfter]
  exact (SashimiEhTighterPower.weight_sound hb ha).2.2 k hk

theorem width_quarter : 4*(tighter.hi-tighter.lo) ≤ original.hi-original.lo := by decide +kernel

#print axioms contains
#print axioms width_quarter
end SashimiEhTighterPowerCase

run_cmd PhysicsLeanAudit.currentWitnesses
