import HzRoundoff

/-! Tight absolute errors include binary64 input conversion and every operation, including sqrt. -/

set_option Elab.async false
set_option exponentiation.threshold 4096
set_option maxHeartbeats 3000000

namespace HzRoundoffExamples

open HzRoundoff FloatLib.Floats

/-- Actual stored values for all inexact decimal inputs in the benchmark. -/
theorem stored_inputs :
    (ExecFloat.Binary.ofRat (3/10) : Binary64).toRat? = some (5404319552844595/18014398509481984) ∧
    (ExecFloat.Binary.ofRat (337/5) : Binary64).toRat? = some (2371426678787277/35184372088832) ∧
    (ExecFloat.Binary.ofRat (63/200) : Binary64).toRat? = some (5674535530486825/18014398509481984) := by
  float_decide

theorem redshift_one : Accurate 70 (3/10) 1 (2/10^14) := by
  certify_hz64 (8672787140647445/70368744177664)

theorem redshift_two : Accurate (337/5) (63/200) 2 (5/10^14) := by
  certify_hz64 (1797245752001877/8796093022208)

theorem lambda_only : Accurate 70 0 11 0 := by
  certify_hz64 70

theorem matter_only : Accurate 70 1 3 0 := by
  certify_hz64 560

theorem today : Accurate (337/5) (63/200) 0 (6/10^15) := by
  certify_hz64 (2371426678787277/35184372088832)

/-- Forward roundoff bounds for the same five inputs. Each result specializes
    the symbolically composed error budget and discharges its hypotheses. -/
theorem forward_redshift_one : Accurate 70 (3/10) 1 (1/10^12) := by
  certify_hz64_forward 70 (3/10) 1 (1745144855606067/562949953421312)

theorem forward_redshift_two : Accurate (337/5) (63/200) 2 (1/10^12) := by
  certify_hz64_forward (337/5) (63/200) 2 (2586755035970929/281474976710656)

theorem forward_lambda_only : Accurate 70 0 11 (1/10^12) := by
  certify_hz64_forward 70 0 11 1

theorem forward_matter_only : Accurate 70 1 3 (2/10^12) := by
  certify_hz64_forward 70 1 3 64

theorem forward_today : Accurate (337/5) (63/200) 0 (1/10^12) := by
  certify_hz64_forward (337/5) (63/200) 0 1

#print axioms forward_redshift_one
#print axioms forward_redshift_two
#print axioms forward_lambda_only
#print axioms forward_matter_only
#print axioms forward_today
#print axioms stored_inputs
#print axioms redshift_one
#print axioms redshift_two
#print axioms lambda_only
#print axioms matter_only
#print axioms today

end HzRoundoffExamples
