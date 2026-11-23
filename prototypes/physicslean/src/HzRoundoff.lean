import HzCertificate
import FloatLib.Floats.Formats.BinaryInterchange.Configured.KernelPath
import FloatLib.Tactic.FloatDecide
import FloatLib.Tactic.Interval
import FloatLib.Tactic.Roundoff

/-!
# Binary64 evaluation of flat ΛCDM

The program rounds each primitive operation to nearest, ties to even. The
symbolic roundoff theorem follows the arithmetic on stored inputs. For a
particular rational input, the executable path also includes the conversion of
all three parameters; its exact decoded output gives a tight absolute bound.
-/

namespace HzRoundoff

open FloatLib.Floats FloatLib.Floats.Formats.BinaryInterchange

/-- Rounded radicand, stored explicitly so its square-root domain and error can be checked. -/
def radicand64 (Ωm z : Binary64) : Binary64 :=
  let a := 1 + z
  (1 - Ωm) + Ωm * ((a * a) * a)

/-- Round the product after the square root separately. -/
def hubble64 (h₀ Ωm z : Binary64) : Binary64 :=
  h₀ * ExecFloat.sqrt (radicand64 Ωm z)

/-- Identify the exact evaluation tree inferred by FloatLib with the Physlib definition. -/
theorem real_reference (h₀ Ωm z : ℝ) :
    h₀ * Real.sqrt ((1 - Ωm) + Ωm * (((1 + z) * (1 + z)) * (1 + z))) =
      HzCertificate.modelH h₀ Ωm z := by
  rw [HzCertificate.modelH_eq_unconditional]
  simp only [pow_succ, pow_zero, one_mul]

/-- Range and domain premises, finiteness, and a mixed bound against Physlib on stored inputs. -/
roundoff mixed theorem hubble64_roundoff using real_reference :=
  fun (h₀ Ωm z : Binary64) => ExecFloat.Binary.toReal (hubble64 h₀ Ωm z)

/-- Round each exact rational input once before evaluating the binary64 program. -/
def evaluate (h₀ Ωm z : ℚ) : Binary64 :=
  hubble64 (ExecFloat.Binary.ofRat h₀) (ExecFloat.Binary.ofRat Ωm) (ExecFloat.Binary.ofRat z)

/-- Convert a rational to nearest-even, including the absolute underflow allowance. -/
theorem input_rounding (q : ℚ) (hf : (ExecFloat.Binary.ofRat q : Binary64).isFinite = true) :
    FloatLib.Numerics.MixedError (ExecFloat.Binary.ofRat q : Binary64).toReal (q : ℝ)
      (Model.unitRoundoffAt FloatFormat.binary64) (Model.underflowErrorAt FloatFormat.binary64) := by
  rw [ExecFloat.Binary.toReal_ofRat_eq_roundAt q hf]
  exact Model.mixedError_roundAt _ _

attribute [local roundoff_rule] input_rounding

/-- End-to-end symbolic bound, including conversion of the three mathematical inputs. -/
roundoff mixed theorem evaluate_roundoff using real_reference :=
  fun (h₀ Ωm z : ℚ) => (evaluate h₀ Ωm z).toReal

/-- Finiteness is explicit because the total observation sends exceptional values to zero. -/
def Accurate (h₀ Ωm z : ℚ) (ε : ℝ) : Prop :=
  (evaluate h₀ Ωm z).isFinite = true ∧
    |(evaluate h₀ Ωm z).toReal - HzCertificate.modelH h₀ Ωm z| ≤ ε

/-- An exact decoding rules out NaN/infinity. The interval is for the original rational inputs. -/
theorem accurate_of_decode {h₀ Ωm z q : ℚ} {ε : ℝ}
    (hd : (evaluate h₀ Ωm z).toRat? = some q)
    (h : (q : ℝ) - ε ≤ HzCertificate.modelH h₀ Ωm z ∧
        HzCertificate.modelH h₀ Ωm z ≤ (q : ℝ) + ε) :
    Accurate h₀ Ωm z ε := by
  refine ⟨Model.Algebraic.isFinite_of_toRat?_eq_some hd, ?_⟩
  rw [ExecFloat.Binary.toReal_of_toRat?_eq_some hd, abs_sub_le_iff]
  exact ⟨sub_le_comm.mp h.1, sub_le_iff_le_add'.mpr h.2⟩

/-- Verify the full execution and both squared endpoint comparisons for a closed input. -/
macro "certify_hz64 " q:term : tactic =>
  `(tactic|
    (apply HzRoundoff.accurate_of_decode (q := $q)
     · float_decide [HzRoundoff.evaluate, HzRoundoff.hubble64, HzRoundoff.radicand64]
     · apply HzCertificate.modelH_bounds <;> norm_num))

private theorem unit64 : Model.unitRoundoffAt .binary64 = 1 / 2^53 := by
  rw [Model.unitRoundoffAt_eq]
  change (2 : ℝ)^(-53 : ℤ) = _
  norm_num

private theorem tiny64 : Model.underflowErrorAt .binary64 = 1 / 2^1075 := by
  rw [Model.underflowErrorAt_eq]
  change (2 : ℝ)^(-1075 : ℤ) = _
  rw [zpow_neg]
  simp

set_option exponentiation.threshold 1024 in
private theorem benchmark_range : (10^9 : ℝ) ≤ Model.toReal (Model.posMaxFinite .binary64) := by
  rw [Model.toReal_posMaxFinite_of_eq .binary64 (significand := 9007199254740991)
    (exponent := 971) (by decide) (by decide)]
  norm_num

/-- Instantiate the forward roundoff theorem for the benchmark inputs. Check the proposed
    rounded radicand, prove all range and sign conditions and compare its mixed budget with
    the requested absolute tolerance. For these benchmarks 10⁹ suffices as a range witness. -/
macro "certify_hz64_forward " h:term:arg m:term:arg z:term:arg rad:term:arg : tactic =>
  `(tactic| (
    have hd : (radicand64 (ExecFloat.Binary.ofRat $m) (ExecFloat.Binary.ofRat $z)).toReal =
        (($rad : ℚ) : ℝ) := by
      apply ExecFloat.Binary.toReal_of_toRat?_eq_some
      float_decide [radicand64]
    unfold radicand64 at hd
    dsimp only at hd
    have bound := evaluate_roundoff $h $m $z
    dsimp only at bound
    simp only [radicand64, hd, unit64, tiny64, HzCertificate.modelH_eq_unconditional,
      FloatLib.Numerics.sqrtError] at bound
    have cert := bound (by
      repeat' apply And.intro
      all_goals first | float_decide | (apply le_trans _ benchmark_range; interval (precision := 96)))
    refine ⟨cert.1, ?_⟩
    rw [HzCertificate.modelH_eq_unconditional]
    apply cert.2.bound.trans
    interval (precision := 96)))

#print axioms hubble64_roundoff
#print axioms evaluate_roundoff
#print axioms accurate_of_decode

end HzRoundoff
