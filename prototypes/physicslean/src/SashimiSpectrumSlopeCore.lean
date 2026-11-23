import SashimiSpectrumWindow
import Mathlib.Analysis.Calculus.ParametricIntervalIntegral

/-! Analytic and integrated log-mass slope of the smooth window. -/
namespace SashimiSpectrum
open Real Set Filter MeasureTheory
open scoped Topology
noncomputable section

/-- Analytic Fourier top-hat derivative at positive argument, in the
trigonometric form used by C (equivalently −3 j₂(x)/x in F). -/
def topHatPrime (x : ℝ) := 9*cos x/x^3 - 9*sin x/x^4 + 3*sin x/x^2

theorem topHat_hasDerivAt (x : ℝ) (hx : 0 < x) : HasDerivAt topHat (topHatPrime x) x := by
  have hd := (((Real.hasDerivAt_sin x).sub ((hasDerivAt_id' x).mul (Real.hasDerivAt_cos x))).const_mul
    (3 : ℝ)).div ((hasDerivAt_id' x).pow 3) (pow_ne_zero 3 hx.ne')
  apply (hd.congr_of_eventuallyEq (f₁ := topHat) ?_).congr_deriv
  · unfold topHatPrime
    dsimp
    field_simp; ring
  · filter_upwards [Ioi_mem_nhds hx] with t ht
    simp [topHat, (show 0 < t from ht).ne']

/-- x = k R(M)/c with M=exp u; B is R(1)/c in the chosen units. -/
def logMassArgument (B u s : ℝ) := exp s * B * exp (u/3)

theorem logMassArgument_eq (ρ c u s : ℝ) (hr : 0 < ρ) :
    logMassArgument (massRadius ρ 1 / c) u s = exp s * massRadius ρ (exp u) / c := by
  unfold massRadius logMassArgument
  rw [mul_one, show 3 * exp u / (4*π*ρ) = (3/(4*π*ρ)) * exp u by ring,
    Real.mul_rpow (by positivity) (exp_pos _).le, ← Real.exp_mul]
  ring_nf

def smoothLogIntegrand (P : ℝ → ℝ) (B u s : ℝ) :=
  delta2 P (exp s) * topHat (logMassArgument B u s)^2

def smoothLogPrime (P : ℝ → ℝ) (B u s : ℝ) :=
  let x := logMassArgument B u s
  delta2 P (exp s) * 2 * topHat x * topHatPrime x * x / 3

theorem smoothLogIntegrand_hasDerivAt (P : ℝ → ℝ) (B u s : ℝ) (hB : 0 < B) :
    HasDerivAt (fun v => smoothLogIntegrand P B v s) (smoothLogPrime P B u s) u := by
  have hx : 0 < logMassArgument B u s := by unfold logMassArgument; positivity
  have hd : HasDerivAt (fun v => logMassArgument B v s) (logMassArgument B u s / 3) u := by
    unfold logMassArgument
    convert ((((hasDerivAt_id' u).div_const 3).exp).const_mul (exp s*B)) using 1; ring
  have h := (((topHat_hasDerivAt _ hx).comp u hd).pow 2).const_mul (delta2 P (exp s))
  convert h using 1
  · rfl
  · unfold smoothLogPrime
    dsimp; ring

end
end SashimiSpectrum
