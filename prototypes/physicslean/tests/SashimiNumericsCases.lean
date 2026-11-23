import SashimiTime
import SashimiIntegrable
import SashimiErf

/-! Kernel checks of a full-cell quadrature, width rejection, and the
numerically distant negative and positive inputs encountered in halo formation. -/
open FloatLib.Numerics FloatLib.Numerics.Interval SashimiNumerics

private def invScale : Expr := .unary .inv (.const 1 + .var 0)

/-- Eight panels enclose the physical elapsed time for unit H in pure Λ,
whose integrand is 1/(1+z) on [0,1]. -/
example : RationalInterval.Contains
    ⟨11388051673/17179869184, 24923587001/34359738368⟩ (SashimiTime.elapsed 1 0 1 0 1) := by
  have h := (integrate_sound_auto (e := invScale) (a := 0) (b := 1) (depth := 3)
    (bits := 32) (tol := 1/10) (degree := 24)
    (I := ⟨11388051673/17179869184, 24923587001/34359738368⟩) (by decide +kernel)).1
  convert h using 1
  apply intervalIntegral.integral_congr
  intro t ht
  have ht' : 0 ≤ t := (by simpa using ht : t ∈ Set.Icc (0:ℝ) 1).1
  rw [SashimiTime.rate, HzCertificate.modelH_eq _ _ _ ht']
  simp [invScale, Expr.eval, UnaryOp.eval, BinaryOp.eval, one_div]

example : integrate? invScale 0 1 3 32 (1/1000) = none := by decide +kernel

example : SashimiTime.enclose 70 0 SashimiTime.gyrFactor 0 none (1/100) = none := by decide +kernel

set_option maxRecDepth 20000 in
example : SashimiErf.erf? ⟨-4/5,-4/5⟩ 56 24 =
    some ⟨-53474010050059287/72057594037927936, -13368502512514821/18014398509481984⟩ := by
  decide +kernel

set_option maxRecDepth 20000 in
set_option maxHeartbeats 8000000 in
example : SashimiErf.erf? ⟨10000,10000⟩ 56 24 =
    some ⟨72057594037817149/72057594037927936, 1⟩ := by decide +kernel
