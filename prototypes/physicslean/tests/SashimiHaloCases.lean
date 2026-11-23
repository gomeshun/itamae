import SashimiHaloRoots
import SashimiHaloWarm
import SashimiDerived

open SashimiHalo FloatLib.Numerics

/-- An explicitly checked numerical bracket around the unique untruncated NFW radius at f(c)=1. -/
example : ∃ c : ℝ, (5 : ℝ) ≤ c ∧ c ≤ 6 ∧ nfw c = 1 ∧ ∀ t : ℝ, 0 ≤ t → nfw t = 1 → t = c := by
  obtain ⟨_, c, hc, _, he, hu⟩ := rootBracket_sound
    (kind := .enclosed) (I := ⟨5,6⟩) (target := .const 1) (tol := 1)
    (bits := 32) (degree := 16) (by decide +kernel)
  exact ⟨c, by simpa using hc.1, by simpa using hc.2, by simpa [rootModel, FloatLib.Numerics.Interval.Expr.eval] using he, by
    simpa [domain,rootModel,FloatLib.Numerics.Interval.Expr.eval] using hu⟩

example : rootBracket? .enclosed ⟨0,1⟩ (.const 10) 1 32 16 = false := by decide +kernel
example : evaluate? "nfw-density" [0,1,2] 64 1 = none := by decide +kernel
example : evaluate? "mass-history" [1,1,1,0,-1] 64 1 = none := by decide +kernel
example : evaluate? "nfw-density" [1,1,2] 64 0 = some (.point (1/2)) := by decide +kernel

/-- Round outward to specified external limits. Every runtime regression also tests acceptance. -/
private def inside (i : Option RationalInterval) (lo hi : ℚ) : Bool :=
  i.any (fun J => lo ≤ J.lo && J.lo ≤ J.hi && J.hi ≤ hi)

#guard inside (evaluate? "cold-m200" [0.315,0.674,10^12,1,0] 96 (1/100) 36)
  574482431960.08 574482431960.10

#guard inside (evaluate? "cold-virial-rate" [0.315,0.674,10^12,1,0] 96 (1/100) 36)
  (-366064752917.98) (-366064752917.97)

#guard inside (evaluate? "cold-m200" [0.315,0.674,10^12,2,1] 96 (1/100) 36)
  195944526417.30 195944526417.31

#guard (virial? (10^12) 7 100 ⟨1,20⟩ 100 40 90 (1/10^15) (1/1000)).isSome
#guard inside (SashimiHaloWarm.concentration? 0.272 0 (.point 4) 64 24 (1/10^8))
  5.8467704701 5.8467704734
#guard (SashimiHaloWarm.concentration? 0.272 0 ⟨-1,1⟩ 64 24 1).isNone

example : SashimiDerived.evaluate? "scale-factor" [1] 16 0 = some (.point (1/2)) := by decide +kernel
example : SashimiDerived.evaluate? "rho-cdm-flat" [0.3,0.4,0.7,1,1] 16 1 = none := by decide +kernel

#guard inside (SashimiDerived.evaluate? "rho-cdm-flat" [0.315,0.049,0.674,0.00000000430091,1]
  96 (1/1000000) 30) 268294924565.057389 268294924565.057390
#guard inside (SashimiDerived.evaluate? "dt-dz-flat" [0.315,0.674,1] 64 (1/10^6))
  (-4.051751) (-4.051750)
#guard inside (SashimiDerived.evaluate? "dynamical-time-flat" [0.315,0.674,1] 64 (1/10^6))
  1.430498 1.430499
#guard inside (SashimiDerived.evaluate? "stripping-phi-csf" [10^12,1,0.674,2,178] 64 (1/10^6))
  1.870714 1.870715
#guard inside (SashimiDerived.evaluate? "stripping-phi-w" [10^12,1,0.674,2,178] 64 (1/10^6))
  2.011429 2.011430
