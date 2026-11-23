import SashimiNumerics
import SashimiGrowth
import SashimiSpectrumFit
import Mathlib.Analysis.Calculus.Deriv.MeanValue
import Mathlib.Analysis.SpecialFunctions.Log.Deriv

/-! # Halo quantities used by SASHIMI-C/SI/W/F
The expressions below run with FloatLib's certified binary-grid evaluator. Inputs may themselves
be expressions, so intermediate enclosures (including variance and growth) can propagate. The
empirical coefficients specify the adopted fits; they carry no claim about calibration error. -/

namespace SashimiHalo

open FloatLib.Numerics FloatLib.Numerics.Interval

local instance (n : Nat) : OfNat Expr n := ⟨.const n⟩
local instance : OfScientific Expr := ⟨fun m s e => .const (OfScientific.ofScientific m s e : ℚ)⟩
local instance : Sub Expr := ⟨.binary .sub⟩
local instance : Div Expr := ⟨.binary .div⟩
local instance : Neg Expr := ⟨.unary .neg⟩
local instance : HPow Expr Nat Expr := ⟨fun x n => .unary (.pow n) x⟩

private def log (x : Expr) := x.unary .log
private def exp (x : Expr) := x.unary .exp
private def sqrt (x : Expr) := x.unary .sqrt
/-- Real powers on strictly positive bases, evaluated as exp(p log x). -/
def posPower (x p : Expr) := exp (p * log x)
private def log10 (x : Expr) := log x / log 10
/-- The exact identity π = 4 atan(1). -/
def piE : Expr := 4 * Expr.unary .atan 1

theorem piE_eval (env : Nat → ℝ) : piE.eval env = Real.pi := by
  change 4 * Real.arctan 1 = _
  rw [Real.arctan_one]
  ring

/-- Dimensionless NFW mass at x = r / r_s. -/
def fc (x : Expr) := log (1 + x) - x / (1 + x)

def density (r rs ρs : Expr) :=
  let x := r / rs
  ρs / (x * (1 + x) ^ 2)

def enclosedMass (r rs ρs : Expr) := 4 * piE * ρs * rs ^ 3 * fc (r / rs)

/-- Potential in (km/s)^2 for Mpc and solar masses; use the analytic limit at r = 0. -/
def potential (r : ℚ) (rs ρs : Expr) :=
  let x := Expr.const r / rs
  (-4) * piE * 0.00000000430091 * ρs * rs ^ 2 * (if r = 0 then 1 else log (1 + x) / x)

/-- Bryan–Norman critical-density virial overdensity, flat matter plus Λ, input Ωm(z). -/
def virialOverdensity (Ωz : Expr) :=
  let x := Ωz - 1
  18 * piE ^ 2 + 82 * x - 39 * x ^ 2

/-- Radius at the named mean overdensity relative to the supplied critical density. -/
def overdensityRadius (mass Δ ρcrit : Expr) := posPower (3 * mass / (4 * piE * Δ * ρcrit)) (1 / 3)

def scaleDensity (c Δ ρcrit : Expr) := Δ * ρcrit * c ^ 3 / (3 * fc c)

/-- In units of r_s, the exact NFW conversion solves meanDensity(c) = Δ rhoCrit. -/
def meanDensity (c ρs : Expr) := 3 * ρs * fc c / c ^ 3

/-- Correa et al. (2015) CDM concentration: specify M200/Msun, z ≥ 0. The calibration
uses one formula through z = 4 and another above; WDM uses a separate Ludlow prescription. -/
def correa200 (mass : Expr) (z : ℚ) : Expr := Id.run do
  let u := 1 + Expr.const z
  let l := log10 mass
  if z ≤ 4 then
    let a := 1.7543 - 0.2766 * u + 0.02039 * u ^ 2
    let b := 0.2753 + 0.00351 * u - 0.3038 * posPower u 0.0269
    let g := -0.01537 + 0.02102 * posPower u (-0.1475)
    return posPower 10 (a + b * l * (1 + g * l ^ 2))
  else
    let a := 1.3081 - 0.1078 * u + 0.00398 * u ^ 2
    let b := 0.0223 - 0.0944 * posPower u (-0.3907)
    return posPower 10 (a + b * l)

/-- Hu–Kravtsov approximation for Mvir given M200, c200 and virial Δ. The residual of
the exact NFW equation is not bounded by certification of the fitting formula. -/
def virialMassFit (mass c Δ : Expr) :=
  let f := Δ / 200 * fc c / c ^ 3
  let p := -0.4283 - 0.00313 * log f - 0.0000352 * (log f) ^ 2
  let x := 1 / sqrt (0.5116 * posPower f (2 * p) + (3 / 4) ^ 2) + 2 * f
  Δ / 200 * mass / (c * x) ^ 3

/-- Correa formation-redshift fit, with a positive formation redshift required for q. -/
def formationRedshift (mass : Expr) :=
  let l := log10 mass
  (-0.0064) * l ^ 2 + 0.0237 * l + 1.8837

def formationMassRatio (mass : Expr) := 4.137 / posPower (formationRedshift mass) 0.9476

/-- Shared mass-history fit for Mzi/Mzzi. mᵢ, Dᵢ, D'ᵢ and S(mᵢ/q,zᵢ)−S(mᵢ,zᵢ) come
from the selected model. Mzi uses zᵢ=0 and Dᵢ=1. Positive variance difference required. -/
def massHistory (mᵢ dz Dᵢ dDᵢ ΔS : Expr) :=
  let f := 1 / sqrt ΔS
  let α := f * (1 + 1.686 * sqrt (2 / piE) * dDᵢ / Dᵢ ^ 2)
  mᵢ * posPower (1 + dz) α * exp (-f * dz)

/-- SASHIMI uses the fitted virial mass in its accretion-rate expression. -/
def massRate (virialMass dz Dᵢ dDᵢ ΔS : Expr) :=
  let f := 1 / sqrt ΔS
  let α := f * (1 + 1.686 * sqrt (2 / piE) * dDᵢ / Dᵢ ^ 2)
  (-f + α / (1 + dz)) * virialMass

/-- W's fitted matter/radiation equality redshift, not the exact Ωm/Ωr equality. -/
def wdmEquality (Ω h : Expr) := 3600 * Ω * h ^ 2 / 0.15 - 1

/-- SASHIMI-W Benson et al. Jeans mass in Msun, thermal relic mass in keV, gx=1.5. -/
def jeansMass (Ω h mkeV : Expr) :=
  306000000 * posPower ((1 + wdmEquality Ω h) / 3000) 1.5 * sqrt (Ω * h ^ 2 / 0.15) / mkeV ^ 4

def wdmCollapse (mass mJeans growth : Expr) :=
  let x := log (mass / mJeans)
  let s := 1 / (1 + exp ((x + 2.4) / 0.1))
  1.686 / growth * (s * (0.04 / exp (2.3 * x)) + (1 - s) * exp (0.31687 / exp (0.809 * x)))

inductive AssemblyOutput where
  | m200 | virialRate

/-- C/SI/F-style Mzi/Mzzi/dMdz assembly with the variance, growth, rate and virial fit explicitly
selected. The closed C/SI default below chooses the Ludlow CDM variance. FDM and WDM use their
own S(M,z) instead. The first stage uses D=1 exactly as the source's Mzi coefficient does. -/
def assemblyUsing (out : AssemblyOutput) (m0 : Expr) (z zi : ℚ)
    (variance : Expr → ℚ → Expr) (growth rate : ℚ → Expr) (toVirial : Expr → Expr) : Expr :=
  let q := formationMassRatio m0
  let varianceDiff := fun m t => variance (m/q) t - variance m t
  let mᵢ := if zi = 0 then m0 else massHistory m0 (.const zi) 1 (rate 0) (varianceDiff m0 0)
  let dt : Expr := .const (z-zi)
  let dS := varianceDiff mᵢ zi
  let m := massHistory mᵢ dt (growth zi) (rate zi) dS
  match out with
  | .m200 => m
  | .virialRate => massRate (toVirial m) dt (growth zi) (rate zi) dS

/-- Complete dark-matter-only baseline for SASHIMI-C (prompt cusps disabled) and SI.
M0 is in Msun. Every intermediate mass, concentration and variance participates in the
interval evaluation. The output rate uses the *fitted* virial conversion, as in the source. -/
def coldAssembly (out : AssemblyOutput) (Ω h m0 z zi : ℚ) : Expr :=
  let D := SashimiGrowth.expression .growth Ω (1-Ω) h
  let dD := SashimiGrowth.expression .growthDerivative Ω (1-Ω) h
  let s := fun m t => (SashimiSpectrum.ludlowExpr (.const (10^10/h)) m * D t)^2
  let Δ := virialOverdensity (.const (Ω*(1+z)^3 / ((1-Ω)+Ω*(1+z)^3)))
  assemblyUsing out (.const m0) z zi s D dD (fun m => virialMassFit m (correa200 m z) Δ)

/-- Guarded entry points for standalone quantities. Masses are in Msun, radii in Mpc,
densities in Msun/Mpc³. M200 and its concentration use the critical-density convention.
The mass-history node consumes the positive difference of variances and the chosen D, D'. -/
def expression? (name : String) (values : List ℚ) : Option Expr :=
  let q := Expr.const
  match name, values with
  | "nfw-f", [x] => if 0 ≤ x then some (fc (q x)) else none
  | "nfw-density", [r, rs, ρs] =>
    if 0 < r ∧ 0 < rs ∧ 0 < ρs then some (density (q r) (q rs) (q ρs)) else none
  | "nfw-mass", [r, rs, ρs] =>
    if 0 ≤ r ∧ 0 < rs ∧ 0 < ρs then some (enclosedMass (q r) (q rs) (q ρs)) else none
  | "nfw-potential", [r, rs, ρs] =>
    if 0 ≤ r ∧ 0 < rs ∧ 0 < ρs then some (potential r (q rs) (q ρs)) else none
  | "virial-overdensity", [Ωz] =>
    if 0 < Ωz ∧ Ωz ≤ 1 then some (virialOverdensity (q Ωz)) else none
  | "overdensity-radius", [m, Δ, ρcrit] =>
    if 0 < m ∧ 0 < Δ ∧ 0 < ρcrit then some (overdensityRadius (q m) (q Δ) (q ρcrit)) else none
  | "nfw-scale-density", [c, Δ, ρcrit] =>
    if 0 < c ∧ 0 < Δ ∧ 0 < ρcrit then some (scaleDensity (q c) (q Δ) (q ρcrit)) else none
  | "nfw-mean-density", [c, ρs] =>
    if 0 < c ∧ 0 < ρs then some (meanDensity (q c) (q ρs)) else none
  | "c200-correa", [m, z] => if 0 < m ∧ 0 ≤ z then some (correa200 (q m) z) else none
  | "mvir-fit", [m, c, Δ] =>
    if 0 < m ∧ 0 < c ∧ 0 < Δ then some (virialMassFit (q m) (q c) (q Δ)) else none
  | "formation-redshift", [m] => if 0 < m then some (formationRedshift (q m)) else none
  | "formation-q", [m] => if 0 < m then some (formationMassRatio (q m)) else none
  | "mass-history", [m, dz, D, dD, dS] =>
    if 0 < m ∧ 0 ≤ dz ∧ 0 < D ∧ 0 < dS then
      some (massHistory (q m) (q dz) (q D) (q dD) (q dS)) else none
  | "mass-rate", [mVir, dz, D, dD, dS] =>
    if 0 < mVir ∧ 0 ≤ dz ∧ 0 < D ∧ 0 < dS then
      some (massRate (q mVir) (q dz) (q D) (q dD) (q dS)) else none
  | "cold-m200", [Ω, h, m, z, zi] =>
    if 0 < Ω ∧ Ω ≤ 1 ∧ 0 < h ∧ 0 < m ∧ 0 ≤ zi ∧ zi ≤ z then
      some (coldAssembly .m200 Ω h m z zi) else none
  | "cold-virial-rate", [Ω, h, m, z, zi] =>
    if 0 < Ω ∧ Ω ≤ 1 ∧ 0 < h ∧ 0 < m ∧ 0 ≤ zi ∧ zi ≤ z then
      some (coldAssembly .virialRate Ω h m z zi) else none
  | "wdm-z-eq", [Ω, h] => if 0 < Ω ∧ Ω ≤ 1 ∧ 0 < h then some (wdmEquality (q Ω) (q h)) else none
  | "wdm-jeans", [Ω, h, m] =>
    if 0 < Ω ∧ Ω ≤ 1 ∧ 0 < h ∧ 0 < m then some (jeansMass (q Ω) (q h) (q m)) else none
  | "wdm-collapse", [m, mj, D] =>
    if 0 < m ∧ 0 < mj ∧ 0 < D then some (wdmCollapse (q m) (q mj) (q D)) else none
  | _, _ => none

def evaluate? (name : String) (values : List ℚ) (bits : Nat) (tol : ℚ) (degree : Nat := 24) :=
  (expression? name values).bind fun e => SashimiNumerics.closed? e bits tol degree

theorem evaluate_sound {name values bits tol degree I} (h : evaluate? name values bits tol degree = some I) :
    ∃ e, expression? name values = some e ∧ I.Contains (e.eval fun _ => 0) ∧ I.hi - I.lo ≤ tol := by
  obtain ⟨e, he, hI⟩ := Option.bind_eq_some_iff.mp h
  exact ⟨e, he, SashimiNumerics.closed_sound hI⟩

end SashimiHalo
