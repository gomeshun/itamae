import SashimiBackground
import SashimiTime
import SashimiHalo
import Physlib.Cosmology.FLRW.Distances

/-! # FLRW and halo scalars consumed by the SASHIMI clients
The native background uses kg-free gravitational units Mpc (km/s)²/Msun, h = H₀/100 and
the native Gyr conversion. The tidal prescriptions retain the distinct numerical fits and
year conventions of SASHIMI-C/SI/F and W. -/
namespace SashimiDerived

open FloatLib.Numerics FloatLib.Numerics.Interval
open SashimiBackground

local instance (n : Nat) : OfNat Expr n := ⟨.const n⟩
local instance : Sub Expr := ⟨.binary .sub⟩
local instance : Div Expr := ⟨.binary .div⟩
local instance : Neg Expr := ⟨.unary .neg⟩

/-- a(t)/a(t₀), for the same emission/observation times as Physlib.redshift. -/
theorem scaleFactor_of_redshift (U : Cosmology.FLRW) (t₀ t : Time) :
    1 / (1 + U.redshift t₀ t) = U.a t / U.a t₀ := by
  rw [U.one_add_redshift]
  simp [one_div]

def scaleFactor (z : ℚ) : ℚ := 1/(1+z)
def lambdaFraction (m z : ℚ) : ℚ := (1-m)/e2 (flat m) z

theorem lambdaFraction_add_matter {m z : ℚ} (hs : e2 (flat m) z ≠ 0) :
    lambdaFraction m z + m*(1+z)^3 / e2 (flat m) z = 1 := by
  unfold lambdaFraction
  rw [← add_div, div_eq_one_iff_eq hs]
  simp [e2,flat]

/-- Ωb refers to present-day baryons. This cold-dark-matter partition assumes that the balance
of the matter component clusters as CDM; a massive-neutrino component would need another term. -/
def cdmDensity (m b h G z : ℚ) : Expr :=
  SashimiBackground.expression .matterDensity (flat (m-b)) (100*h) G z

theorem cdm_add_baryon {m b h G z : ℚ} :
    (cdmDensity m b h G z).eval (fun _ => 0) +
      value .matterDensity (flat b) (100*h) G z =
      value .matterDensity (flat m) (100*h) G z := by
  rw [cdmDensity, expression_value]
  simp only [value,flat]
  push_cast
  ring

/-- The signed cosmic-time derivative in Gyr per increasing redshift. -/
def signedTime (m h z : ℚ) : Expr :=
  -((1 : Expr) / (.const (1+z) * expression .hubble (flat m) (100*h) 1 z *
    .const SashimiTime.gyrFactor))

theorem signedTime_eq (m h z : ℚ) : (signedTime m h z).eval (fun _ => 0) =
    -SashimiTime.rate (100*h : ℚ) m SashimiTime.gyrFactor z := by
  change -(1 / ((↑(1+z : ℚ) : ℝ) * (expression _ _ _ _ _).eval (fun _ => 0) * _)) = _
  rw [expression_value]
  unfold SashimiTime.rate value HzCertificate.modelH flat
  push_cast
  rfl

/-- On the physical redshift domain the result differentiates the improper FLRW age integral.
At redshift zero the assertion is the right derivative. Positive matter ensures finite age. -/
theorem signedTime_age {m h z : ℚ} (hm : 0 < m) (hm₁ : m ≤ 1) (hh : 0 < h) (hz : 0 ≤ z) :
    HasDerivWithinAt (SashimiTime.age (100*h : ℚ) m SashimiTime.gyrFactor)
      ((signedTime m h z).eval (fun _ => 0)) (Set.Ici 0) z := by
  rw [signedTime_eq]
  have hH : (0:ℝ) < (↑(100*h : ℚ) : ℝ) := by exact_mod_cast (show (0:ℚ) < 100*h by positivity)
  have hm' : (0:ℝ) < m := by exact_mod_cast hm
  have hm₁' : (m:ℝ) ≤ 1 := by exact_mod_cast hm₁
  have hz' : (0:ℝ) ≤ z := by exact_mod_cast hz
  have hc : (0:ℝ) < SashimiTime.gyrFactor := by norm_num [SashimiTime.gyrFactor]
  apply (SashimiTime.age_deriv hH hm' hm₁' hc hz').hasDerivWithinAt.congr_of_mem
  · exact fun t ht => SashimiTime.age_eq hH hm' hm₁' hc ht
  · exact hz'

/-- Source formula for the halo dynamical time, in billions of the source profile's years.
E is H/H₀. Δ is the virial overdensity relative to critical density. -/
def dynamicalTime (h E Δ : Expr) :=
  .const 1.628 / h / Expr.unary .sqrt (Δ/178) / E

def dynamicalFlat (m h z : ℚ) := dynamicalTime (.const h) (.unary .sqrt (.const (e2 (flat m) z)))
  (expression .virialOverdensity (flat m) (100*h) 1 z)

inductive TidalFit where | csf | warm deriving DecidableEq, Repr

/-- Base-10 amplitude and exponent of the dark-matter-only stripping model. The fits take
the current host virial mass in Msun. They are inputs to the ODE and perturbation series. -/
def aCoeff (fit : TidalFit) (mass z : Expr) : Expr :=
  let L := Expr.unary .log mass / Expr.unary .log 10
  let p := match fit with
    | .csf => (.const (-0.0003)*L + .const 0.02)*z + (.const 0.011*L - .const 0.354)
    | .warm => (.const (-0.0019)*L + .const 0.045)*z + (.const 0.0097*L - .const 0.313)
  SashimiHalo.posPower 10 p

def zetaCoeff (fit : TidalFit) (mass z : Expr) : Expr :=
  let L := Expr.unary .log mass / Expr.unary .log 10
  match fit with
  | .csf => (.const 0.00012*L - .const 0.0033)*z + (.const (-0.0011)*L + .const 0.026)
  | .warm => (.const (-0.0000555)*L + .const 0.00143)*z + (.const 0.000334*L - .const 0.00811)

/-- Convert source km/s/Mpc to its reciprocal Gyr using its own adopted Mpc and year.
The W solver writes 86400*365; C, SI and F set yr=31557600 seconds. -/
def sourceGyrFactor (fit : TidalFit) : ℚ :=
  let seconds : ℚ := match fit with | .csf => 31557600 | .warm => 31536000
  seconds * 10^9 / (3086 * 10^16)

/-- Positive dm/dz prefactor A/[tdyn H (1+z)] for evolving from high to low redshift. -/
def strippingPhi (fit : TidalFit) (mass z h E Δ : Expr) : Expr :=
  aCoeff fit mass z / (dynamicalTime h E Δ * (100*h*E) * .const (sourceGyrFactor fit) * (1+z))

def expression? (name : String) (values : List ℚ) : Option Expr :=
  let q := Expr.const
  match name, values with
  | "scale-factor", [z] => if 0 ≤ z then some (q (scaleFactor z)) else none
  | "omega-lambda-flat", [m,z] =>
    if 0 ≤ m ∧ m ≤ 1 ∧ 0 ≤ z then some (q (lambdaFraction m z)) else none
  | "rho-cdm-flat", [m,b,h,G,z] =>
    if 0 ≤ b ∧ b ≤ m ∧ m ≤ 1 ∧ 0 < h ∧ 0 < G ∧ 0 ≤ z then some (cdmDensity m b h G z) else none
  | "dt-dz-flat", [m,h,z] =>
    if 0 ≤ m ∧ m ≤ 1 ∧ 0 < h ∧ 0 ≤ z then some (signedTime m h z) else none
  | "dynamical-time-flat", [m,h,z] =>
    if 0 < m ∧ m ≤ 1 ∧ 0 < h ∧ 0 ≤ z then some (dynamicalFlat m h z) else none
  | "dynamical-time", [h,E,Δ] =>
    if 0 < h ∧ 0 < E ∧ 0 < Δ then some (dynamicalTime (q h) (q E) (q Δ)) else none
  | "stripping-a-csf", [m,z] | "stripping-a-w", [m,z] =>
    if 0 < m ∧ 0 ≤ z then some (aCoeff (if name = "stripping-a-w" then .warm else .csf) (q m) (q z)) else none
  | "stripping-zeta-csf", [m,z] | "stripping-zeta-w", [m,z] =>
    if 0 < m ∧ 0 ≤ z then some (zetaCoeff (if name = "stripping-zeta-w" then .warm else .csf) (q m) (q z)) else none
  | "stripping-phi-csf", [m,z,h,E,Δ] | "stripping-phi-w", [m,z,h,E,Δ] =>
    if 0 < m ∧ 0 ≤ z ∧ 0 < h ∧ 0 < E ∧ 0 < Δ then
      some (strippingPhi (if name = "stripping-phi-w" then .warm else .csf) (q m) (q z) (q h) (q E) (q Δ))
    else none
  | _, _ => none

def evaluate? (name : String) (values : List ℚ) (bits : Nat) (tol : ℚ) (degree : Nat := 24) :=
  (expression? name values).bind fun e => SashimiNumerics.closed? e bits tol degree

theorem evaluate_sound {name values bits tol degree I} (h : evaluate? name values bits tol degree = some I) :
    ∃ e, expression? name values = some e ∧ I.Contains (e.eval fun _ => 0) ∧ I.hi - I.lo ≤ tol := by
  obtain ⟨e, he, hI⟩ := Option.bind_eq_some_iff.mp h
  exact ⟨e, he, SashimiNumerics.closed_sound hI⟩

#print axioms evaluate_sound
#print axioms signedTime_age
#print axioms scaleFactor_of_redshift

end SashimiDerived
