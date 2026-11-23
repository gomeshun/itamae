import HzCertificate
import Physlib.Cosmology.FLRW.Distances
import SashimiNumerics

/-! FLRW background functions used by the SASHIMI clients. The four density
coefficients are independent so the WMAP7 legacy convention can be stated
exactly; only a normalized set of coefficients has E(0) = 1. -/
namespace SashimiBackground

open FloatLib.Numerics FloatLib.Numerics.Interval
open Cosmology.FLRW.FriedmannEquation

/-- Coefficients of Λ, nonrelativistic matter, radiation and curvature at the reference epoch. -/
structure Densities where
  lambda : ℚ
  matter : ℚ
  radiation : ℚ := 0
  curvature : ℚ := 0
deriving Repr

def flat (m : ℚ) : Densities := ⟨1 - m, m, 0, 0⟩

/-- In W's legacy g, ΩΛ = 1 − Ωc − Ωr with Ωc = Ωm − Ωb. -/
def wLegacy (m b r : ℚ) : Densities := ⟨1 - (m - b) - r, m, r, 0⟩

/-- W's historical Hz omits the radiation term but uses the same Λ coefficient. -/
def wLegacyHz (m b r : ℚ) : Densities := { wLegacy m b r with radiation := 0 }

/-- The second entry of historical pOmega is Ωr, although Omegaz multiplies
it by (1+z)². Store that coefficient in the quadratic slot to replay this call. -/
def wLegacyOmegaz (m b r : ℚ) : Densities := { wLegacyHz m b r with curvature := r }

def e2 (d : Densities) (z : ℚ) : ℚ :=
  d.lambda + d.matter * (1 + z)^3 + d.radiation * (1 + z)^4 + d.curvature * (1 + z)^2

/-- The historical W coefficients give E²(0)=1+Ωb in the g function. -/
theorem wLegacy_g_zero (m b r : ℚ) : e2 (wLegacy m b r) 0 = 1+b := by
  simp [e2, wLegacy]
  ring

private theorem inversePower (a : ℝ) (n : ℕ) :
    (1 / a) ^ (-(n : ℝ)) = a ^ n := by
  rw [Real.rpow_neg_eq_inv_rpow, Real.rpow_natCast]
  simp [one_div]

/-- This polynomial, with independent curvature and radiation, is the
radicand of Physlib's exact reduced Hubble function. -/
theorem reducedHubble_eq (d : Densities) (z : ℚ) :
    reducedHubble d.lambda d.matter d.radiation d.curvature (1 / (1 + (z : ℝ))) =
      Real.sqrt (e2 d z : ℚ) := by
  unfold reducedHubble e2
  push_cast
  rw [show (-3 : ℝ) = -((3 : ℕ) : ℝ) by norm_num, inversePower,
      show (-4 : ℝ) = -((4 : ℕ) : ℝ) by norm_num, inversePower,
      show (-2 : ℝ) = -((2 : ℕ) : ℝ) by norm_num, inversePower]

inductive Quantity where
  | expansionSquared | hubble | omegaMatter | criticalDensity | matterDensity | virialOverdensity
deriving Repr

/-- Hubble and gravity scalars must be in compatible units. ITAMAE uses
km/s/Mpc for H, and Mpc (km/s)²/Msun for G, hence Msun/Mpc³ for densities. -/
noncomputable def value (q : Quantity) (d : Densities) (h₀ G z : ℚ) : ℝ :=
  let s : ℝ := e2 d z
  let mass : ℝ := d.matter * (1 + z : ℚ)^3
  let x := mass / s - 1
  match q with
  | .expansionSquared => s
  | .hubble => (h₀ : ℝ) * reducedHubble d.lambda d.matter d.radiation d.curvature (1 / (1 + (z : ℝ)))
  | .omegaMatter => mass / s
  | .criticalDensity => 3 * (h₀ : ℝ)^2 * s / (8 * Real.pi * G)
  | .matterDensity => 3 * (h₀ : ℝ)^2 * mass / (8 * Real.pi * G)
  | .virialOverdensity => 18 * Real.pi^2 + 82 * x - 39 * x^2

/-- At any FLRW time with this Hubble rate, the scalar critical density agrees
with Physlib's density definition, using the same Newton constant. -/
theorem criticalDensity_eq_physlib {U : Cosmology.FLRW} (t : Time) (d : Densities)
    (h₀ G z : ℚ) (hs : 0 ≤ e2 d z)
    (hH : hubbleConstant U.a t = value .hubble d h₀ G z) :
    value .criticalDensity d h₀ G z = U.criticalDensity G t := by
  simp only [value, Cosmology.FLRW.criticalDensity, hH, reducedHubble_eq, mul_pow]
  rw [Real.sq_sqrt (by exact_mod_cast hs)]
  ring

/-- Bryan–Norman uses x=Ωm(z)−1. This fits the dimensionless virial overdensity
relative to critical density. It is separate from the EPS barrier 1.686/D. -/
def expression (q : Quantity) (d : Densities) (h₀ G z : ℚ) : Expr :=
  let s := e2 d z
  let mass := d.matter * (1 + z)^3
  let x := mass / s - 1
  let pi : Expr := .const 4 * .unary .atan (.const 1)
  match q with
  | .expansionSquared => .const s
  | .hubble => .const h₀ * .unary .sqrt (.const s)
  | .omegaMatter => .const (mass / s)
  | .criticalDensity => .binary .div (.const (3 * h₀^2 * s)) (.const (8 * G) * pi)
  | .matterDensity => .binary .div (.const (3 * h₀^2 * mass)) (.const (8 * G) * pi)
  | .virialOverdensity => .const 18 * .unary (.pow 2) pi + .const (82 * x - 39 * x^2)

private theorem eval_mul (a b : Expr) (v : Nat → ℝ) :
    (a*b).eval v = a.eval v * b.eval v := rfl
private theorem eval_add (a b : Expr) (v : Nat → ℝ) :
    (a+b).eval v = a.eval v + b.eval v := rfl

theorem expression_value (q : Quantity) (d : Densities) (h₀ G z : ℚ) :
    (expression q d h₀ G z).eval (fun _ => 0) = value q d h₀ G z := by
  cases q <;> simp only [expression, value, eval_mul, eval_add, Expr.eval, UnaryOp.eval, BinaryOp.eval,
    reducedHubble_eq] <;> push_cast <;> norm_num [Real.arctan_one] <;> ring

/-- Failure includes invalid physical inputs or an enclosure wider than the requested tolerance. -/
def enclose (q : Quantity) (d : Densities) (h₀ G z tol : ℚ)
    (bits : Nat := 64) (degree : Nat := 24) : Option RationalInterval :=
  if 0 < h₀ ∧ 0 < G ∧ 0 ≤ z ∧ 0 ≤ d.matter ∧ 0 ≤ d.radiation ∧ 0 < e2 d z then
    SashimiNumerics.closed? (expression q d h₀ G z) bits tol degree
  else none

theorem enclose_sound {q d h₀ G z tol bits degree I}
    (h : enclose q d h₀ G z tol bits degree = some I) :
    I.Contains (value q d h₀ G z) ∧ I.hi - I.lo ≤ tol := by
  unfold enclose at h
  split at h
  · exact expression_value q d h₀ G z ▸ SashimiNumerics.closed_sound h
  · contradiction

open Cosmology Cosmology.FLRW in
/-- Transport the accepted numerical interval to a Friedmann solution with matter,
radiation, curvature and Λ. Input matching includes the reference epoch and units. -/
theorem friedmann_hubble_enclosure {U : FLRW} {F Fm Fr : CosmicFluid} {Λ : ℝ}
    {c : SpeedOfLight} {d : Densities} {h₀ G z tol : ℚ} {bits degree : Nat}
    {I : RationalInterval} {t t₀ : Time} (sol : U.FriedmannSolution F Λ G c)
    (splitDensity : ∀ s, F.ρ s = Fm.ρ s + Fr.ρ s)
    (ht : t ∈ U.I) (ht₀ : t₀ ∈ U.I)
    (expanding₀ : 0 < hubbleConstant U.a t₀) (expanding : 0 ≤ hubbleConstant U.a t)
    (matter : Fm.ρ t = Fm.ρ t₀ * (U.a t₀ / U.a t)^3)
    (radiation : Fr.ρ t = Fr.ρ t₀ * (U.a t₀ / U.a t)^4)
    (refH : hubbleConstant U.a t₀ = h₀) (redshift : U.redshift t₀ t = z)
    (hΛ : U.lambdaDensityParameter Λ c t₀ = d.lambda)
    (hm : U.densityParameter Fm G t₀ = d.matter)
    (hr : U.densityParameter Fr G t₀ = d.radiation)
    (hk : U.curvatureDensityParameter c t₀ = d.curvature)
    (cert : enclose .hubble d h₀ G z tol bits degree = some I) :
    I.Contains (hubbleConstant U.a t) ∧ I.hi - I.lo ≤ tol := by
  have solH := sol.hubbleConstant_eq_reducedHubble splitDensity ht ht₀ expanding₀
    expanding matter radiation
  have scale : 1 / (1 + (z : ℝ)) = U.a t / U.a t₀ := by
    rw [← redshift, U.one_add_redshift, one_div, inv_div]
  have h : hubbleConstant U.a t = value .hubble d h₀ G z := by
    simpa only [value, refH, hΛ, hm, hr, hk, scale] using solH
  rw [h]
  exact enclose_sound cert

#print axioms friedmann_hubble_enclosure
#print axioms criticalDensity_eq_physlib
#print axioms reducedHubble_eq
#print axioms enclose_sound

end SashimiBackground
