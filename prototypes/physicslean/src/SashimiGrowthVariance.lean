import SashimiGrowth
import SashimiMassVariance

/-! Redshift the supplied reference-epoch variance (or its signed mass derivative)
with the real Carroll–Press–Turner growth function. The input is an enclosure
of the z=0 integral, before applying a growth factor. -/
namespace SashimiGrowthVariance
open FloatLib.Numerics FloatLib.Numerics.Interval

inductive Quantity where | variance | sigma | derivative | legacyFDerivative
deriving DecidableEq, Repr

noncomputable def value (q : Quantity) (D source : ℝ) : ℝ :=
  match q with
  | .variance | .derivative => SashimiSpectrum.variance D source
  | .sigma => SashimiSpectrum.sigma D source
  | .legacyFDerivative => D * source

private def expression (q : Quantity) : Expr :=
  let d := Expr.var 0
  let raw := Expr.var 1
  match q with
  | .variance | .derivative => .unary (.pow 2) d * raw
  | .sigma => d * .unary .sqrt raw
  | .legacyFDerivative => d * raw

/-- Interval multiplication keeps signed derivative endpoints in their proper order.
Taking a square root requires a certified nonnegative source variance. -/
def combine? (q : Quantity) (D raw : RationalInterval) (tol : ℚ) (bits : Nat) :
    Option RationalInterval := do
  guard (0 ≤ D.lo ∧ (q = .sigma → 0 ≤ raw.lo))
  let out ← (expression q).gridEnclosureWith? [D,raw] bits 0
  SashimiNumerics.within? tol out

theorem combine_sound {q D raw tol bits I} {d s : ℝ} (hd : D.Contains d) (hs : raw.Contains s)
    (h : combine? q D raw tol bits = some I) : I.Contains (value q d s) ∧ I.hi-I.lo ≤ tol := by
  unfold combine? at h
  obtain ⟨_, _, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨J, hJ, hw⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨rfl, hw⟩ := SashimiNumerics.within_sound hw
  let env : Nat → ℝ := fun | 0 => d | _ => s
  have v := Expr.contains_gridEnclosureWith? (values := env) (h := hJ) (by
    intro i box hi
    cases i with
    | zero => simp at hi; subst box; exact hd
    | succ i => cases i with
      | zero => simp at hi; subst box; exact hs
      | succ i => simp at hi)
  exact ⟨by cases q <;> exact v, hw⟩

/-- The dimensionless CPT growth depends only on m, Λ and z. The unused h
argument of the common SashimiGrowth API is set to one. -/
def atRedshift? (q : Quantity) (m l z : ℚ) (raw : RationalInterval) (tol growthTol : ℚ)
    (bits degree : Nat) : Option RationalInterval := do
  let D ← SashimiGrowth.enclose .growth m l 1 z growthTol bits degree
  combine? q D raw tol bits

theorem atRedshift_sound {q m l z raw tol growthTol bits degree I} {s : ℝ}
    (hs : raw.Contains s) (h : atRedshift? q m l z raw tol growthTol bits degree = some I) :
    I.Contains (value q (SashimiGrowth.growth m l z) s) ∧ I.hi-I.lo ≤ tol := by
  obtain ⟨D, hD, hc⟩ := Option.bind_eq_some_iff.mp h
  exact combine_sound (SashimiGrowth.enclose_sound hD).1 hs hc

/-- Compose the mass integral computed at unit growth with the redshift
calculation. The spectrum refinement is the same one used by the mass API. -/
theorem integral_atRedshift_sound {q mode placement ratio cal rho mass rows a b firstTol depth
    bits degree raw observable m l z tol growthTol I}
    (P : SashimiTable.Cell → ℝ → ℝ)
    (hp : ∀ c ∈ SashimiTable.cells rows, ∀ x : Expr, ∀ v, 0 < x.eval v →
      (SashimiTable.powerExpr mode placement ratio c x).eval v = P c (x.eval v))
    (hr : SashimiMassVariance.evaluate? q mode placement ratio cal rho mass 1 rows a b firstTol
      depth bits degree = some raw)
    (hz : atRedshift? observable m l z raw tol growthTol bits degree = some I) :
    I.Contains (value observable (SashimiGrowth.growth m l z)
      (SashimiMassVariance.massBand q (fun c k => SashimiMassVariance.kernel q (P c) rho mass 1 k)
        (cal.value / SashimiSpectrum.massRadius rho mass) a b (SashimiTable.cells rows))) ∧
      I.hi-I.lo ≤ tol := by
  apply atRedshift_sound _ hz
  simpa using (SashimiMassVariance.evaluate_sound P hp hr).1

/-- At fixed redshift D has no mass dependence, so the signed variance
derivative uses the squared growth factor. -/
theorem variance_derivative {m l z M s' : ℝ} {source : ℝ → ℝ} (h : HasDerivAt source s' M) :
    HasDerivAt (fun x => value .variance (SashimiGrowth.growth m l z) (source x))
      (value .derivative (SashimiGrowth.growth m l z) s') M :=
  h.const_mul _

#print axioms atRedshift_sound
#print axioms integral_atRedshift_sound
#print axioms variance_derivative

end SashimiGrowthVariance
