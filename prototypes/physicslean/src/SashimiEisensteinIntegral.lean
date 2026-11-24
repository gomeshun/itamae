import SashimiAnalyticWindow
import SashimiEisensteinHu
import SashimiMassVariance

/-! Enclosures for normalized finite-band variance of the continuous EH98 CDM+baryon
transfer, with k in h/Mpc and a Lagrangian top-hat radius in Mpc/h. -/
namespace SashimiEisensteinIntegral
open FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory SashimiSpectrum
open SashimiEisensteinHu SashimiAnalyticWindow

structure Cosmology where
  h : ℚ
  om : ℚ
  ob : ℚ
  temperature : ℚ
  ns : ℚ
  sigma8 : ℚ
  criticalDensity : ℚ
deriving Repr

/-- rho_crit,0 in (Msun/h)/(Mpc/h)^3. -/
def matched : Cosmology := ⟨0.674,0.315,0.049,2.7255,0.965,0.811,277536627245.708⟩

/-- Multiplicative constants independent of R cancel in the anchor-normalized ratio. -/
def weighted (c : Cosmology) : Expr :=
  (.unary .exp (.const (2+c.ns) * .unary .log (.var 0))) *
    .unary (.pow 2) variableTransfer

noncomputable def coefficients (c : Cosmology) : Coefficients ℝ :=
  prepareReal c.h c.om c.ob c.temperature

noncomputable def raw (c : Cosmology) (a b : ℚ) (r : ℝ) : ℝ :=
  ∫ k in (a:ℝ)..b, Real.exp (((2+c.ns : ℚ):ℝ) * Real.log k) *
    (transferReal (coefficients c) k)^2 * (topHat (k*r))^2

/-- Number multiplying the normalized spectrum for a mass in physical Msun. -/
noncomputable def radius (c : Cosmology) (m : ℚ) : ℝ :=
  massRadius ((c.criticalDensity*c.om:ℚ):ℝ) ((m*c.h:ℚ):ℝ)

def radiusBox? (c : Cosmology) (m : ℚ) (bits degree : Nat) : Option RationalInterval :=
  if 0 < c.criticalDensity*c.om ∧ 0 < m*c.h then
    (SashimiMassVariance.radiusExpr (c.criticalDensity*c.om) (m*c.h)).gridEnclosure? bits degree
  else none

theorem radiusBox_sound {c m bits degree R} (h : radiusBox? c m bits degree = some R) :
    R.Contains (radius c m) := by
  unfold radiusBox? at h
  split at h
  · rename_i hd
    simpa [SashimiMassVariance.radius_eval _ _ hd.1 hd.2, radius] using Expr.contains_gridEnclosure? h
  · contradiction

theorem raw_eq (c : Cosmology) (a b : ℚ) (r : ℝ) :
    raw c a b r = ∫ k in (a:ℝ)..b,
      integrand (weighted c) (fun i => ((coefficients c).toList[i]?).getD 0) r k := by
  rfl

/-- Positive k interpretation of the factor k^(2+ns) in the linear spectrum integral. -/
theorem tilt_eq (c : Cosmology) {k : ℝ} (hk : 0 < k) :
    Real.exp (((2+c.ns : ℚ):ℝ) * Real.log k) = k^(2+(c.ns:ℝ)) := by
  rw [Real.rpow_def_of_pos hk]
  push_cast
  rw [mul_comm]

def boxes? (c : Cosmology) (bits degree : Nat) := coefficientBoxes? c.h c.om c.ob c.temperature bits degree

theorem raw_sound {c m bits degree cb cd boxes R a ends depth J}
    (hboxes : boxes? c cb cd = some boxes) (hR : radiusBox? c m bits degree = some R)
    (h : band? (weighted c) boxes R depth bits degree a ends = some J) :
    J.Contains (raw c a (ends.getLastD a) (radius c m)) := by
  rw [raw_eq]
  exact (band_sound (fun _ _ hi => coefficientBoxes_sound hboxes hi) (radiusBox_sound hR) ends h).2

theorem anchor_sound {c bits degree cb cd boxes a ends depth J}
    (hboxes : boxes? c cb cd = some boxes)
    (h : band? (weighted c) boxes (.point 8) depth bits degree a ends = some J) :
    J.Contains (raw c a (ends.getLastD a) 8) := by
  rw [raw_eq]
  exact (band_sound (fun _ _ hi => coefficientBoxes_sound hboxes hi)
    (show (RationalInterval.point 8).Contains (8 : ℝ) by norm_num [RationalInterval.Contains, RationalInterval.point]) ends h).2

/-- Independent numerator/anchor uncertainty, including the cube root for the physical mass. -/
theorem normalized_sound {c m bits degree cb cd boxes R a ends depth J A result tol}
    (hboxes : boxes? c cb cd = some boxes) (hR : radiusBox? c m bits degree = some R)
    (hj : band? (weighted c) boxes R depth bits degree a ends = some J)
    (ha : band? (weighted c) boxes (.point 8) depth bits degree a ends = some A)
    (hn : SashimiNormalization.apply? c.sigma8 tol J A = some result) :
    result.Contains (normalize8 c.sigma8 (raw c a (ends.getLastD a) 8)
      (raw c a (ends.getLastD a) (radius c m))) ∧ result.hi-result.lo ≤ tol :=
  SashimiNormalization.apply_sound (raw_sound hboxes hR hj) (anchor_sound hboxes ha) hn

def variance? (c : Cosmology) (m a : ℚ) (ends : List ℚ) (depth bits degree : Nat)
    (tol : ℚ) : Option RationalInterval := do
  let boxes ← boxes? c bits degree
  let R ← radiusBox? c m bits degree
  let N ← band? (weighted c) boxes R depth bits degree a ends
  let A ← band? (weighted c) boxes (.point 8) depth bits degree a ends
  SashimiNormalization.apply? c.sigma8 tol N A

theorem variance_sound {c m a ends depth bits degree tol I}
    (h : variance? c m a ends depth bits degree tol = some I) :
    I.Contains (normalize8 c.sigma8 (raw c a (ends.getLastD a) 8)
      (raw c a (ends.getLastD a) (radius c m))) ∧ I.hi-I.lo ≤ tol := by
  obtain ⟨B,hB,h1⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨R,hR,h2⟩ := Option.bind_eq_some_iff.mp h1
  obtain ⟨J,hJ,h3⟩ := Option.bind_eq_some_iff.mp h2
  obtain ⟨A,hA,hN⟩ := Option.bind_eq_some_iff.mp h3
  exact normalized_sound hB hR hJ hA hN

theorem variance_contains_of_check {c m a ends depth bits degree tol L U}
    (h : ((variance? c m a ends depth bits degree tol).any
      fun I => decide (L ≤ I.lo ∧ I.hi ≤ U)) = true) :
    (⟨L,U⟩ : RationalInterval).Contains
      (normalize8 c.sigma8 (raw c a (ends.getLastD a) 8)
        (raw c a (ends.getLastD a) (radius c m))) := by
  cases hI : variance? c m a ends depth bits degree tol with
  | none => simp [hI] at h
  | some I =>
    simp [hI] at h
    have hh := (variance_sound hI).1
    exact ⟨(by exact_mod_cast h.1 : (L:ℝ) ≤ I.lo).trans hh.1,
      hh.2.trans (by exact_mod_cast h.2)⟩

#print axioms variance_sound

end SashimiEisensteinIntegral
