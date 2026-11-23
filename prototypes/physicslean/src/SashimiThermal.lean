import SashimiNumerics
import Mathlib.Analysis.Calculus.Deriv.Pow
import Mathlib.Tactic.Positivity

/-! Saikawa & Shirai, arXiv:1803.01038v2, Appendix C, Eqs. C.1–C.9.
Temperature is in GeV, and densities below use natural units. The analytic fits and
the published tabulated uncertainties are distinct inputs. A certificate for a fit
is unconditional mathematics; statements about the physical plasma name their hypotheses.
-/
namespace Cosmology.EffectiveDegreesOfFreedom

open FloatLib.Numerics FloatLib.Numerics.Interval Set

abbrev Temperature := ℝ
abbrev GRho := Temperature → ℝ
abbrev GS := Temperature → ℝ

inductive Quantity where | radiation | entropy deriving DecidableEq, Repr

private def poly (cs : List ℚ) (x : Expr) : Expr :=
  cs.foldr (fun c remainder => .const c + x * remainder) (.const 0)

-- Table 1, coefficients in increasing order, t = ln(T/GeV).
private def a : List ℚ := [1, 1.11724, 0.312672, -0.0468049, -0.0265004, -0.00119760,
  0.000182812, 0.000136436, 8.55051e-5, 1.22840e-5, 3.82259e-7, -6.87035e-9]
private def b : List ℚ := [0.0143382, 0.0137559, 0.00292108, -5.38533e-4, -1.62496e-4,
  -2.87906e-5, -3.84278e-6, 2.78776e-6, 7.40342e-7, 1.17210e-7, 3.72499e-9, -6.74107e-11]
private def c : List ℚ := [1, 0.607869, -0.154485, -0.224034, -0.0282147, 0.0290620,
  0.00686778, -0.00100005, -1.69104e-4, 1.06301e-5, 1.69528e-6, -9.33311e-8]
private def d : List ℚ := [70.7388, 91.8011, 33.1892, -1.39779, -1.52558, -0.0197857,
  -0.160146, 8.22615e-5, 0.0202651, -1.82134e-5, 7.83943e-5, 7.13518e-5]

private def high (q : Quantity) (T : Expr) : Expr :=
  let t := Expr.unary .log T
  let rho := Expr.binary .div (poly a t) (poly b t)
  match q with
  | .radiation => rho
  | .entropy => .binary .div rho (.const 1 + .binary .div (poly c t) (poly d t))

private def boltzmann (k x₁ x₂ x₃ : ℚ) (x : Expr) : Expr :=
  .unary .exp (.const (-k) * x) * poly [1,x₁,x₂,x₃] x

private def low (q : Quantity) (T : Expr) : Expr :=
  let x := fun (mass : ℚ) => Expr.binary .div (.const mass) T
  let se := .const 1 + .const (7/4) * boltzmann 1.0419 1.034 0.456426 0.0595249 (x 0.000511)
  let (offset, coeff, e, mu, bosons, fermion, boson) := match q with
    | .radiation => ((2.030:ℚ), (1.353:ℚ), (3.495:ℚ), (3.446:ℚ),
        [(1.05:ℚ), 2.08, 4.165, 30.55, 89.4, 8209],
        boltzmann 1.04855 1.03757 0.508630 0.0893988, boltzmann 1.03149 1.03317 0.398264 0.0648056)
    | .entropy => ((2.008:ℚ), (1.923:ℚ), (3.442:ℚ), (3.468:ℚ),
        [(1.034:ℚ), 2.068, 4.16, 30.55, 90, 6209],
        boltzmann 1.04190 1.03400 0.456426 0.0595248, boltzmann 1.03365 1.03397 0.342548 0.0506182)
  let neutrino := match q with
    | .radiation => Expr.unary .exp (.const (4/3) * .unary .log se)
    | .entropy => se
  let base := .const offset + .const coeff * neutrino + .const e * fermion (x 0.000511)
    + .const mu * fermion (x 0.1056)
  (bosons.zip [0.135, 0.140, 0.5, 0.77, 1.2, 2]).foldl
    (fun s (weight, mass) => s + .const weight * boson (x mass)) base

private def fitExpr (q : Quantity) (below : Bool) : Expr :=
  if below then low q (.var 0) else high q (.var 0)

def FitDomain (T : Temperature) : Prop := 0 < T ∧ T ≤ 10^16

noncomputable def fit (q : Quantity) (T : Temperature) : ℝ :=
  (fitExpr q (decide (T < ((0.12 : ℚ) : ℝ)))).eval (fun _ => T)

noncomputable def saikawaShiraiGRhoFit : GRho := fit .radiation
noncomputable def saikawaShiraiGSFit : GS := fit .entropy

/-- Enclose the Appendix C formula at an exact rational temperature.
    Out-of-domain inputs, backend failure and missed width budgets return none. -/
def encloseFit (q : Quantity) (T tol : ℚ) (bits : Nat := 64) (degree : Nat := 28) :
    Option RationalInterval :=
  if 0 < T ∧ T ≤ 10^16 then
    (SashimiNumerics.exprRange (fitExpr q (decide (T < 0.12))) bits degree T T).bind
      (SashimiNumerics.within? tol)
  else none

theorem encloseFit_sound {q T tol bits degree I} (h : encloseFit q T tol bits degree = some I) :
    FitDomain T ∧ I.Contains (fit q T) ∧ I.hi - I.lo ≤ tol := by
  unfold encloseFit at h
  split at h
  · rename_i valid
    obtain ⟨J, hJ, hwidth⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨rfl, hw⟩ := SashimiNumerics.within_sound hwidth
    have test : decide ((T : ℝ) < ((0.12 : ℚ) : ℝ)) = decide (T < (0.12 : ℚ)) := by
      simp only [Rat.cast_lt]
    refine ⟨⟨by exact_mod_cast valid.1, by exact_mod_cast valid.2⟩, ?_, hw⟩
    unfold fit; rw [test]
    exact SashimiNumerics.exprRange_sound _ _ _ _ _ _ hJ _ ⟨le_rfl, le_rfl⟩
  · contradiction

/-- Columns of standardmodel2018.dat: temperature, means and absolute error columns. -/
structure Knot where
  T : ℚ
  rho : ℚ
  rhoErr : ℚ
  entropy : ℚ
  entropyErr : ℚ
deriving Repr

structure Cell where
  left : Knot
  right : Knot
deriving Repr

def Knot.mean (p : Knot) (q : Quantity) := match q with
  | .radiation => p.rho | .entropy => p.entropy
def Knot.err (p : Knot) (q : Quantity) := match q with
  | .radiation => p.rhoErr | .entropy => p.entropyErr

def Cell.valid (p : Cell) : Prop := 0 < p.left.T ∧ p.left.T < p.right.T ∧
  0 ≤ p.left.rhoErr ∧ 0 ≤ p.right.rhoErr ∧ 0 ≤ p.left.entropyErr ∧ 0 ≤ p.right.entropyErr
instance (p : Cell) : Decidable p.valid := by unfold Cell.valid; infer_instance

/-- Chosen interpolation coordinate is T in GeV. Each cell has a one-sided slope at its endpoints. -/
def Cell.slope (p : Cell) (f : Knot → ℚ) : ℚ := (f p.right - f p.left) / (p.right.T - p.left.T)
def Cell.linear (p : Cell) (f : Knot → ℚ) (T : ℝ) : ℝ :=
  f p.left + (T - p.left.T) * (p.slope f : ℝ)

theorem Cell.linear_derivative (p : Cell) (f : Knot → ℚ) (T : ℝ) :
    HasDerivAt (p.linear f) (p.slope f : ℝ) T := by
  unfold Cell.linear
  simpa only [one_mul, id_eq] using (((hasDerivAt_id T).sub_const (p.left.T : ℝ)).mul_const
    (p.slope f : ℝ)).const_add (f p.left : ℝ)

def Cell.mean (p : Cell) (q : Quantity) := p.linear (·.mean q)
def Cell.error (p : Cell) (q : Quantity) := p.linear (·.err q)

private def lerpQ (p : Cell) (f : Knot → ℚ) (T : ℚ) : ℚ :=
  f p.left + (T - p.left.T) * p.slope f

/-- Exact rational band for a selected adjacent cell. Values outside the segment
    return none, including extrapolations past the first and last source knots. -/
def tableBand (p : Cell) (q : Quantity) (T : ℚ) : Option RationalInterval :=
  if p.valid ∧ p.left.T ≤ T ∧ T ≤ p.right.T then
    let center := lerpQ p (·.mean q) T
    let error := lerpQ p (·.err q) T
    some ⟨center - error, center + error⟩
  else none

/-- The tabulated uncertainty is centered on the tabulated mean. Recenter on the analytic
    fit by adding their absolute discrepancy, which can be nonzero also at tabulated nodes. -/
noncomputable def fitError (p : Cell) (q : Quantity) (T : Temperature) : ℝ :=
  p.error q T + |p.mean q T - fit q T|
noncomputable def saikawaShiraiGRhoError (p : Cell) := fitError p .radiation
noncomputable def saikawaShiraiGSError (p : Cell) := fitError p .entropy

/-- Empirical premise on a declared source cell: both physical functions lie in the
    interpolated bands at every temperature in that cell. This assumption is stronger
    than the endpoint data and has no assigned frequentist coverage. -/
structure SaikawaShiraiDoF (p : Cell) (gρ : GRho) (gs : GS) : Prop where
  valid : p.valid
  radiation : ∀ T ∈ Icc (p.left.T : ℝ) (p.right.T : ℝ), |gρ T - p.mean .radiation T| ≤ p.error .radiation T
  entropy : ∀ T ∈ Icc (p.left.T : ℝ) (p.right.T : ℝ), |gs T - p.mean .entropy T| ≤ p.error .entropy T

theorem SaikawaShiraiDoF.fit_bounds {p gρ gs} (h : SaikawaShiraiDoF p gρ gs) {T}
    (hT : T ∈ Icc (p.left.T : ℝ) (p.right.T : ℝ)) (_ : FitDomain T) :
    |gρ T - saikawaShiraiGRhoFit T| ≤ saikawaShiraiGRhoError p T ∧
    |gs T - saikawaShiraiGSFit T| ≤ saikawaShiraiGSError p T := by
  constructor
  · exact (abs_sub_le _ (p.mean .radiation T) _).trans (add_le_add (h.radiation T hT) le_rfl)
  · exact (abs_sub_le _ (p.mean .entropy T) _).trans (add_le_add (h.entropy T hT) le_rfl)

theorem tableBand_sound {p gρ gs} (h : SaikawaShiraiDoF p gρ gs) {q T I}
    (hr : tableBand p q T = some I) :
    I.Contains (match q with | .radiation => gρ T | .entropy => gs T) := by
  unfold tableBand at hr
  split at hr
  · rename_i hv
    cases hr
    have ht : (T : ℝ) ∈ Icc (p.left.T : ℝ) (p.right.T : ℝ) :=
      ⟨by exact_mod_cast hv.2.1, by exact_mod_cast hv.2.2⟩
    cases q with
    | radiation =>
      have he := abs_sub_le_iff.mp (h.radiation T ht)
      simp only [RationalInterval.Contains, lerpQ, Rat.cast_sub, Rat.cast_add, Rat.cast_mul]
      unfold Cell.error Cell.mean Cell.linear at he
      constructor <;> linarith [he.1, he.2]
    | entropy =>
      have he := abs_sub_le_iff.mp (h.entropy T ht)
      simp only [RationalInterval.Contains, lerpQ, Rat.cast_sub, Rat.cast_add, Rat.cast_mul]
      unfold Cell.error Cell.mean Cell.linear at he
      constructor <;> linarith [he.1, he.2]
  · contradiction

/-- Natural units: energy density in GeV⁴ and entropy density in GeV³. -/
noncomputable def radiationDensity (gρ : GRho) (T : Temperature) := (Real.pi^2 * T^4 / 30) * gρ T
noncomputable def entropyDensity (gs : GS) (T : Temperature) := (2 * Real.pi^2 * T^3 / 45) * gs T

theorem density_errors {p gρ gs} (h : SaikawaShiraiDoF p gρ gs) {T}
    (hT : T ∈ Icc (p.left.T : ℝ) (p.right.T : ℝ)) (hr : FitDomain T) :
    |radiationDensity gρ T - radiationDensity saikawaShiraiGRhoFit T| ≤
      (Real.pi^2 * T^4 / 30) * saikawaShiraiGRhoError p T ∧
    |entropyDensity gs T - entropyDensity saikawaShiraiGSFit T| ≤
      (2 * Real.pi^2 * T^3 / 45) * saikawaShiraiGSError p T := by
  have hpos := hr.1
  constructor
  · unfold radiationDensity; rw [← mul_sub, abs_mul, abs_of_nonneg (by positivity : 0 ≤ Real.pi^2 * T^4 / 30)]
    exact mul_le_mul_of_nonneg_left (h.fit_bounds hT hr).1 (by positivity)
  · unfold entropyDensity; rw [← mul_sub, abs_mul, abs_of_nonneg (by positivity : 0 ≤ 2 * Real.pi^2 * T^3 / 45)]
    exact mul_le_mul_of_nonneg_left (h.fit_bounds hT hr).2 (by positivity)

/-- From s(T) a³ conserved between two epochs; the scale-factor ratio cubed is useful
    directly for quantities proportional to a_kd⁻³ such as the prompt-cusp f_max. -/
theorem adiabatic_scale_cube {gs : GS} {a a₀ T T₀ : ℝ} (ha : a₀ ≠ 0) (ht : T ≠ 0)
    (hs : gs T ≠ 0) (h : entropyDensity gs T * a^3 = entropyDensity gs T₀ * a₀^3) :
    (a/a₀)^3 = (gs T₀ / gs T) * (T₀/T)^3 := by
  unfold entropyDensity at h
  field_simp [ha, ht, hs]
  linear_combination (norm := (field_simp [Real.pi_ne_zero]; ring)) (45 / (2 * Real.pi^2)) * h

/-- Radiation-dominated flat Friedmann contribution in GeV² for Newton's constant
    in GeV⁻². Nonnegative G preserves the model error bound. -/
noncomputable def radiationHubbleSquared (G : ℝ) (gρ : GRho) (T : Temperature) :=
  (8 * Real.pi * G / 3) * radiationDensity gρ T

theorem expansion_error {p gρ gs} (h : SaikawaShiraiDoF p gρ gs) {T G : ℝ}
    (hG : 0 ≤ G) (ht : T ∈ Icc (p.left.T : ℝ) (p.right.T : ℝ)) (hr : FitDomain T) :
    |radiationHubbleSquared G gρ T - radiationHubbleSquared G saikawaShiraiGRhoFit T| ≤
      (8 * Real.pi * G / 3) * ((Real.pi^2 * T^4 / 30) * saikawaShiraiGRhoError p T) := by
  unfold radiationHubbleSquared
  rw [← mul_sub, abs_mul, abs_of_nonneg (by positivity : 0 ≤ 8 * Real.pi * G / 3)]
  exact mul_le_mul_of_nonneg_left (density_errors h ht hr).1 (by positivity)

/-- The adiabatic radiation dilution relative to a⁻⁴, cubed to avoid fractional powers. -/
noncomputable def radiationTransferCube (gρs gs₀ gs gρ₀ : ℝ) := (gρs / gρ₀)^3 * (gs₀ / gs)^4

theorem radiation_transfer {gρ : GRho} {gs : GS} {a a₀ T T₀ : ℝ}
    (ha : a₀ ≠ 0) (ht : T ≠ 0) (ht₀ : T₀ ≠ 0) (hs : gs T ≠ 0) (hr₀ : gρ T₀ ≠ 0)
    (h : entropyDensity gs T * a^3 = entropyDensity gs T₀ * a₀^3) :
    ((radiationDensity gρ T / radiationDensity gρ T₀) * (a/a₀)^4)^3 =
      radiationTransferCube (gρ T) (gs T₀) (gs T) (gρ T₀) := by
  have hc := adiabatic_scale_cube ha ht hs h
  rw [mul_pow, ← pow_mul, Nat.mul_comm 4 3, pow_mul, hc]
  unfold radiationTransferCube radiationDensity
  field_simp [Real.pi_ne_zero, ht, ht₀, hr₀, hs]

private def transferExpr : Expr :=
  .unary (.pow 3) (.binary .div (.var 0) (.var 3)) *
    .unary (.pow 4) (.binary .div (.var 2) (.var 1))

/-- Propagate the two tabulated uncertainty bands at each of two epochs.
    The result is the cubed radiation dilution factor. This is the thermal
    background factor, not the halo matter-power transfer or a tensor-mode ODE. -/
def transferBand (p p₀ : Cell) (T T₀ tol : ℚ) (bits : Nat := 56) : Option RationalInterval := do
  let r ← tableBand p .radiation T
  let s ← tableBand p .entropy T
  let s₀ ← tableBand p₀ .entropy T₀
  let r₀ ← tableBand p₀ .radiation T₀
  let out ← transferExpr.gridEnclosureWith? [r,s,s₀,r₀] bits 0
  SashimiNumerics.within? tol out

theorem transferBand_sound {p p₀ T T₀ tol bits I gρ gs}
    (h : SaikawaShiraiDoF p gρ gs) (h₀ : SaikawaShiraiDoF p₀ gρ gs)
    (hr : transferBand p p₀ T T₀ tol bits = some I) :
    I.Contains (radiationTransferCube (gρ T) (gs T₀) (gs T) (gρ T₀)) ∧ I.hi-I.lo ≤ tol := by
  unfold transferBand at hr
  obtain ⟨r, hRad, tail⟩ := Option.bind_eq_some_iff.mp hr
  obtain ⟨s, hs, tail⟩ := Option.bind_eq_some_iff.mp tail
  obtain ⟨s₀, hs₀, tail⟩ := Option.bind_eq_some_iff.mp tail
  obtain ⟨r₀, hr₀, tail⟩ := Option.bind_eq_some_iff.mp tail
  obtain ⟨J, hj, hw⟩ := Option.bind_eq_some_iff.mp tail
  obtain ⟨rfl, hw⟩ := SashimiNumerics.within_sound hw
  let v : Nat → ℝ := fun | 0 => gρ T | 1 => gs T | 2 => gs T₀ | _ => gρ T₀
  refine ⟨Expr.contains_gridEnclosureWith? (values := v) (h := hj) ?_, hw⟩
  intro i A ha
  cases i with
  | zero => simp at ha; subst A; exact tableBand_sound h hRad
  | succ i => cases i with
    | zero => simp at ha; subst A; exact tableBand_sound h hs
    | succ i => cases i with
      | zero => simp at ha; subst A; exact tableBand_sound h₀ hs₀
      | succ i => cases i with
        | zero => simp at ha; subst A; exact tableBand_sound h₀ hr₀
        | succ i => simp at ha

/-- The energy/entropy value bands do not bound derivatives. A physical entropy slope
    must be supplied separately, e.g. for temperature-time or freeze-out calculations. -/
theorem entropy_slope {gs : GS} {T D : ℝ} (hg : HasDerivAt gs D T) :
    HasDerivAt (entropyDensity gs) ((2*Real.pi^2/45) * (3*T^2*gs T + T^3*D)) T := by
  unfold entropyDensity
  convert (((hasDerivAt_id T).pow 3).const_mul (2*Real.pi^2/45)).mul hg using 1 <;>
    first | (ext t; dsimp; ring) | (dsimp; ring)

/-- The actual cooling rate under comoving entropy conservation. Its input slope D is an
    additional physical hypothesis. A band for g_s at each point does not fix this slope. -/
theorem adiabatic_cooling {gs : GS} {a T : ℝ → ℝ} {t D H rate S : ℝ}
    (hD : HasDerivAt gs D (T t)) (hT : HasDerivAt T rate t) (ha : HasDerivAt a (H * a t) t)
    (hcons : ∀ u, entropyDensity gs (T u) * (a u)^3 = S)
    (hat : a t ≠ 0) (htt : T t ≠ 0) (hden : 3 * gs (T t) + T t * D ≠ 0) :
    rate = -(3 * H * T t * gs (T t)) / (3 * gs (T t) + T t * D) := by
  have law := ((entropy_slope hD).comp t hT).mul (ha.pow 3)
  change HasDerivAt (fun u => entropyDensity gs (T u) * (a u)^3) _ t at law
  simp only [hcons] at law
  have zero := law.unique (hasDerivAt_const t S)
  simp only [entropyDensity, Pi.pow_apply, Function.comp_apply] at zero
  apply (eq_div_iff hden).2
  linear_combination (norm := (field_simp [Real.pi_ne_zero, hat, htt]; ring))
    (45 / (2 * Real.pi^2 * (T t)^2 * (a t)^3)) * zero

/-- Decidable certificate for a concrete rounded enclosure of a fit. -/
theorem fit_at {q T tol bits degree lo hi} (checked :
    (encloseFit q T tol bits degree).any (fun I => decide (lo ≤ I.lo ∧ I.hi ≤ hi)) = true) :
    (lo : ℝ) ≤ fit q T ∧ fit q T ≤ (hi : ℝ) := by
  cases h : encloseFit q T tol bits degree with
  | none => simp [h] at checked
  | some I =>
    simp [h] at checked
    have hs := (encloseFit_sound h).2.1
    exact ⟨(by exact_mod_cast checked.1 : (lo : ℝ) ≤ I.lo).trans hs.1,
      hs.2.trans (by exact_mod_cast checked.2)⟩

-- The first pair samples the prompt-cusp kinetic decoupling temperature;
-- the second uses the high-temperature polynomial branch.
set_option maxRecDepth 16000 in
theorem fit_30MeV :
    12.5681 ≤ saikawaShiraiGRhoFit (3/100) ∧ saikawaShiraiGRhoFit (3/100) ≤ 12.5684 ∧
    12.3551 ≤ saikawaShiraiGSFit (3/100) ∧ saikawaShiraiGSFit (3/100) ≤ 12.3553 := by
  have r := fit_at (q := .radiation) (T := 3/100) (bits := 40) (degree := 20) (tol := 1/1000)
    (lo := 12.5681) (hi := 12.5684) (by decide +kernel)
  have s := fit_at (q := .entropy) (T := 3/100) (bits := 40) (degree := 20) (tol := 1/1000)
    (lo := 12.3551) (hi := 12.3553) (by decide +kernel)
  norm_num [saikawaShiraiGRhoFit, saikawaShiraiGSFit] at r s ⊢
  exact ⟨r.1, r.2, s⟩

theorem fit_1GeV :
    69.7437 ≤ saikawaShiraiGRhoFit 1 ∧ saikawaShiraiGRhoFit 1 ≤ 69.7438 ∧
    68.7715 ≤ saikawaShiraiGSFit 1 ∧ saikawaShiraiGSFit 1 ≤ 68.7716 := by
  have r := fit_at (q := .radiation) (T := 1) (bits := 32) (degree := 14) (tol := 1/1000)
    (lo := 69.7437) (hi := 69.7438) (by decide +kernel)
  have s := fit_at (q := .entropy) (T := 1) (bits := 32) (degree := 14) (tol := 1/1000)
    (lo := 68.7715) (hi := 68.7716) (by decide +kernel)
  norm_num [saikawaShiraiGRhoFit, saikawaShiraiGSFit] at r s ⊢
  exact ⟨r.1, r.2, s⟩

/-- Four source rows used in the transfer regression. Data-row numbers exclude the header. -/
def cell30MeV : Cell := -- 1841 and 1842
  ⟨⟨0.029992236, 12.568909, 0.24392569, 12.353211, 0.21288038⟩,
   ⟨0.03014942, 12.58898, 0.24730491, 12.371662, 0.21590458⟩⟩
def cell1GeV : Cell := -- 2511 and 2512
  ⟨⟨0.99537081, 69.716185, 5.4362231, 68.68934, 5.2698531⟩,
   ⟨1.0005874, 69.767232, 5.45202, 68.743712, 5.2843238⟩⟩

/-- The cubed radiation dilution from 1 GeV relative to 30 MeV, conditional on the
    two cells' empirical coverage. The range includes the uncertainty of both epochs. -/
theorem transfer_1GeV_to_30MeV {gρ gs} (hh : SaikawaShiraiDoF cell1GeV gρ gs)
    (hc : SaikawaShiraiDoF cell30MeV gρ gs) :
    0.091 ≤ radiationTransferCube (gρ 1) (gs (3/100)) (gs 1) (gρ (3/100)) ∧
    radiationTransferCube (gρ 1) (gs (3/100)) (gs 1) (gρ (3/100)) ≤ 0.350 := by
  have h := (transferBand_sound hh hc (T := 1) (T₀ := 3/100) (tol := 2) (bits := 32)
    (I := ⟨196501177/2147483648, 1501264839/4294967296⟩) (by decide +kernel)).1
  norm_num [RationalInterval.Contains] at h ⊢
  constructor <;> linarith [h.1,h.2]

example : transferBand cell1GeV cell30MeV 1 (3/100) (1/1000) 32 = none := by decide +kernel
example : tableBand cell30MeV .entropy 1 = none := by decide +kernel
example : encloseFit .entropy 0 1 = none ∧ encloseFit .radiation (10^16+1) 1 = none := by
  decide +kernel

#print axioms transfer_1GeV_to_30MeV
#print axioms fit_30MeV
#print axioms fit_1GeV
#print axioms adiabatic_cooling
#print axioms tableBand_sound
#print axioms transferBand_sound
#print axioms expansion_error
#print axioms encloseFit_sound
#print axioms SaikawaShiraiDoF.fit_bounds
#print axioms density_errors
#print axioms radiation_transfer
#print axioms entropy_slope

end Cosmology.EffectiveDegreesOfFreedom
