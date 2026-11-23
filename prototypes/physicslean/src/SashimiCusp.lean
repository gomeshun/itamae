import SashimiTable
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts

/-! SASHIMI-C prompt_cusps.py uses a dimensionless, redshift-31 spectrum in Mpc⁻¹.
The ordinate is already Δ². Interpolation is log–log, the cutoff multiplies the
power once, and the peak moments use a finite band with no mass window. -/
namespace SashimiCusp
open SashimiTable FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory Real

/-- Piecewise exact redshift-31 integrand in dk. Its Gaussian acts on the power. -/
def integrand (j : Nat) (kfs : ℚ) (c : Cell) : Expr :=
  let k := Expr.var 0
  .binary .div ((k.unary (.pow (2*j))) *
    powerExpr .logLog .continuous (SashimiSpectrum.cutoffExpr (.const kfs)) c k) k

theorem integrand_eval (j : Nat) (kfs : ℚ) (c : Cell) (k : ℝ) :
    (integrand j kfs c).eval (fun _ => k) =
      k^(2*j) * (interpReal .logLog c.left.k c.right.k c.left.power c.right.power k *
        Real.exp (-(k/kfs)^2)) / k := rfl

/-- The finite logarithmic-band integral used by the reference code. -/
theorem log_measure (P : ℝ → ℝ) (j : Nat) {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    (∫ s in log a..log b, (exp s)^(2*j) * P (exp s)) =
      ∫ k in a..b, k^(2*j) * P k / k := by
  let g := fun k : ℝ => k^(2*j) * P k / k
  have h := intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
    (a := log a) (b := log b) (g := g) continuous_exp.continuousOn
    (fun s _ => hasDerivAt_exp s) (fun s _ => (exp_pos s).le)
  rw [exp_log ha,exp_log hb] at h
  apply Eq.trans _ h
  apply intervalIntegral.integral_congr
  intro x _
  dsimp [g]
  field_simp

/-- Adjacent data for each integral. Outside the table, SciPy interp1d's
'extrapolate' continues the first/last log–log slope. -/
inductive EdgePolicy where | inside | extrapolate deriving DecidableEq, Repr
structure Segment where
  lo : ℚ
  hi : ℚ
  cell : Cell

def segments (policy : EdgePolicy) (rows : List Knot) (a b : ℚ) : List Segment :=
  let cs := cells rows
  let inside := cs.map fun c => Segment.mk (max a c.left.k) (min b c.right.k) c
  match policy, cs.head?, cs.getLast? with
  | .extrapolate, some l, some r =>
    ⟨a, min b l.left.k, l⟩ :: inside ++ [⟨max a r.right.k, b, r⟩]
  | _, _, _ => inside

def valid (p : EdgePolicy) (rows : List Knot) (a b kfs : ℚ) : Bool :=
  match rows.head?, rows.getLast? with
  | some l, some r => accepts .logLog rows l.k r.k &&
    decide (0 < a ∧ a < b ∧ 0 < kfs ∧ (p = .inside → l.k ≤ a ∧ b ≤ r.k))
  | _, _ => false

noncomputable def integrateSegments (f : Cell → Expr) : List Segment → ℝ
  | [] => 0
  | p :: ps => (if p.lo < p.hi then
      ∫ k in (p.lo : ℝ)..(p.hi : ℝ), (f p.cell).eval fun _ => k else 0) + integrateSegments f ps

def segmentBoxes (f : Cell → Expr) (depth bits degree : Nat) :
    List Segment → Option RationalInterval
  | [] => some (.point 0)
  | p :: ps => do
    let tail ← segmentBoxes f depth bits degree ps
    if p.lo < p.hi then
      let box ← SashimiNumerics.rectangles? (SashimiNumerics.exprRange (f p.cell) bits degree)
        p.lo p.hi depth
      return box.add tail
    else return tail

theorem segmentBoxes_sound {f depth bits degree ps I}
    (h : segmentBoxes f depth bits degree ps = some I) : I.Contains (integrateSegments f ps) := by
  induction ps generalizing I with
  | nil =>
    simp [segmentBoxes] at h
    subst I
    simp [integrateSegments, RationalInterval.Contains, RationalInterval.point]
  | cons p ps ih =>
    obtain ⟨tail, htail, hrest⟩ := Option.bind_eq_some_iff.mp h
    by_cases hp : p.lo < p.hi
    · simp only [hp, ite_true] at hrest
      obtain ⟨B, hb, hI⟩ := Option.bind_eq_some_iff.mp hrest
      cases hI
      have box := SashimiNumerics.rectangles_sound (SashimiNumerics.exprRange_sound _ _ _) hp.le
        (SashimiNumerics.rectangles_integrable (SashimiNumerics.measurable_eval _) (SashimiNumerics.exprRange_sound _ _ _) hp.le hb) hb
      simpa only [integrateSegments, hp, ite_true] using RationalInterval.contains_add box (ih htail)
    · simp only [hp, ite_false, pure, Option.some.injEq] at hrest
      subst I
      simpa [integrateSegments, hp] using ih htail

/-- Evolution of dimensionless power from a₃₁=1/32. -/
noncomputable def growth (z : ℝ) := (32/(1+z)) ^ (901/500 : ℝ)

def growthExpr (z : ℚ) : Expr :=
  SashimiSpectrum.positivePowExpr (.const (32/(1+z))) (.const (901/500))

theorem growthExpr_eval {z : ℚ} (hz : -1 < z) (env : Nat → ℝ) :
    (growthExpr z).eval env = growth z := by
  have hp : 0 < (32 : ℚ)/(1+z) := div_pos (by decide) (by linarith)
  have h := SashimiSpectrum.eval_positivePowExpr (x := .const (32/(1+z))) (y := .const (901/500))
    (env := env) (show 0 < ((32/(1+z) : ℚ) : ℝ) by exact_mod_cast hp)
  simpa [growthExpr, growth, Expr.eval] using h

theorem growth_at_reference : growth 31 = 1 := by norm_num [growth]

inductive Quantity where | variance (order : Fin 3) | sigma (order : Fin 3) | radius | gamma
deriving DecidableEq, Repr

private def outputExpr (q : Quantity) (z : ℚ) : Expr :=
  let v := fun (i : Nat) => growthExpr z * Expr.var i
  let s := fun (i : Nat) => Expr.unary .sqrt (v i)
  match q with
  | .variance j => v j
  | .sigma j => s j
  | .radius => .binary .div (.unary .sqrt (.const 3) * s 1) (s 2)
  | .gamma => .binary .div ((s 1).unary (.pow 2)) (s 0 * s 2)

noncomputable def statistic (q : Quantity) (z : ℝ) (m : Nat → ℝ) : ℝ :=
  let v := fun i => growth z * m i
  let s := fun i => Real.sqrt (v i)
  match q with
  | .variance j => v j
  | .sigma j => s j
  | .radius => Real.sqrt 3 * s 1 / s 2
  | .gamma => (s 1)^2 / (s 0*s 2)

/-- The scale-independent growth of this small-scale prescription cancels in
the two BBKS spectral-shape parameters on a fixed band with a fixed cutoff. -/
theorem redshift_shape {z : ℝ} (hz : -1 < z) (m : Nat → ℝ) :
    statistic .radius z m = statistic .radius 31 m ∧
      statistic .gamma z m = statistic .gamma 31 m := by
  have hg : 0 < growth z := Real.rpow_pos_of_pos (div_pos (by norm_num) (by linarith)) _
  have hn := (Real.sqrt_pos.2 hg).ne'
  simp only [statistic, growth_at_reference, one_mul, Real.sqrt_mul hg.le, mul_pow]
  constructor <;> field_simp [hn] <;> ring

theorem redshift_sigma {z : ℝ} (hz : -1 < z) (m : Nat → ℝ) (j : Fin 3) :
    statistic (.sigma j) z m = Real.sqrt (growth z) * statistic (.sigma j) 31 m := by
  have hg : 0 < growth z := Real.rpow_pos_of_pos (div_pos (by norm_num) (by linarith)) _
  simp only [statistic, growth_at_reference, one_mul, Real.sqrt_mul hg.le]

private theorem outputExpr_eval {z : ℚ} (hz : -1 < z) (q : Quantity) (env : Nat → ℝ) :
    (outputExpr q z).eval env = statistic q z env := by
  have h := growthExpr_eval hz env
  cases q <;> simp [statistic, outputExpr, Expr.eval, BinaryOp.eval, UnaryOp.eval, h]

/-- Integrate the three powers of k once. Reuse the same boxes for each
redshift and derived observable, which have different physical units. -/
def moments? (p : EdgePolicy) (rows : List Knot) (a b kfs : ℚ) (depth bits degree : Nat) :
    Option (List RationalInterval) := do
  guard (valid p rows a b kfs)
  let box := fun j => segmentBoxes (integrand j kfs) depth bits degree (segments p rows a b)
  let a₀ ← box 0
  let a₁ ← box 1
  let a₂ ← box 2
  return [a₀, a₁, a₂]

noncomputable def moment (p : EdgePolicy) (rows : List Knot) (a b kfs : ℚ) (j : Nat) :=
  integrateSegments (integrand j kfs) (segments p rows a b)

def peak? (q : Quantity) (z tol : ℚ) (m : List RationalInterval) (bits degree : Nat) :
    Option RationalInterval := do
  guard (-1 < z ∧ m.length = 3)
  let box ← (outputExpr q z).gridEnclosureWith? m bits degree
  SashimiNumerics.within? tol box

theorem peak_sound {p rows a b kfs depth bits degree m q z tol bits' degree' I}
    (hm : moments? p rows a b kfs depth bits degree = some m)
    (h : peak? q z tol m bits' degree' = some I) :
    I.Contains (statistic q z (moment p rows a b kfs)) ∧ I.hi - I.lo ≤ tol := by
  unfold moments? at hm
  obtain ⟨_, _, hm⟩ := Option.bind_eq_some_iff.mp hm
  obtain ⟨A, hA, hm⟩ := Option.bind_eq_some_iff.mp hm
  obtain ⟨B, hB, hm⟩ := Option.bind_eq_some_iff.mp hm
  obtain ⟨C, hC, hm⟩ := Option.bind_eq_some_iff.mp hm
  cases hm
  unfold peak? at h
  obtain ⟨_, hg, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨J, hJ, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨rfl, hw⟩ := SashimiNumerics.within_sound h
  have hz : -1 < z := (show -1 < z ∧ [A,B,C].length = 3 by simpa [guard] using hg).1
  rw [← outputExpr_eval hz]
  refine ⟨Expr.contains_gridEnclosureWith? (h := hJ) ?_, hw⟩
  intro i box hx
  cases i with
  | zero => simp at hx; subst box; exact segmentBoxes_sound hA
  | succ i => cases i with
    | zero => simp at hx; subst box; exact segmentBoxes_sound hB
    | succ i => cases i with
      | zero => simp at hx; subst box; exact segmentBoxes_sound hC
      | succ i => simp at hx

#print axioms peak_sound
#print axioms log_measure

end SashimiCusp
