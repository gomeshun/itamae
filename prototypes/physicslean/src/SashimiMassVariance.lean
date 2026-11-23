import SashimiTable
import SashimiSpectrumSlopeCore

/-! Variance at a physical mass on a supplied finite power-spectrum band.
The sharp-k endpoint is a real expression including π and a cube root.
Each clipped cell is reparameterized on the rational unit interval; no
rounded endpoint replaces the real cutoff. All powers and radii use Mpc,
solar masses, and physical comoving wavenumber in inverse Mpc. -/
namespace SashimiMassVariance
open SashimiTable SashimiSpectrum FloatLib.Numerics FloatLib.Numerics.Interval
open Real Set MeasureTheory

inductive Calibration where
  | fixed (c : ℚ)
  | sashimiW
deriving DecidableEq

/-- Legacy family cutoff coefficient; includes the extra factor (9π/2)^(1/3). -/
abbrev Calibration.wFamilyLegacy : Calibration := .sashimiW

/-- Public W main uses the Schneider coefficient c=2.5 in k_max=c/R_th. -/
abbrev Calibration.wMain : Calibration := .fixed (5/2)

def Calibration.expression : Calibration → Expr
  | .fixed c => .const c
  | .sashimiW => .const (5/2) * positivePowExpr (.const (9/2) * piExpr) (.const (1/3))

noncomputable def Calibration.value : Calibration → ℝ
  | .fixed c => c
  | .sashimiW => (5/2) * (9*π/2)^(1/3 : ℝ)

theorem Calibration.eval (c : Calibration) (v : Nat → ℝ) : c.expression.eval v = c.value := by
  cases c with
  | fixed => rfl
  | sashimiW =>
    have h : 0 < (Expr.const (9/2) * piExpr).eval v := by
      change 0 < ((9/2 : ℚ) : ℝ) * piExpr.eval v; rw [eval_piExpr]; positivity
    change ((5/2 : ℚ) : ℝ) * _ = _
    rw [eval_positivePowExpr h]
    change ((5/2 : ℚ) : ℝ) * (((9/2 : ℚ) : ℝ) * piExpr.eval v)^(((1/3 : ℚ) : ℝ)) = _
    rw [eval_piExpr]
    norm_num [Calibration.value]; ring_nf

def radiusExpr (rho mass : ℚ) : Expr :=
  positivePowExpr (.binary .div (.const (3*mass)) (.const (4*rho) * piExpr)) (.const (1/3))

theorem radius_eval (rho mass : ℚ) (hr : 0 < rho) (hm : 0 < mass) (v : Nat → ℝ) :
    (radiusExpr rho mass).eval v = massRadius rho mass := by
  have hr' : 0 < (rho:ℝ) := by exact_mod_cast hr
  have hm' : 0 < (mass:ℝ) := by exact_mod_cast hm
  have he : (Expr.binary .div (.const (3*mass)) (.const (4*rho) * piExpr)).eval v =
      3*(mass:ℝ)/(4*π*(rho:ℝ)) := by
    change ((3*mass:ℚ):ℝ)/(((4*rho:ℚ):ℝ)*piExpr.eval v) = _
    rw [eval_piExpr]; push_cast; ring
  unfold radiusExpr massRadius
  have hp : 0 < (Expr.binary .div (.const (3*mass)) (.const (4*rho) * piExpr)).eval v := by
    rw [he]; positivity
  rw [eval_positivePowExpr hp, he]
  norm_num [Expr.eval]

def cutoff (c : Calibration) (rho mass : ℚ) : Expr :=
  .binary .div c.expression (radiusExpr rho mass)

theorem cutoff_eval (c : Calibration) {rho mass : ℚ} (hr : 0 < rho) (hm : 0 < mass) (v : Nat → ℝ) :
    (cutoff c rho mass).eval v = c.value / massRadius rho mass := by
  change c.expression.eval v / (radiusExpr rho mass).eval v = _
  rw [c.eval, radius_eval _ _ hr hm]

/-- The integrand for dS/dM at fixed growth D, including the signed W'. -/
def topHatPrimeExpr (x : Expr) : Expr :=
  .binary .add (.binary .sub (.binary .div (.const 9 * .unary .cos x) (.unary (.pow 3) x))
      (.binary .div (.const 9 * .unary .sin x) (.unary (.pow 4) x)))
    (.binary .div (.const 3 * .unary .sin x) (.unary (.pow 2) x))

theorem eval_prime (x : Expr) (v : Nat → ℝ) : (topHatPrimeExpr x).eval v = topHatPrime (x.eval v) := by
  simp [topHatPrimeExpr, topHatPrime, Expr.eval, UnaryOp.eval, BinaryOp.eval]

inductive Quantity where | sharp | topHat | topHatDerivative deriving DecidableEq

def integrand (q : Quantity) (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (rho mass growth : ℚ) (cell : Cell) (k : Expr) : Expr :=
  let base := .binary .div (.const (growth^2) * k.unary (.pow 2) * powerExpr mode placement ratio cell k)
    (.const 2 * piExpr.unary (.pow 2))
  let x := k * radiusExpr rho mass
  match q with
  | .sharp => base
  | .topHat => base * (topHatExpr x).unary (.pow 2)
  | .topHatDerivative => .binary .div
      (base * .const 2 * topHatExpr x * topHatPrimeExpr x * x) (.const (3*mass))

noncomputable def kernel (q : Quantity) (P : ℝ → ℝ) (rho mass growth k : ℝ) :=
  let base := growth^2 * k^2 * P k / (2*π^2)
  let x := k * massRadius rho mass
  match q with
  | .sharp => base
  | .topHat => base * topHat x ^ 2
  | .topHatDerivative => base * 2 * topHat x * topHatPrime x * x / (3*mass)

theorem eval_integrand (q mode placement ratio cell k) (rho mass growth : ℚ)
    (hr : 0 < rho) (hm : 0 < mass) (v : Nat → ℝ) (hk : 0 < (k : Expr).eval v) :
    (integrand q mode placement ratio rho mass growth cell k).eval v =
      kernel q (fun _ => (powerExpr mode placement ratio cell k).eval v) rho mass growth (k.eval v) := by
  have hR : 0 < massRadius rho mass := by
    have rr : 0 < (rho:ℝ) := by exact_mod_cast hr
    have mm : 0 < (mass:ℝ) := by exact_mod_cast hm
    unfold massRadius; positivity
  have hx : (k * radiusExpr rho mass).eval v ≠ 0 := by
    change k.eval v * (radiusExpr rho mass).eval v ≠ 0
    rw [radius_eval _ _ hr hm]; exact (mul_pos hk hR).ne'
  cases q <;> simp only [integrand, kernel, Expr.eval, UnaryOp.eval, BinaryOp.eval,
    eval_piExpr, eval_topHatExpr hx, eval_prime, radius_eval _ _ hr hm] <;> push_cast <;> ring

/-- Real change of variables for a clipped cell; covers collapsed cells too. -/
theorem unit_integral (g : ℝ → ℝ) {l u : ℝ} (h : l ≤ u) :
    (∫ t in (0:ℝ)..1, g (l + t*(u-l)) * (u-l)) = ∫ k in l..u, g k := by
  let f := fun t : ℝ => l+t*(u-l)
  have hf (t : ℝ) : HasDerivAt f (u-l) t := by
    simpa [f] using ((hasDerivAt_id' t).mul_const (u-l)).const_add l
  have ht := intervalIntegral.integral_comp_mul_deriv_of_deriv_nonneg
    (a := (0:ℝ)) (b := 1) (f := f) (g := g) (f' := fun _ => u-l)
    (by dsimp [f]; fun_prop) (fun t _ => hf t) (fun _ _ => sub_nonneg.mpr h)
  simpa only [f, Function.comp_apply, zero_mul, add_zero, one_mul, add_sub_cancel] using ht

def clip (q : Quantity) (K : Expr) (l r : ℚ) : Expr :=
  if q = .sharp then .binary .max (.const l) (.binary .min (.const r) K) else .const r

noncomputable def upper (q : Quantity) (K l r : ℝ) := if q = .sharp then max l (min r K) else r

def unitExpr (g : Expr → Expr) (q : Quantity) (K : Expr) (l r : ℚ) : Expr :=
  let d := Expr.binary .sub (clip q K l r) (.const l)
  g (.const l + .var 0 * d) * d

/-- A finite sum over the adjacent real cells, physically clipped before integration. -/
noncomputable def massBand (q : Quantity) (g : Cell → ℝ → ℝ) (K : ℝ) (a b : ℚ) : List Cell → ℝ
  | [] => 0
  | cell::cs =>
    let l := max a cell.left.k
    let r := min b cell.right.k
    (if l < r then ∫ k in (l:ℝ)..upper q K l r, g cell k else 0) + massBand q g K a b cs

def panels? (q : Quantity) (g : Cell → Expr → Expr) (K : Expr) (a b : ℚ) (depth bits degree : Nat) :
    List Cell → Option RationalInterval
  | [] => some (.point 0)
  | cell::cs => do
    let l := max a cell.left.k
    let r := min b cell.right.k
    if l < r then
      let e := unitExpr (g cell) q K l r
      let I ← SashimiNumerics.rectangles? (SashimiNumerics.exprRange e bits degree) 0 1 depth
      return I.add (← panels? q g K a b depth bits degree cs)
    else panels? q g K a b depth bits degree cs

theorem panels_sound {q g K} {a b : ℚ} {depth bits degree cs I} (G : Cell → ℝ → ℝ) (kc : ℝ)
    (hcut : ∀ v, (K : Expr).eval v = kc)
    (hg : ∀ cell ∈ cs, ∀ x : Expr, ∀ v, (a : ℝ) ≤ x.eval v →
      (g cell x).eval v = G cell (x.eval v))
    (h : panels? q g K a b depth bits degree cs = some I) : I.Contains (massBand q G kc a b cs) := by
  induction cs generalizing I with
  | nil => cases h; simp [massBand, RationalInterval.Contains, RationalInterval.point]
  | cons cell cs ih =>
    let l := max a cell.left.k
    let r := min b cell.right.k
    simp only [panels?] at h
    split at h
    · rename_i hlr
      have hlr' : (l:ℝ) ≤ (r:ℝ) := by exact_mod_cast hlr.le
      let u := upper q kc l r
      have hlu : (l:ℝ) ≤ u := by dsimp [u, upper]; split <;> first | exact le_max_left _ _ | exact hlr'
      have hclip (v : ℕ → ℝ) : (clip q K l r).eval v = u := by
        dsimp [u, clip, upper]; split <;> simp only [Expr.eval, BinaryOp.eval, hcut]
      obtain ⟨A,hA,h'⟩ := Option.bind_eq_some_iff.mp h
      obtain ⟨R,hR,hI⟩ := Option.bind_eq_some_iff.mp h'
      cases hI
      have hh := SashimiNumerics.exprRange_sound (unitExpr (g cell) q K l r) bits degree
      have ham := SashimiNumerics.rectangles_integrable (SashimiNumerics.measurable_eval _) hh (by norm_num : (0:ℚ) ≤ 1) hA
      have hac := SashimiNumerics.rectangles_sound hh (by norm_num : (0:ℚ) ≤ 1) ham hA
      have heq : (∫ t in (0:ℝ)..1, (unitExpr (g cell) q K l r).eval fun _ => t) =
          ∫ k in (l:ℝ)..u, G cell k := by
        rw [← unit_integral (G cell) hlu]
        apply intervalIntegral.integral_congr
        intro t ht
        simp only [uIcc_of_le (zero_le_one : (0:ℝ) ≤ 1)] at ht
        dsimp only [unitExpr, Expr.eval, BinaryOp.eval]
        rw [hg cell (by simp), hclip]
        · change G cell ((l:ℝ) + t*((clip q K l r).eval _-(l:ℝ))) * _ = _
          rw [hclip]
        · change (a : ℝ) ≤ (l:ℝ) + t*((clip q K l r).eval _ - (l:ℝ))
          rw [hclip]
          exact (by exact_mod_cast le_max_left a cell.left.k : (a:ℝ) ≤ (l:ℝ)).trans
            (le_add_of_nonneg_right (mul_nonneg ht.1 (sub_nonneg.mpr hlu)))
      norm_num only [Rat.cast_zero, Rat.cast_one] at hac
      rw [heq] at hac
      simp only [massBand, ite_eq_left hlr]
      exact RationalInterval.contains_add hac (ih (fun c hc => hg c (by simp [hc])) hR)
    · rename_i hlr
      simp only [massBand, ite_eq_right hlr, zero_add]
      exact ih (fun c hc => hg c (by simp [hc])) h

def evaluate? (q : Quantity) (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (cal : Calibration) (rho mass growth : ℚ) (rows : List Knot) (a b tol : ℚ)
    (depth bits degree : Nat) : Option RationalInterval := do
  guard (0 < a ∧ 0 < rho ∧ 0 < mass ∧ 0 ≤ growth)
  guard (accepts mode rows a b)
  let I ← panels? q (integrand q mode placement ratio rho mass growth)
    (cutoff cal rho mass) a b depth bits degree (cells rows)
  SashimiNumerics.within? tol I

/-- A pointwise refinement of the power interpolant is enough; no numerical
integral is assumed. The accepted tolerance includes radius and endpoint error. -/
theorem evaluate_sound {q mode placement ratio cal rho mass growth rows a b tol depth bits degree I}
    (P : Cell → ℝ → ℝ)
    (hp : ∀ cell ∈ cells rows, ∀ x : Expr, ∀ v, (0:ℝ) < x.eval v →
      (powerExpr mode placement ratio cell x).eval v = P cell (x.eval v))
    (h : evaluate? q mode placement ratio cal rho mass growth rows a b tol depth bits degree = some I) :
    I.Contains (massBand q (fun c k => kernel q (P c) rho mass growth k)
      (cal.value / massRadius rho mass) a b (cells rows)) ∧ I.hi - I.lo ≤ tol := by
  unfold evaluate? at h
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  have hd : 0 < a ∧ 0 < rho ∧ 0 < mass ∧ 0 ≤ growth := by by_contra hn; simp [guard, hn] at hg
  obtain ⟨_,_,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨J,hJ,hw⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨rfl,hw⟩ := SashimiNumerics.within_sound hw
  refine ⟨panels_sound _ _ (cutoff_eval cal hd.2.1 hd.2.2.1) ?_ hJ, hw⟩
  intro c hc x v hax
  have hx := (by exact_mod_cast hd.1 : (0:ℝ) < a).trans_le hax
  rw [eval_integrand _ _ _ _ _ _ _ _ _ hd.2.1 hd.2.2.1 v hx, hp c hc x v hx]
  cases q <;> rfl

/-- Within a strictly interior table cell the moving sharp cutoff contributes
the fundamental-theorem boundary term. A bracket touching a band edge or
crossing a knot is rejected by this implementation. -/
def sharpDerivative? (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (cal : Calibration) (rho mass growth : ℚ) (rows : List Knot) (a b tol : ℚ)
    (bits degree : Nat) : Option RationalInterval := do
  guard (0 < a ∧ 0 < rho ∧ 0 < mass ∧ 0 ≤ growth)
  guard (accepts mode rows a b)
  let K := cutoff cal rho mass
  let B ← K.gridEnclosure? bits degree
  if B.hi < a ∨ b < B.lo then
    SashimiNumerics.within? tol (.point 0)
  else
    guard (a < B.lo ∧ B.hi < b)
    let c ← (cells rows).find? (fun c => decide (c.left.k < B.lo ∧ B.hi < c.right.k))
    let e := Expr.binary .div (.const (-(growth^2)) * K.unary (.pow 3) * powerExpr mode placement ratio c K)
      (.const (6*mass) * piExpr.unary (.pow 2))
    SashimiNumerics.closed? e bits tol degree

theorem sharpDerivative_sound {mode placement ratio cal rho mass growth rows a b tol bits degree I}
    (P : Cell → ℝ → ℝ)
    (hp : ∀ cell ∈ cells rows, ∀ x : Expr, ∀ v, (0:ℝ) < x.eval v →
      (powerExpr mode placement ratio cell x).eval v = P cell (x.eval v))
    (h : sharpDerivative? mode placement ratio cal rho mass growth rows a b tol bits degree = some I) :
    let kc := cal.value / massRadius rho mass
    I.hi - I.lo ≤ tol ∧
      (((kc < (a:ℝ) ∨ (b:ℝ) < kc) ∧ I.Contains 0) ∨
      ((a:ℝ) < kc ∧ kc < (b:ℝ) ∧ ∃ c ∈ cells rows, (c.left.k:ℝ) < kc ∧ kc < (c.right.k:ℝ) ∧
        I.Contains (-((growth:ℝ)^2)*kc^3*P c kc / (6*(mass:ℝ)*π^2)))) := by
  dsimp only
  unfold sharpDerivative? at h
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  have hd : 0 < a ∧ 0 < rho ∧ 0 < mass ∧ 0 ≤ growth := by by_contra hn; simp [guard, hn] at hg
  obtain ⟨_,_,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨B,hB,h⟩ := Option.bind_eq_some_iff.mp h
  have hb := Expr.contains_gridEnclosure? hB
  rw [cutoff_eval _ hd.2.1 hd.2.2.1] at hb
  split at h
  · rename_i hout
    obtain ⟨rfl,hw⟩ := SashimiNumerics.within_sound h
    refine ⟨hw,Or.inl ⟨?_, by simp [RationalInterval.Contains, RationalInterval.point]⟩⟩
    rcases hout with hl|hr
    · exact Or.inl (hb.2.trans_lt (by exact_mod_cast hl))
    · exact Or.inr ((by exact_mod_cast hr : (b:ℝ) < B.lo).trans_le hb.1)
  · obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
    have hi : a < B.lo ∧ B.hi < b := by by_contra hn; simp [guard, hn] at hg
    obtain ⟨c,hc,h⟩ := Option.bind_eq_some_iff.mp h
    have mem := List.mem_of_find?_eq_some hc
    have hrow : c.left.k < B.lo ∧ B.hi < c.right.k := by
      simpa only [decide_eq_true_eq] using (List.find?_some (p := fun c : Cell => decide (c.left.k < B.lo ∧ B.hi < c.right.k)) hc)
    have hlo := (by exact_mod_cast hi.1 : (a:ℝ) < B.lo).trans_le hb.1
    have hhi := hb.2.trans_lt (by exact_mod_cast hi.2 : (B.hi:ℝ) < b)
    obtain ⟨heq, hw⟩ := SashimiNumerics.closed_sound h
    refine ⟨hw,Or.inr ⟨hlo,hhi,c,mem,(by exact_mod_cast hrow.1 : (c.left.k:ℝ) < B.lo).trans_le hb.1,
      hb.2.trans_lt (by exact_mod_cast hrow.2),?_⟩⟩
    have hk : 0 < (cutoff cal rho mass).eval (fun _ => 0) := by
      rw [cutoff_eval _ hd.2.1 hd.2.2.1]
      exact (by exact_mod_cast hd.1 : (0:ℝ) < a).trans hlo
    simpa only [Expr.eval, UnaryOp.eval, BinaryOp.eval, hp c mem _ _ hk,
      cutoff_eval _ hd.2.1 hd.2.2.1, eval_piExpr, Rat.cast_neg, Rat.cast_pow,
      Rat.cast_mul, Rat.cast_ofNat] using heq

#print axioms evaluate_sound
#print axioms sharpDerivative_sound

end SashimiMassVariance
