import SashimiIntegrable
import SashimiSpectrumVariance

/-! Certified use of a supplied linear-power table. Knots are exact rationals in
Mpc⁻¹ and Mpc³. Endpoints are clipped in k, so neither exponentiating an endpoint
nor sampling across a sharp cutoff changes the requested integration band. -/

namespace SashimiTable

open FloatLib.Numerics FloatLib.Numerics.Interval MeasureTheory

structure Knot where
  k : ℚ
  power : ℚ
deriving Repr, DecidableEq

inductive Interpolation where | linear | logLog deriving Repr, DecidableEq
inductive Placement where | continuous | atKnots deriving Repr, DecidableEq

structure Cell where
  left : Knot
  right : Knot
deriving Repr

/-- Convert a CAMB row in h/Mpc and (Mpc/h)³ to canonical physical units. -/
def fromHUnits (h : ℚ) (row : Knot) : Knot := ⟨row.k * h, row.power / h^3⟩

/-- Every cell spans adjacent data, including the two samples bracketing each query. -/
def cells (rows : List Knot) : List Cell :=
  (rows.zip rows.tail).map (fun (l,r) => ⟨l,r⟩)

def accepts (mode : Interpolation) (rows : List Knot) (a b : ℚ) : Bool :=
  match rows.head?, rows.getLast? with
  | some first, some last =>
    decide (first.k ≤ a ∧ a ≤ b ∧ b ≤ last.k ∧ rows.length ≥ 2) &&
      rows.all (fun r => decide (0 < r.k ∧ 0 ≤ r.power ∧ (mode = .logLog → 0 < r.power))) &&
      (cells rows).all (fun c => decide (c.left.k < c.right.k))
  | _, _ => false

private def div (x y : Expr) := Expr.binary .div x y
private def sub (x y : Expr) := Expr.binary .sub x y
private def log (x : Expr) := Expr.unary .log x
private def exp (x : Expr) := Expr.unary .exp x

/-- The SASHIMI-W interpolant is affine in k and P. The SASHIMI-F interpolant
is affine in log k and log P. A zero knot is legal in linear interpolation. -/
def interpExpr (mode : Interpolation) (l r yl yr x : Expr) : Expr :=
  match mode with
  | .linear => yl + div (sub x l) (sub r l) * sub yr yl
  | .logLog => exp (log yl + div (log (div x l)) (log (div r l)) * log (div yr yl))

noncomputable def interpReal (mode : Interpolation) (l r yl yr x : ℝ) : ℝ :=
  match mode with
  | .linear => yl + (x-l)/(r-l)*(yr-yl)
  | .logLog => Real.exp (Real.log yl + Real.log (x/l)/Real.log (r/l)*Real.log (yr/yl))

theorem eval_interp (mode : Interpolation) (l r p q x : Expr) (v : Nat → ℝ) :
    (interpExpr mode l r p q x).eval v =
      interpReal mode (l.eval v) (r.eval v) (p.eval v) (q.eval v) (x.eval v) := by
  cases mode <;> rfl

/-- Apply a power ratio continuously (e.g. WDM), or to the tabulated ordinates
before interpolation (e.g. log-interpolating a modified FDM spectrum). -/
def powerExpr (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (c : Cell) (x : Expr) : Expr :=
  let l := Expr.const c.left.k
  let r := Expr.const c.right.k
  let p := Expr.const c.left.power
  let q := Expr.const c.right.power
  match placement with
  | .continuous => interpExpr mode l r p q x * ratio x
  | .atKnots => interpExpr mode l r (p * ratio l) (q * ratio r) x

noncomputable def powerReal (mode : Interpolation) (placement : Placement) (ratio : ℝ → ℝ)
    (c : Cell) (x : ℝ) : ℝ :=
  let l : ℝ := c.left.k
  let r : ℝ := c.right.k
  let p : ℝ := c.left.power
  let q : ℝ := c.right.power
  match placement with
  | .continuous => interpReal mode l r p q x * ratio x
  | .atKnots => interpReal mode l r (p * ratio l) (q * ratio r) x

theorem eval_power (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (c : Cell) (x : Expr) (env : Nat → ℝ) (f : ℝ → ℝ)
    (hx : (ratio x).eval env = f (x.eval env))
    (hl : (ratio (.const c.left.k)).eval env = f c.left.k)
    (hr : (ratio (.const c.right.k)).eval env = f c.right.k) :
    (powerExpr mode placement ratio c x).eval env = powerReal mode placement f c (x.eval env) := by
  have mulE : ∀ (a b : Expr), (a*b).eval env = a.eval env * b.eval env := fun _ _ => rfl
  cases placement <;> simp only [powerExpr, powerReal, mulE, eval_interp, hx, hl, hr] <;> rfl

/-- Look up the specified table, without extrapolation. -/
def lookup? (mode : Interpolation) (placement : Placement) (ratio : Expr → Expr)
    (rows : List Knot) (k tol : ℚ) (bits degree : Nat) : Option RationalInterval := do
  guard (accepts mode rows k k)
  let c ← (cells rows).find? (fun c => decide (c.left.k ≤ k ∧ k ≤ c.right.k))
  SashimiNumerics.closed? (powerExpr mode placement ratio c (.const k)) bits tol degree

/-- k²P(k)W(k)²/(2π²) dk is the variance measure in linear wavenumber. -/
def moment (j : Nat) (mode : Interpolation) (placement : Placement) (ratio window : Expr → Expr)
    (c : Cell) : Expr :=
  let x := Expr.var 0
  let pi := Expr.const 4 * Expr.unary .atan (.const 1)
  div ((x.unary (.pow (2 + 2*j))) * powerExpr mode placement ratio c x *
    (window x).unary (.pow 2)) (Expr.const 2 * pi.unary (.pow 2))

/-- Real integrand in dk, with the same normalization as the log-k variance.
For j = 0 the change of measure is `SashimiSpectrum.windowVariance_dk`. -/
noncomputable def momentReal (j : Nat) (mode : Interpolation) (placement : Placement)
    (ratio window : ℝ → ℝ) (c : Cell) (k : ℝ) :=
  k ^ (2 + 2*j) * powerReal mode placement ratio c k * window k ^ 2 / (2 * Real.pi ^ 2)

theorem eval_moment (j : Nat) (mode : Interpolation) (placement : Placement)
    (ratio window : Expr → Expr) (f w : ℝ → ℝ) (c : Cell) (k : ℝ)
    (hx : (ratio (.var 0)).eval (fun _ => k) = f k)
    (hl : (ratio (.const c.left.k)).eval (fun _ => k) = f c.left.k)
    (hr : (ratio (.const c.right.k)).eval (fun _ => k) = f c.right.k)
    (hw : (window (.var 0)).eval (fun _ => k) = w k) :
    (moment j mode placement ratio window c).eval (fun _ => k) =
      momentReal j mode placement f w c k := by
  change k ^ (2 + 2*j) * (powerExpr mode placement ratio c (.var 0)).eval (fun _ => k) *
    ((window (.var 0)).eval (fun _ => k))^2 /
    ((2 : ℚ) * (SashimiSpectrum.piExpr.eval (fun _ => k)) ^ 2) = _
  rw [eval_power mode placement ratio c (.var 0) (fun _ => k) f hx hl hr, hw,
    SashimiSpectrum.eval_piExpr]
  norm_num [momentReal, Expr.eval]

/-- Integrate each adjacent cell only over its intersection with the requested band. -/
def piecewise? (f : Cell → Expr) (a b : ℚ) (depth bits degree : Nat) :
    List Cell → Option RationalInterval
  | [] => some (.point 0)
  | c :: cs => do
    let l := max a c.left.k
    let r := min b c.right.k
    if l < r then
      let box ← SashimiNumerics.rectangles? (SashimiNumerics.exprRange (f c) bits degree) l r depth
      return box.add (← piecewise? f a b depth bits degree cs)
    else piecewise? f a b depth bits degree cs

noncomputable def band (f : Cell → Expr) (a b : ℚ) : List Cell → ℝ
  | [] => 0
  | c :: cs =>
    let l := max a c.left.k
    let r := min b c.right.k
    if l < r then (∫ k in (l : ℝ)..(r : ℝ), (f c).eval fun _ => k) + band f a b cs
    else band f a b cs

def integrate? (f : Cell → Expr) (mode : Interpolation) (rows : List Knot)
    (a b tol : ℚ) (depth bits degree : Nat) : Option RationalInterval := do
  guard (accepts mode rows a b)
  let I ← piecewise? f a b depth bits degree (cells rows)
  SashimiNumerics.within? tol I

theorem lookup_sound {mode placement ratio rows k tol bits degree I}
    (h : lookup? mode placement ratio rows k tol bits degree = some I) :
    ∃ c, (cells rows).find? (fun c => decide (c.left.k ≤ k ∧ k ≤ c.right.k)) = some c ∧
      I.Contains ((powerExpr mode placement ratio c (.const k)).eval fun _ => 0) ∧
        I.hi - I.lo ≤ tol := by
  unfold lookup? at h
  obtain ⟨_, _, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨c, hc, h⟩ := Option.bind_eq_some_iff.mp h
  exact ⟨c, hc, SashimiNumerics.closed_sound h⟩

theorem piecewise_sound {f a b depth bits degree cs I}
    (h : piecewise? f a b depth bits degree cs = some I) : I.Contains (band f a b cs) := by
  induction cs generalizing I with
  | nil =>
    simp [piecewise?] at h
    subst I
    simp [band, RationalInterval.Contains, RationalInterval.point]
  | cons c cs ih =>
    simp only [piecewise?, band] at *
    split at h
    · rename_i hab
      rw [ite_eq_left hab]
      obtain ⟨A, hA, hrest⟩ := Option.bind_eq_some_iff.mp h
      obtain ⟨B, hB, hI⟩ := Option.bind_eq_some_iff.mp hrest
      cases hI
      have hr := SashimiNumerics.exprRange_sound (f c) bits degree
      have hf := SashimiNumerics.rectangles_integrable (SashimiNumerics.measurable_eval (f c))
        hr hab.le hA
      exact RationalInterval.contains_add (SashimiNumerics.rectangles_sound hr hab.le hf hA) (ih hB)
    · rename_i hab
      rw [ite_eq_right hab]
      exact ih h

/-- The successful numerical output encloses the sum of the actual clipped cell integrals.
The table ordinates are supplied data; every quadrature error is covered by the same width. -/
theorem integrate_sound {f mode rows a b tol depth bits degree I}
    (h : integrate? f mode rows a b tol depth bits degree = some I) :
    I.Contains (band f a b (cells rows)) ∧ I.hi - I.lo ≤ tol := by
  unfold integrate? at h
  obtain ⟨_, _, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨J, hJ, hj⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨rfl, hw⟩ := SashimiNumerics.within_sound hj
  exact ⟨piecewise_sound hJ, hw⟩

#print axioms integrate_sound

end SashimiTable
