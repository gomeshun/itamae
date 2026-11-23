import SashimiHalo
import SashimiErf

/-! # SASHIMI-W Ludlow concentration grid
The input variance difference is S(0.02 M)-S(M) for the selected mass bin of the W
top-hat table with its σ₈ normalization. This is the discrete concentration calculation,
using 100 trial values 10^(4i/99), the source exponent 0.3333 and A=650/200. -/
namespace SashimiHaloWarm

open FloatLib.Numerics FloatLib.Numerics.Interval

local instance (n : Nat) : OfNat Expr n := ⟨.const n⟩
local instance : Sub Expr := ⟨.binary .sub⟩
local instance : Div Expr := ⟨.binary .div⟩
local instance : HPow Expr Nat Expr := ⟨fun x n => .unary (.pow n) x⟩

def trial (i : Nat) : Expr := SashimiHalo.posPower 10 (.const (4 * (i : ℚ) / 99))
def massFraction (c : Expr) : Expr := SashimiHalo.fc 1 / SashimiHalo.fc c

/-- Base of (1+z₂)=base^0.3333. A negative value is removed by the legacy finite mask.
At zero the growth expression encounters z₂=-1. -/
def formationBase (Ω z : ℚ) (c m2 : Expr) : Expr :=
  let Λ := 1-Ω
  let e2 := Ω * (1+z)^3 + Λ
  ((c^3 * m2 / .const (650/200)) * .const e2 - .const Λ) / .const Ω

private def twoInputs (x y : ℝ) : Nat → ℝ := fun i => if i = 0 then x else y

private theorem gridOne {e : Expr} {A : RationalInterval} {x : ℝ} {bits degree J}
    (hx : A.Contains x) (h : e.gridEnclosureWith? [A] bits degree = some J) :
    J.Contains (e.eval fun _ => x) := by
  apply Expr.contains_gridEnclosureWith? (h := h)
  intro i K hK
  cases i with
  | zero => simp at hK; subst K; exact hx
  | succ i => simp at hK

private theorem gridTwo {e : Expr} {A B : RationalInterval} {x y : ℝ} {bits degree J}
    (hx : A.Contains x) (hy : B.Contains y) (h : e.gridEnclosureWith? [A,B] bits degree = some J) :
    J.Contains (e.eval (twoInputs x y)) := by
  apply Expr.contains_gridEnclosureWith? (h := h)
  intro i K hK
  rcases i with _|_|i
  · simp at hK; subst K; exact hx
  · simp at hK; subst K; exact hy
  · simp at hK

/-- The barrier 1.686 / D(z), with the W concentration routine's internal *flat*
ΩΛ=1−Ωm convention and normalized CPT growth. Input u = 1+z. -/
def barrier (Ω : ℚ) (u : Expr) : Expr :=
  let l : Expr := .const (1-Ω)
  let g := fun x y : Expr =>
    (.const (5/2) * x) / (SashimiHalo.posPower x (.const (4/7)) - y + (1 + x/2) * (1 + y/70))
  let y := l / (l + .const Ω * u^3)
  .const 1.686 * g (.const Ω) l * u / g (1-y) y

def erfArg (Ω z : ℚ) (base ΔS : Expr) : Expr :=
  let u2 := SashimiHalo.posPower base (.const 0.3333)
  (barrier Ω u2 - barrier Ω (.const (1+z))) / .unary .sqrt (2*ΔS)

structure Sample (α : Type) where
  concentration : α
  residual : α
deriving Repr

def Sample.Covers (I : Sample RationalInterval) (x : Sample ℝ) : Prop :=
  I.concentration.Contains x.concentration ∧ I.residual.Contains x.residual

noncomputable def realSample (erf : ℝ → ℝ) (Ω z : ℚ) (dS : ℝ) (i : Nat) : Option (Sample ℝ) := Id.run do
  let c := (trial i).eval (fun _ => 0)
  let m := (massFraction (.var 0)).eval (fun _ => c)
  let b := (formationBase Ω z (.var 0) (.var 1)).eval (twoInputs c m)
  if b ≤ 0 then return none
  let v := (erfArg Ω z (.var 0) (.var 1)).eval (twoInputs b dS)
  return some ⟨c, m + erf v - 1⟩

/-- Outer none means undecided. Inner none is a provably invalid point of the source grid. -/
def sample? (erf? : RationalInterval → Option RationalInterval) (Ω z : ℚ) (dS : RationalInterval)
    (bits degree i : Nat) : Option (Option (Sample RationalInterval)) := do
  let C ← (trial i).gridEnclosure? bits degree
  let M ← (massFraction (.var 0)).gridEnclosureWith? [C] bits degree
  let B ← (formationBase Ω z (.var 0) (.var 1)).gridEnclosureWith? [C,M] bits degree
  if B.hi ≤ 0 then return none
  guard (0 < B.lo)
  let X ← (erfArg Ω z (.var 0) (.var 1)).gridEnclosureWith? [B,dS] bits degree
  let E ← erf? X
  return some ⟨C, (M.add E).sub (.point 1)⟩

theorem sample_sound {erf? erf}
    (sound : ∀ {I x J}, I.Contains x → erf? I = some J → RationalInterval.Contains J (erf x))
    {Ω z dS D bits degree i result} (hD : D.Contains dS)
    (h : sample? erf? Ω z D bits degree i = some result) :
    Option.Rel Sample.Covers result (realSample erf Ω z dS i) := by
  unfold sample? at h
  obtain ⟨C,hC,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨M,hM,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨B,hB,h⟩ := Option.bind_eq_some_iff.mp h
  have hc := Expr.contains_gridEnclosure? hC
  have hm := gridOne hc hM
  have hb := gridTwo hc hm hB
  split at h
  · cases h
    have hinvalid := hb.2.trans (by exact_mod_cast ‹B.hi ≤ 0› : (B.hi : ℝ) ≤ 0)
    simp [realSample, hinvalid]
  · obtain ⟨_, hg, h⟩ := Option.bind_eq_some_iff.mp h
    have hpos : 0 < B.lo := by simpa [guard] using hg
    have hlegal := (by exact_mod_cast hpos : (0:ℝ) < B.lo).trans_le hb.1
    obtain ⟨X,hX,h⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨E,hE,h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    simp only [realSample, ite_eq_right hlegal.not_ge, Id.run, pure]
    exact .some ⟨hc, RationalInterval.contains_sub (RationalInterval.contains_add hm
      (sound (gridTwo hb hD hX) hE)) (RationalInterval.contains_point 1)⟩

/-- A trial has a certified sign only when the interval excludes zero, or is exactly zero.
An uncertain sign can change the selected interpolation cell and is rejected. -/
private def nonnegative? (J : RationalInterval) : Option Bool :=
  if 0 ≤ J.lo then some true else if J.hi < 0 then some false else none

private theorem nonnegative_sound {J : RationalInterval} {x sign} (hx : J.Contains x)
    (h : nonnegative? J = some sign) : sign = decide (0 ≤ x) := by
  unfold nonnegative? at h
  split at h
  · cases h
    simp [(by exact_mod_cast ‹0 ≤ J.lo› : (0:ℝ) ≤ J.lo).trans hx.1]
  · split at h
    · cases h
      simp [(hx.2.trans_lt (by exact_mod_cast ‹J.hi < 0› : (J.hi : ℝ) < 0)).not_ge]
    · contradiction

abbrev Tagged (α : Type) := Bool × Sample α

private def row? (erf? : RationalInterval → Option RationalInterval) (Ω z : ℚ) (D : RationalInterval)
    (bits degree i : Nat) : Option (Option (Tagged RationalInterval)) := do
  let row ← sample? erf? Ω z D bits degree i
  match row with
  | none => return none
  | some r => return some (← nonnegative? r.residual, r)

noncomputable def realRow (erf : ℝ → ℝ) (Ω z : ℚ) (dS : ℝ) (i : Nat) : Option (Tagged ℝ) :=
  (realSample erf Ω z dS i).map (fun r => (decide (0 ≤ r.residual), r))

private def Related (r : Tagged RationalInterval) (s : Tagged ℝ) : Prop := r.1 = s.1 ∧ r.2.Covers s.2

private theorem row_sound {erf? erf}
    (sound : ∀ {I x J}, I.Contains x → erf? I = some J → RationalInterval.Contains J (erf x))
    {Ω z dS D bits degree i result} (hD : D.Contains dS)
    (h : row? erf? Ω z D bits degree i = some result) :
    Option.Rel Related result (realRow erf Ω z dS i) := by
  obtain ⟨r, hr, h⟩ := Option.bind_eq_some_iff.mp h
  have hs := sample_sound sound hD hr
  unfold realRow
  cases hReal : realSample erf Ω z dS i with
  | none =>
    rw [hReal] at hs
    cases hs
    cases h
    exact .none
  | some r' =>
    rw [hReal] at hs
    cases hs with
    | some hcov =>
      obtain ⟨s, hsign, h⟩ := Option.bind_eq_some_iff.mp h
      cases h
      exact .some ⟨nonnegative_sound hcov.2 hsign, hcov⟩

private def collect (row : Nat → Option (Option (Tagged RationalInterval))) :
    List Nat → Option (List (Tagged RationalInterval))
  | [] => some []
  | i :: is => do
    let r ← row i
    let tail ← collect row is
    return r.toList ++ tail

private theorem collect_sound {row rows}
    (hr : ∀ {i r}, row i = some r → Option.Rel Related r (rows i))
    {is out} (h : collect row is = some out) : List.Forall₂ Related out (is.filterMap rows) := by
  induction is generalizing out with
  | nil => cases h; exact .nil
  | cons i is ih =>
    obtain ⟨r, hr', h⟩ := Option.bind_eq_some_iff.mp h
    obtain ⟨tail, ht, h⟩ := Option.bind_eq_some_iff.mp h
    cases h
    have hs := hr hr'
    cases he : rows i with
    | none => rw [he] at hs; cases hs; simp [he]; exact ih ht
    | some s =>
      rw [he] at hs
      cases hs with
      | some hcov => simp [he]; exact ⟨hcov, ih ht⟩

inductive Pick (α : Type) where | point (r : Sample α) | between (left right : Sample α)

/-- Clamp outside the table and interpolate at the unique negative-to-nonnegative crossing.
An additional sign change violates the lookup's ordering requirement. -/
private def select {α} (previous : Option (Sample α)) : List (Tagged α) → Option (Pick α)
  | [] => previous.map Pick.point
  | (false,r) :: tail => select (some r) tail
  | (true,r) :: tail =>
    if tail.all Prod.fst then some (match previous with | none => .point r | some p => .between p r)
    else none

private inductive PickRel : Pick RationalInterval → Pick ℝ → Prop
  | point {a b} (h : a.Covers b) : PickRel (.point a) (.point b)
  | between {a b c d} (h : a.Covers c) (g : b.Covers d) : PickRel (.between a b) (.between c d)

private theorem signs_agree {xs ys} (h : List.Forall₂ Related xs ys) :
    xs.all Prod.fst = ys.all Prod.fst := by
  induction h with
  | nil => rfl
  | cons hr _ ih => simp [hr.1,ih]

private theorem select_sound {xs ys} (h : List.Forall₂ Related xs ys)
    {a b} (hp : Option.Rel Sample.Covers a b) :
    Option.Rel PickRel (select a xs) (select b ys) := by
  induction h generalizing a b with
  | nil =>
    cases hp with
    | none => exact .none
    | some hx => exact .some (.point hx)
  | @cons x y xs ys hr ht ih =>
    rcases x with ⟨flag,r⟩
    rcases y with ⟨flag',s⟩
    obtain ⟨rfl, hs⟩ := hr
    cases flag with
    | false => exact ih (.some hs)
    | true =>
      simp only [select, signs_agree ht]
      split
      · cases hp with
        | none => exact .some (.point hs)
        | some hx => exact .some (.between hx hs)
      · exact .none

/-- y coordinate at residual zero for a linear interpolation between the bracketing samples. -/
private def lineExpr : Expr :=
  let l : Expr := .var 0
  let r : Expr := .var 1
  let f : Expr := .var 2
  let g : Expr := .var 3
  l + (r-l) * (.unary .neg f / (g-f))

private def interpolate? (bits degree : Nat) : Pick RationalInterval → Option RationalInterval
  | .point a => some a.concentration
  | .between a b =>
    lineExpr.gridEnclosureWith? [a.concentration,b.concentration,a.residual,b.residual] bits degree

noncomputable def interpReal : Pick ℝ → ℝ
  | .point a => a.concentration
  | .between a b =>
    a.concentration + (b.concentration-a.concentration)*(-a.residual/(b.residual-a.residual))

private theorem interp_sound {p q : _} (h : PickRel p q) {bits degree I}
    (he : interpolate? bits degree p = some I) : I.Contains (interpReal q) := by
  cases h with
  | point ha => cases he; exact ha.1
  | @between a b c d hA hB =>
    let env := fun i : Nat => match i with
      | 0 => c.concentration | 1 => d.concentration | 2 => c.residual | _ => d.residual
    apply Expr.contains_gridEnclosureWith? (values := env) (h := he)
    intro i J hj
    rcases i with _|_|_|_|n
    · simp at hj; subst J; exact hA.1
    · simp at hj; subst J; exact hB.1
    · simp at hj; subst J; exact hA.2
    · simp at hj; subst J; exact hB.2
    · simp at hj

noncomputable def realGrid (erf : ℝ → ℝ) (Ω z : ℚ) (dS : ℝ) : Option ℝ :=
  (select none ((List.range 100).filterMap (realRow erf Ω z dS))).map interpReal

def grid? (erf? : RationalInterval → Option RationalInterval) (Ω z : ℚ) (dS : RationalInterval)
    (bits degree : Nat) (tol : ℚ) : Option RationalInterval := do
  guard (0 < Ω ∧ Ω ≤ 1 ∧ 0 ≤ z ∧ 0 < dS.lo ∧ dS.lo ≤ dS.hi)
  let rows ← collect (row? erf? Ω z dS bits degree) (List.range 100)
  let p ← select none rows
  let I ← interpolate? bits degree p
  SashimiNumerics.within? tol I

theorem grid_sound {erf? erf}
    (sound : ∀ {I x J}, I.Contains x → erf? I = some J → RationalInterval.Contains J (erf x))
    {Ω z dS D bits degree tol I} (hD : D.Contains dS) (h : grid? erf? Ω z D bits degree tol = some I) :
    ∃ c, realGrid erf Ω z dS = some c ∧ I.Contains c ∧ I.hi - I.lo ≤ tol := by
  obtain ⟨_,_,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨rows,hr,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨p,hp,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨J,hj,hw⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨rfl, hw⟩ := SashimiNumerics.within_sound hw
  have hrows := collect_sound (row_sound sound hD) hr
  have hs := select_sound hrows (Option.Rel.none (r := Sample.Covers))
  rw [hp] at hs
  cases he : select none _ with
  | none => rw [he] at hs; cases hs
  | some q =>
    rw [he] at hs
    cases hs with
    | some hcov => exact ⟨_, by simp [realGrid,he], interp_sound hcov hj, hw⟩


/-- Concentration conditional on the top-hat variance difference at the selected mass bin.
All 100 trial concentrations, the validity mask, signs and interpolation are checked. -/
def concentration? (Ω z : ℚ) (dS : RationalInterval) (bits degree : Nat) (tol : ℚ) :=
  grid? (fun I => SashimiErf.erf? I bits degree) Ω z dS bits degree tol

theorem concentration_sound {Ω z dS D bits degree tol I} (hD : D.Contains dS)
    (h : concentration? Ω z D bits degree tol = some I) :
    ∃ c, realGrid SashimiErf.erf Ω z dS = some c ∧ I.Contains c ∧ I.hi - I.lo ≤ tol :=
  grid_sound (fun hx hI => SashimiErf.erf_sound hx hI) hD h

#print axioms concentration_sound

end SashimiHaloWarm
