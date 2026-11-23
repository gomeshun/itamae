import SashimiHalo

/-! # Inverting the two NFW mass relations
A candidate bracket may be supplied by any numerical root finder. Two outward evaluations at
its endpoints suffice. The equations have a unique nonnegative or positive root, respectively. -/

namespace SashimiHalo

open FloatLib.Numerics FloatLib.Numerics.Interval Set

noncomputable def nfw (x : ℝ) := Real.log (1 + x) - x / (1 + x)

theorem eval_fc (x : Expr) (env : Nat → ℝ) : (fc x).eval env = nfw (x.eval env) := by rfl

private theorem nfw_deriv {x : ℝ} (hx : 0 ≤ x) :
    HasDerivAt nfw (x / (1 + x)^2) x := by
  have h1 := (hasDerivAt_const x (1 : ℝ)).add (hasDerivAt_id x)
  have hn : 1 + x ≠ 0 := by positivity
  convert (h1.log hn).sub ((hasDerivAt_id x).div h1 hn) using 1 <;> try rfl
  simp only [Pi.add_apply, id_eq]
  field_simp
  ring

theorem nfw_strictMono : StrictMonoOn nfw (Ici 0) := by
  apply strictMonoOn_of_deriv_pos (convex_Ici 0)
  · exact fun x hx => (nfw_deriv hx).continuousAt.continuousWithinAt
  · intro x hx
    have hx' : 0 < x := by simpa using hx
    rw [(nfw_deriv hx'.le).deriv]
    positivity

/-- Dividing the NFW mass function by c³ yields a strictly decreasing positive-radius function. -/
theorem meanProfile_strictAnti : StrictAntiOn (fun x : ℝ => nfw x / x^3) (Ioi 0) := by
  let gap := fun x : ℝ => nfw x - x^2 / (3 * (1+x)^2)
  have dgap {x : ℝ} (hx : 0 ≤ x) :
      HasDerivAt gap (x * (1+3*x) / (3*(1+x)^3)) x := by
    have h1 := (hasDerivAt_const x (1 : ℝ)).add (hasDerivAt_id x)
    have hd : 3 * (1+x)^2 ≠ 0 := by positivity
    convert (nfw_deriv hx).sub (((hasDerivAt_id x).pow 2).div
      (((h1.pow 2).const_mul 3)) hd) using 1 <;> try rfl
    simp only [Pi.add_apply, Pi.pow_apply, id_eq]
    field_simp
    ring
  have hgap : StrictMonoOn gap (Ici 0) := by
    apply strictMonoOn_of_deriv_pos (convex_Ici 0)
    · exact fun x hx => (dgap hx).continuousAt.continuousWithinAt
    · intro x hx
      have hx' : 0 < x := by simpa using hx
      rw [(dgap hx'.le).deriv]
      positivity
  have df {x : ℝ} (hx : 0 < x) :=
    (nfw_deriv hx.le).div ((hasDerivAt_id x).pow 3) (show x ^ 3 ≠ 0 by positivity)
  apply strictAntiOn_of_deriv_neg (convex_Ioi 0)
  · exact fun x hx => (df hx).continuousAt.continuousWithinAt
  · intro x hx
    have hx' : 0 < x := by simpa using hx
    have hg := hgap (by simp : (0 : ℝ) ∈ Ici 0) hx'.le hx'
    have h0 : gap 0 = 0 := by simp [gap,nfw]
    rw [h0] at hg
    change 0 < nfw x - x^2 / (3*(1+x)^2) at hg
    have hnum := (div_lt_iff₀ (show 0 < 3*(1+x)^2 by positivity)).mp (sub_pos.mp hg)
    change deriv (nfw / id ^ 3) x < 0
    rw [(df hx').deriv]
    simp only [Pi.pow_apply, id_eq]
    apply div_neg_of_neg_of_pos _ (by positivity)
    simp only [Nat.cast_ofNat, Nat.reduceSub, mul_one]
    rw [div_mul_eq_mul_div, sub_neg, div_lt_iff₀ (by positivity)]
    nlinarith [mul_lt_mul_of_pos_right hnum (sq_pos_of_pos hx')]

inductive RootKind where
  | enclosed | mean
  deriving DecidableEq, Repr

def rootExpr : RootKind → Expr → Expr
  | .enclosed, c => fc c
  | .mean, c => .binary .div (fc c) (.unary (.pow 3) c)

noncomputable def rootModel : RootKind → ℝ → ℝ
  | .enclosed, x => nfw x
  | .mean, x => nfw x / x^3

theorem rootExpr_eval (kind : RootKind) (q : ℚ) :
    (rootExpr kind (.const q)).eval (fun _ => 0) = rootModel kind q := by
  cases kind <;> rfl

def domain (kind : RootKind) (x : ℝ) : Prop :=
  match kind with | .enclosed => 0 ≤ x | .mean => 0 < x

private theorem enclosedInversion {l u y : ℝ} (hl : 0 ≤ l) (hlu : l ≤ u)
    (hy : nfw l ≤ y ∧ y ≤ nfw u) :
    ∃ x ∈ Icc l u, 0 ≤ x ∧ nfw x = y ∧ ∀ t, 0 ≤ t → nfw t = y → t = x := by
  have hc : ContinuousOn nfw (Icc l u) :=
    fun t ht => (nfw_deriv (hl.trans ht.1)).continuousAt.continuousWithinAt
  obtain ⟨x, hx, he⟩ := intermediate_value_Icc hlu hc hy
  exact ⟨x, hx, hl.trans hx.1, he, fun t ht hy' => nfw_strictMono.injOn ht (hl.trans hx.1) (hy'.trans he.symm)⟩

private theorem meanInversion {l u y : ℝ} (hl : 0 < l) (hlu : l ≤ u)
    (hy : nfw u / u^3 ≤ y ∧ y ≤ nfw l / l^3) :
    ∃ x ∈ Icc l u, 0 < x ∧ nfw x / x^3 = y ∧ ∀ t, 0 < t → nfw t / t^3 = y → t = x := by
  have hc : ContinuousOn (fun t => nfw t / t^3) (Icc l u) := by
    intro t ht
    have ht' := hl.trans_le ht.1
    exact (nfw_deriv ht'.le).continuousAt.continuousWithinAt.div (continuousAt_id.pow 3).continuousWithinAt
      (pow_ne_zero _ ht'.ne')
  obtain ⟨x, hx, he⟩ := intermediate_value_Icc' hlu hc hy
  exact ⟨x, hx, hl.trans_le hx.1, he, fun t ht hy' =>
    meanProfile_strictAnti.injOn ht (hl.trans_le hx.1) (hy'.trans he.symm)⟩

/-- Endpoints, exact target expression and tolerance of a positive NFW concentration bracket.
For virial conversion the target is Δρcrit/(3ρs). -/
def rootBracket? (kind : RootKind) (I : RationalInterval) (target : Expr) (tol : ℚ)
    (bits degree : Nat) : Bool :=
  if 0 ≤ I.lo ∧ I.lo ≤ I.hi ∧ I.hi - I.lo ≤ tol ∧ (kind = .mean → 0 < I.lo) then
    match (rootExpr kind (.const I.lo)).gridEnclosure? bits degree,
          (rootExpr kind (.const I.hi)).gridEnclosure? bits degree, target.gridEnclosure? bits degree with
    | some A, some B, some T =>
      match kind with
      | .enclosed => decide (A.hi ≤ T.lo ∧ T.hi ≤ B.lo)
      | .mean => decide (B.hi ≤ T.lo ∧ T.hi ≤ A.lo)
    | _, _, _ => false
  else false

theorem rootBracket_sound {kind I target tol bits degree} (h : rootBracket? kind I target tol bits degree) :
    I.hi - I.lo ≤ tol ∧ ∃ x, I.Contains x ∧ domain kind x ∧ rootModel kind x = target.eval (fun _ => 0) ∧
      ∀ t, domain kind t → rootModel kind t = target.eval (fun _ => 0) → t = x := by
  unfold rootBracket? at h
  split at h
  · rename_i valid
    obtain ⟨hl, hlu, hw, hm⟩ := valid
    refine ⟨hw, ?_⟩
    split at h <;> try contradiction
    rename_i A B T ha hb ht
    have hA := Expr.contains_gridEnclosure? ha
    have hB := Expr.contains_gridEnclosure? hb
    have hT := Expr.contains_gridEnclosure? ht
    rw [rootExpr_eval] at hA hB
    have hlu' : (I.lo : ℝ) ≤ I.hi := by exact_mod_cast hlu
    cases kind with
    | enclosed =>
      obtain ⟨hlow, hhigh⟩ := of_decide_eq_true h
      exact enclosedInversion (by exact_mod_cast hl) hlu'
        ⟨hA.2.trans ((by exact_mod_cast hlow : (A.hi : ℝ) ≤ T.lo).trans hT.1),
         hT.2.trans ((by exact_mod_cast hhigh : (T.hi : ℝ) ≤ B.lo).trans hB.1)⟩
    | mean =>
      obtain ⟨hlow, hhigh⟩ := of_decide_eq_true h
      exact meanInversion (by exact_mod_cast hm rfl) hlu'
        ⟨hB.2.trans ((by exact_mod_cast hlow : (B.hi : ℝ) ≤ T.lo).trans hT.1),
         hT.2.trans ((by exact_mod_cast hhigh : (T.hi : ℝ) ≤ A.lo).trans hA.1)⟩
  · contradiction

/-- Propose a narrower bracket, retaining the original when midpoint and target overlap.
The public solver checks both endpoints of the proposal independently. -/
private def bisect (kind : RootKind) (T : RationalInterval) (bits degree : Nat) (tol : ℚ) :
    Nat → RationalInterval → RationalInterval
  | 0, I => I
  | n+1, I => Id.run do
    if I.hi - I.lo ≤ tol then return I
    let m := (I.lo + I.hi) / 2
    let some J := (rootExpr kind (.const m)).gridEnclosure? bits degree | return I
    if J.hi ≤ T.lo then
      return bisect kind T bits degree tol n (if kind = .enclosed then ⟨m, I.hi⟩ else ⟨I.lo, m⟩)
    if T.hi ≤ J.lo then
      return bisect kind T bits degree tol n (if kind = .enclosed then ⟨I.lo, m⟩ else ⟨m, I.hi⟩)
    return I

def invert? (kind : RootKind) (target : Expr) (initial : RationalInterval) (bits degree fuel : Nat)
    (tol : ℚ) : Option RationalInterval := do
  let T ← target.gridEnclosure? bits degree
  let I := bisect kind T bits degree tol fuel initial
  if rootBracket? kind I target tol bits degree then some I else none

theorem invert_sound {kind target initial bits degree fuel tol I}
    (h : invert? kind target initial bits degree fuel tol = some I) :
    I.hi - I.lo ≤ tol ∧ ∃ x, I.Contains x ∧ domain kind x ∧ rootModel kind x = target.eval (fun _ => 0) ∧
      ∀ t, domain kind t → rootModel kind t = target.eval (fun _ => 0) → t = x := by
  unfold invert? at h
  obtain ⟨T, _, hI⟩ := Option.bind_eq_some_iff.mp h
  dsimp at hI
  split at hI
  · cases hI
    exact rootBracket_sound (by assumption)
  · contradiction


/-- The NFW virial equation after cancelling rhoCrit, marked by its overdensity Δ. -/
def virialTarget (c200 Δ : ℚ) : Expr :=
  .binary .div (.const (Δ/200) * fc (.const c200)) (.const (c200^3))

/-- Mvir = M200 f(cvir)/f(c200), with cvir supplied by its certified interval. -/
def massAtRoot (m c200 : ℚ) : Expr :=
  .binary .div (.const m * fc (.var 0)) (fc (.const c200))

structure VirialResult where
  cvir : RationalInterval
  mass : RationalInterval
deriving Repr

def virial? (m c200 Δ : ℚ) (initial : RationalInterval) (bits degree fuel : Nat)
    (ctol mtol : ℚ) : Option VirialResult := do
  guard (0 < m ∧ 0 < c200 ∧ 0 < Δ)
  let c ← invert? .mean (virialTarget c200 Δ) initial bits degree fuel ctol
  let J ← (massAtRoot m c200).gridEnclosureWith? [c] bits degree
  let mass ← SashimiNumerics.within? mtol J
  return ⟨c, mass⟩

theorem virial_sound {m c200 Δ initial bits degree fuel ctol mtol R}
    (h : virial? m c200 Δ initial bits degree fuel ctol mtol = some R) :
    R.cvir.hi - R.cvir.lo ≤ ctol ∧ R.mass.hi - R.mass.lo ≤ mtol ∧
    ∃ c : ℝ, R.cvir.Contains c ∧ 0 < c ∧ nfw c / c^3 = (Δ : ℝ)/200 * nfw c200 / (c200 : ℝ)^3 ∧
      (∀ t : ℝ, 0 < t → nfw t / t^3 = (Δ : ℝ)/200 * nfw c200 / (c200 : ℝ)^3 → t = c) ∧
      R.mass.Contains ((m : ℝ) * nfw c / nfw c200) := by
  unfold virial? at h
  obtain ⟨_, _, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨C, hc, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨J, hj, h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨W, hw, h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  obtain ⟨rfl, hwidth⟩ := SashimiNumerics.within_sound hw
  obtain ⟨hwc,c,hcbox,hpos,heq,huniq⟩ := invert_sound hc
  have ht : (virialTarget c200 Δ).eval (fun _ => 0) = (Δ : ℝ)/200 * nfw c200 / (c200 : ℝ)^3 := by
    change (↑(Δ/200 : ℚ) : ℝ) * nfw c200 / (↑(c200^3 : ℚ) : ℝ) = _
    push_cast
    rfl
  rw [ht] at heq huniq
  refine ⟨hwc,hwidth,c,hcbox,hpos,heq,huniq,?_⟩
  apply Expr.contains_gridEnclosureWith? (values := fun _ => c) (h := hj)
  intro i I hI
  cases i with
  | zero => simp at hI; subst I; exact hcbox
  | succ i => simp at hI

#print axioms invert_sound
#print axioms virial_sound

end SashimiHalo
