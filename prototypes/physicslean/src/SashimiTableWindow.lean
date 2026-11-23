import SashimiTable
import SashimiSpectrumWindow
import SashimiSpectrumFdmWindow
import FloatLib.Numerics.Enclosure.Rational.Converges

/-! Whole-cell enclosures for the variance of a tabulated (linear or log-log)
spectrum. The supplied rational interval for R encloses every kR used in a
cell. Low-k and high-k analytic bounds avoid subtractive cancellation and
resolve oscillatory tails without sampling them. -/
namespace SashimiTableWindow
open FloatLib.Numerics FloatLib.Numerics.Interval SashimiTable SashimiSpectrum MeasureTheory Set Real

private def clamp (A B : RationalInterval) : RationalInterval := ⟨max A.lo B.lo, min A.hi B.hi⟩
private theorem clamp_sound {A B : RationalInterval} {x : ℝ} (a : A.Contains x) (b : B.Contains x) :
    (clamp A B).Contains x := by
  simpa [clamp, RationalInterval.Contains] using And.intro (max_le a.1 b.1) (le_min a.2 b.2)

private def rounded (bits : Nat) (I : RationalInterval) : RationalInterval :=
  let s : ℚ := 2^bits
  ⟨⌊I.lo*s⌋/s, ⌈I.hi*s⌉/s⟩

private theorem round_sound (bits : Nat) {I : RationalInterval} {x : ℝ} (h : I.Contains x) :
    (rounded bits I).Contains x := by
  let s : ℚ := 2^bits
  have hs : 0 < s := by dsimp [s]; positivity
  have hl : (⌊I.lo*s⌋ : ℚ)/s ≤ I.lo := (div_le_iff₀ hs).2 (Int.floor_le _)
  have hu : I.hi ≤ (⌈I.hi*s⌉ : ℚ)/s := (le_div_iff₀ hs).2 (Int.le_ceil _)
  exact ⟨(by exact_mod_cast hl : (((rounded bits I).lo : ℚ) : ℝ) ≤ I.lo).trans h.1,
    h.2.trans (by exact_mod_cast hu : (I.hi : ℝ) ≤ ((rounded bits I).hi : ℝ))⟩

private def analytic (l u : ℚ) : RationalInterval :=
  let decay := if 0 < l then min 1 ((3*(1+u)/l^3)^2) else 1
  ⟨(max 0 (1-u^2/10))^2, decay⟩

private theorem analytic_sound {x : ℝ} {l u : ℚ} (hl : 0 ≤ l) (hx : (l:ℝ) ≤ x)
    (hu : x ≤ (u:ℝ)) : (analytic l u).Contains ((topHat x)^2) := by
  have hp : 0 ≤ x := (by exact_mod_cast hl : (0:ℝ) ≤ l).trans hx
  obtain ⟨ha,hb⟩ := topHat_bounds x hp
  have h2 : (topHat x)^2 ≤ 1 := by nlinarith [mul_nonneg (by linarith : 0 ≤ 1-topHat x) (by linarith : 0 ≤ 1+topHat x)]
  constructor
  · dsimp [analytic]
    push_cast
    rcases le_total (1-(u:ℝ)^2/10) 0 with h | h
    · rw [max_eq_left h]; simpa using sq_nonneg (topHat x)
    · rw [max_eq_right h]
      apply pow_le_pow_left₀ h
      have h' := pow_le_pow_left₀ hp hu 2
      linarith [topHat_taylor_lower x hp]
  · unfold analytic; dsimp
    split
    · rename_i hl'
      have hlr : (0:ℝ) < l := by exact_mod_cast hl'
      have hxr := hlr.trans_le hx
      have hd := topHat_decay x hxr
      have hU := hp.trans hu
      have hub : |topHat x| ≤ 3*(1+(u:ℝ))/(l:ℝ)^3 := hd.trans (by gcongr)
      push_cast
      apply le_min h2
      simpa only [sq_abs] using pow_le_pow_left₀ (abs_nonneg _) hub 2
    · simpa using h2

/-- Intersect the analytic bounds with direct trigonometry in the middle of
the spectrum. Outside the middle band the algebraic bounds avoid expensive
period reduction and preserve the small-k limit. -/
def windowSq (l u : ℚ) (bits degree : Nat) : RationalInterval :=
  if 1/2 ≤ l ∧ u ≤ 32 then
    match SashimiNumerics.exprRange ((topHatExpr (.var 0)).unary (.pow 2)) bits degree l u with
    | some I => clamp (analytic l u) I
    | none => analytic l u
  else analytic l u

theorem windowSq_sound {x : ℝ} {l u : ℚ} (hl : 0 ≤ l) (hx : (l:ℝ) ≤ x)
    (hu : x ≤ (u:ℝ)) (bits degree : Nat) : (windowSq l u bits degree).Contains ((topHat x)^2) := by
  have h := analytic_sound hl hx hu
  unfold windowSq
  split
  · rename_i hmid
    split
    · rename_i I hi
      apply clamp_sound h
      have hpos : 0 < x := lt_of_lt_of_le (by exact_mod_cast (lt_of_lt_of_le (by norm_num : (0:ℚ)<1/2) hmid.1)) hx
      have hh := SashimiNumerics.exprRange_sound _ _ _ _ _ _ hi x ⟨hx,hu⟩
      change I.Contains (((topHatExpr (.var 0)).eval (fun _ => x))^2) at hh
      rwa [eval_topHatExpr (show ((Expr.var 0).eval fun _ => x) ≠ 0 from hpos.ne')] at hh
    · exact h
  · exact h

/-- Convex-hull enclosure. If the left logarithmic ordinate could vanish,
the extra upper endpoint covers the total real convention log(0)=0. -/
private def hull (mode : Interpolation) (L U : RationalInterval) : RationalInterval :=
  ⟨min L.lo U.lo, if mode = .logLog ∧ L.lo ≤ 0 then max 1 (max L.hi U.hi) else max L.hi U.hi⟩

private theorem between {l r z w : ℝ} (hl : z ≤ l) (hr : z ≤ r)
    (hw : 0 ≤ w) (hw' : w ≤ 1) : z ≤ l + w*(r-l) := by
  nlinarith [mul_nonneg (sub_nonneg.mpr hl) (sub_nonneg.mpr hw'), mul_nonneg (sub_nonneg.mpr hr) hw]

private theorem interp_between {l r p q x : ℝ} (hl : 0 < l) (hlr : l < r) (hx : x ∈ Icc l r)
    (mode : Interpolation) (hp : 0 < p) (hq : 0 < q) :
    min p q ≤ interpReal mode l r p q x ∧ interpReal mode l r p q x ≤ max p q := by
  cases mode with
  | linear =>
    have h0 : 0 ≤ (x-l)/(r-l) := div_nonneg (sub_nonneg.mpr hx.1) (sub_nonneg.mpr hlr.le)
    have h1 : (x-l)/(r-l) ≤ 1 := (div_le_one (sub_pos.mpr hlr)).2 (by linarith [hx.2])
    exact ⟨between (min_le_left _ _) (min_le_right _ _) h0 h1,
      by have h := between (neg_le_neg (le_max_left p q)) (neg_le_neg (le_max_right p q)) h0 h1
         unfold interpReal; linarith⟩
  | logLog =>
    have hx0 := hl.trans_le hx.1
    have hb : 0 < log (r/l) := log_pos ((one_lt_div hl).2 hlr)
    have w0 : 0 ≤ log (x/l) / log (r/l) :=
      div_nonneg (log_nonneg ((one_le_div hl).2 hx.1)) hb.le
    have w1 : log (x/l) / log (r/l) ≤ 1 := (div_le_one hb).2
      (log_le_log (div_pos hx0 hl) (div_le_div_of_nonneg_right hx.2 hl.le))
    simp only [interpReal, log_div hq.ne' hp.ne']
    rw [← Real.log_le_log_iff (lt_min hp hq) (exp_pos _), ← Real.log_le_log_iff (exp_pos _) (hp.trans_le (le_max_left p q)), log_exp]
    constructor
    · exact between (log_le_log (lt_min hp hq) (min_le_left _ _)) (log_le_log (lt_min hp hq) (min_le_right _ _)) w0 w1
    · have h := between (neg_le_neg (log_le_log hp (le_max_left p q))) (neg_le_neg (log_le_log hq (le_max_right p q))) w0 w1
      linarith

private theorem hull_sound {l r p q x : ℝ} {L U : RationalInterval} (mode : Interpolation)
    (hl : 0 < l) (hr : l < r) (hx : x ∈ Icc l r) (p0 : 0 ≤ p) (q0 : 0 ≤ q)
    (hL : L.Contains p) (hU : U.Contains q) :
    (hull mode L U).Contains (interpReal mode l r p q x) := by
  have hmid : min (L.lo:ℝ) U.lo ≤ min p q := min_le_min hL.1 hU.1
  have ht : max p q ≤ max (L.hi:ℝ) U.hi := max_le_max hL.2 hU.2
  have hhi : max (L.hi:ℝ) U.hi ≤ ((hull mode L U).hi:ℝ) := by
    unfold hull; dsimp; split <;> push_cast <;> simp
  have between' (hh : min p q ≤ interpReal mode l r p q x ∧ interpReal mode l r p q x ≤ max p q) :
      (hull mode L U).Contains (interpReal mode l r p q x) :=
    ⟨by simpa only [hull, Rat.cast_min] using hmid.trans hh.1, hh.2.trans (ht.trans hhi)⟩
  cases mode with
  | linear =>
    have w0 := div_nonneg (sub_nonneg.mpr hx.1) (sub_nonneg.mpr hr.le)
    have w1 : (x-l)/(r-l) ≤ 1 := (div_le_one (sub_pos.mpr hr)).2 (by linarith [hx.2])
    apply between'
    exact ⟨between (min_le_left _ _) (min_le_right _ _) w0 w1,
      by have h := between (neg_le_neg (le_max_left p q)) (neg_le_neg (le_max_right p q)) w0 w1
         unfold interpReal; linarith⟩
  | logLog =>
    rcases p0.eq_or_lt with hp | hp
    · subst p
      have hlr : L.lo ≤ 0 := by exact_mod_cast hL.1
      have hl' := (min_le_left (L.lo:ℝ) (U.lo:ℝ)).trans hL.1
      simpa [interpReal, RationalInterval.Contains, hull, hlr] using
        (show min (L.lo:ℝ) U.lo ≤ 1 ∧ (1:ℝ) ≤ max 1 (max (L.hi:ℝ) U.hi) from ⟨hl'.trans zero_le_one, le_max_left _ _⟩)
    · rcases q0.eq_or_lt with hq | hq
      · subst q
        apply between'
        simp [interpReal, exp_log hp]
      · exact between' (interp_between hl hr hx .logLog hp hq)

noncomputable def cellValue (j : Nat) (mode : Interpolation) (Y : Knot → ℝ)
    (R : ℝ) (c : Cell) (k : ℝ) : ℝ :=
  k^(2+2*j) * interpReal mode c.left.k c.right.k (Y c.left) (Y c.right) k * topHat (k*R)^2

private theorem measurable_cell (j : Nat) (mode : Interpolation) (Y : Knot → ℝ) (R : ℝ) (c : Cell) : Measurable (cellValue j mode Y R c) := by
  have hw : Measurable topHat := by
    unfold topHat
    exact Measurable.ite (measurableSet_singleton 0) measurable_const (by fun_prop)
  have hp : Measurable (interpReal mode c.left.k c.right.k (Y c.left) (Y c.right)) := by
    cases mode <;> unfold interpReal <;> fun_prop
  exact ((measurable_id.pow_const _).mul hp).mul ((hw.comp (measurable_id.mul_const R)).pow_const _)

private def panel? (j : Nat) (mode : Interpolation) (R : RationalInterval) (c : Cell)
    (L U : RationalInterval) (bits degree : Nat) (l u : ℚ) : Option RationalInterval :=
  if 0 < c.left.k ∧ c.left.k < c.right.k ∧ c.left.k ≤ l ∧ u ≤ c.right.k ∧ 0 ≤ R.lo then
    let K : RationalInterval := ⟨l^(2+2*j), u^(2+2*j)⟩
    let W := windowSq (l*R.lo) (u*R.hi) bits degree
    some ((K.mul (hull mode L U)).mul W)
  else none

private theorem panel_sound {Y : Knot → ℝ} {R L U : RationalInterval} {v : ℝ} {c : Cell}
    (hR : R.Contains v) (hL : L.Contains (Y c.left)) (hU : U.Contains (Y c.right))
    (pl : 0 ≤ Y c.left) (pu : 0 ≤ Y c.right)
    (j mode bits degree l u I)
    (h : panel? j mode R c L U bits degree l u = some I) (k : ℝ) (hk : k ∈ Icc (l:ℝ) u) :
    I.Contains (cellValue j mode Y v c k) := by
  unfold panel? at h
  split at h
  · rename_i hv
    obtain ⟨hl,hr,hll,hur,hRl⟩ := hv
    cases h
    have hp : 0 ≤ (l:ℝ) := by exact_mod_cast (hl.le.trans hll)
    have hk0 := hp.trans hk.1
    have hR0 : 0 ≤ (R.lo:ℝ) := by exact_mod_cast hRl
    have hc : k ∈ Icc (c.left.k:ℝ) c.right.k :=
      ⟨(by exact_mod_cast hll : (c.left.k:ℝ) ≤ l).trans hk.1, hk.2.trans (by exact_mod_cast hur)⟩
    apply RationalInterval.contains_mul
    · apply RationalInterval.contains_mul
      · constructor <;> push_cast <;> gcongr <;> first | exact hk.1 | exact hk.2
      · exact hull_sound mode (by exact_mod_cast hl) (by exact_mod_cast hr) hc pl pu hL hU
    · apply windowSq_sound (mul_nonneg (hl.le.trans hll) hRl)
      · push_cast
        exact mul_le_mul hk.1 hR.1 hR0 hk0
      · push_cast
        exact mul_le_mul hk.2 hR.2 (hR0.trans hR.1) (hk0.trans hk.2)
  · contradiction

/-- Per-knot evaluation is supplied explicitly, so positive spectra, Gaussian,
WDM and FDM power transfers can share the same quadrature proof. -/
def sumCells? (j : Nat) (mode : Interpolation) (power : Knot → Option RationalInterval)
    (R : RationalInterval) (a b : ℚ) (depth bits degree : Nat) :
    List Cell → Option RationalInterval
  | [] => some (.point 0)
  | c :: cs => do
    let l := max a c.left.k
    let u := min b c.right.k
    if l < u then
      let L ← power c.left
      let U ← power c.right
      let I ← SashimiNumerics.rectangles? (panel? j mode R c L U bits degree) l u depth
      let T ← sumCells? j mode power R a b depth bits degree cs
      return (rounded bits I).add T
    else sumCells? j mode power R a b depth bits degree cs

/-- Combine independently checked groups of adjacent cells without re-evaluating
their quadratures or transfer functions in the kernel. -/
theorem sumCells?_append_of {j mode power R a b depth bits degree X Z cs ds}
    (hX : sumCells? j mode power R a b depth bits degree cs = some X)
    (hZ : sumCells? j mode power R a b depth bits degree ds = some Z) :
    sumCells? j mode power R a b depth bits degree (cs ++ ds) = some (X.add Z) := by
  induction cs generalizing X with
  | nil =>
    simp only [sumCells?, Option.some.injEq] at hX
    subst X
    simpa [RationalInterval.add_def, RationalInterval.point] using hZ
  | cons c cs ih =>
    simp only [List.cons_append, sumCells?] at hX ⊢
    split at hX
    · rename_i hz
      rw [ite_eq_left hz]
      obtain ⟨L,hL,rest⟩ := Option.bind_eq_some_iff.mp hX
      obtain ⟨U,hU,rest⟩ := Option.bind_eq_some_iff.mp rest
      obtain ⟨B,hB,rest⟩ := Option.bind_eq_some_iff.mp rest
      obtain ⟨T,hT,h⟩ := Option.bind_eq_some_iff.mp rest
      cases h
      simp [hL, hU, ih hT, hB, RationalInterval.add_def, add_assoc]
    · rename_i hz
      rw [ite_eq_right hz]
      exact ih hX

/-- Sum of genuine clipped interval integrals in k, prior to division by 2π². -/
noncomputable def sumCells (j : Nat) (mode : Interpolation) (Y : Knot → ℝ) (R : ℝ) (a b : ℚ) :
    List Cell → ℝ
  | [] => 0
  | c :: cs =>
    let l := max a c.left.k
    let u := min b c.right.k
    if l < u then (∫ k in (l:ℝ)..(u:ℝ), cellValue j mode Y R c k) + sumCells j mode Y R a b cs
    else sumCells j mode Y R a b cs

theorem sumCells_sound {j mode power R a b depth bits degree cs I} {Y : Knot → ℝ} {v : ℝ}
    (hR : R.Contains v) (hp : ∀ k B, power k = some B → 0 ≤ Y k ∧ B.Contains (Y k))
    (h : sumCells? j mode power R a b depth bits degree cs = some I) :
    I.Contains (sumCells j mode Y v a b cs) := by
  induction cs generalizing I with
  | nil =>
    simp [sumCells?] at h
    cases h
    simp [sumCells, RationalInterval.Contains, RationalInterval.point]
  | cons c cs ih =>
    simp only [sumCells?, sumCells] at h ⊢
    split at h
    · rename_i hpos
      rw [ite_eq_left hpos]
      obtain ⟨L, hL, hrest⟩ := Option.bind_eq_some_iff.mp h
      obtain ⟨U, hU, hrest⟩ := Option.bind_eq_some_iff.mp hrest
      obtain ⟨J, hJ, hrest⟩ := Option.bind_eq_some_iff.mp hrest
      obtain ⟨T, hT, hI⟩ := Option.bind_eq_some_iff.mp hrest
      cases hI
      obtain ⟨pl, hl⟩ := hp _ _ hL
      obtain ⟨pu, hu⟩ := hp _ _ hU
      have rb := panel_sound hR hl hu pl pu j mode bits degree
      have hi := SashimiNumerics.rectangles_integrable (measurable_cell _ _ _ _ _) rb hpos.le hJ
      exact RationalInterval.contains_add (round_sound bits (SashimiNumerics.rectangles_sound rb hpos.le hi hJ)) (ih hT)
    · rename_i hfalse
      rw [ite_eq_right hfalse]
      exact ih h

private def normalizer : Expr := .binary .div (.var 0)
  (.binary .mul (.const 2) (.unary (.pow 2) (.binary .mul (.const 4) (.unary .atan (.const 1)))))

/-- Normalize a certified raw sum and check the requested width. -/
def finish? (S : RationalInterval) (tol : ℚ) (bits degree : Nat) : Option RationalInterval := do
  let I ← normalizer.gridEnclosureWith? [S] bits degree
  SashimiNumerics.within? tol I

def integrate? (j : Nat) (mode : Interpolation) (power : Knot → Option RationalInterval)
    (rows : List Knot) (radius : RationalInterval) (a b tol : ℚ)
    (depth bits degree : Nat) : Option RationalInterval := do
  guard (accepts mode rows a b && decide (0 ≤ radius.lo ∧ radius.lo ≤ radius.hi))
  let S ← sumCells? j mode power radius a b depth bits degree (cells rows)
  finish? S tol bits degree

theorem integrate_of_parts {j mode power rows radius a b tol depth bits degree S I}
    (ha : (accepts mode rows a b && decide (0 ≤ radius.lo ∧ radius.lo ≤ radius.hi)) = true)
    (hs : sumCells? j mode power radius a b depth bits degree (cells rows) = some S)
    (hI : finish? S tol bits degree = some I) :
    integrate? j mode power rows radius a b tol depth bits degree = some I := by
  unfold integrate?
  rw [ha]
  simpa [hs] using hI

theorem integrate_sound {j mode power rows radius a b tol depth bits degree I} {Y : Knot → ℝ} {v : ℝ}
    (hR : radius.Contains v) (hp : ∀ k B, power k = some B → 0 ≤ Y k ∧ B.Contains (Y k))
    (h : integrate? j mode power rows radius a b tol depth bits degree = some I) :
    I.Contains (sumCells j mode Y v a b (cells rows) / (2*π^2)) ∧ I.hi - I.lo ≤ tol := by
  unfold integrate? at h
  obtain ⟨_, _, hrest⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨S, hS, hrest⟩ := Option.bind_eq_some_iff.mp hrest
  unfold finish? at hrest
  obtain ⟨J, hJ, hw⟩ := Option.bind_eq_some_iff.mp hrest
  obtain ⟨rfl,hw'⟩ := SashimiNumerics.within_sound hw
  have hs := sumCells_sound hR hp hS
  have hj := Expr.contains_gridEnclosureWith? (values := fun _ => sumCells j mode Y v a b (cells rows))
    (h := hJ) (by intro i U hU; cases i <;> simp_all)
  change J.Contains (_ / (2 * (4 * arctan 1)^2)) at hj
  rw [arctan_one, show 4*(π/4)=π by ring] at hj
  exact ⟨hj, hw'⟩

/-- Positive observed ordinates with an independently verified nonnegative transfer. -/
def weighted? (ratio : ℚ → Option RationalInterval) (k : Knot) : Option RationalInterval :=
  if 0 ≤ k.power then (ratio k.k).map (·.scaleNonnegative k.power) else none

theorem weighted_sound {ratio : ℚ → Option RationalInterval} {T : ℝ → ℝ}
    (hT : ∀ k J, ratio k = some J → 0 ≤ T k ∧ J.Contains (T k))
    (k : Knot) (I : RationalInterval) (h : weighted? ratio k = some I) :
    0 ≤ (k.power : ℝ) * T k.k ∧ I.Contains ((k.power:ℝ) * T k.k) := by
  unfold weighted? at h
  split at h
  · rename_i hk
    obtain ⟨U,hu,rfl⟩ := Option.map_eq_some_iff.mp h
    obtain ⟨ht, hbox⟩ := hT _ _ hu
    exact ⟨mul_nonneg (by exact_mod_cast hk) ht, RationalInterval.contains_scaleNonnegative hbox hk⟩
  · contradiction

/-- Variance of SASHIMI-F's tabulated FDM power multiplier, applied at the knots
before positive-k log-log interpolation. -/
theorem fdm_integrate_sound {m rows radius a b tol depth bits degree I} {v : ℝ}
    (hR : radius.Contains v)
    (h : integrate? 0 .logLog (weighted? (fun k => fdm? m k bits degree)) rows radius a b tol depth bits degree = some I) :
    I.Contains (sumCells 0 .logLog (fun k => (k.power:ℝ) * fdmPowerRatio m k.k) v a b (cells rows) / (2*π^2)) ∧
      I.hi - I.lo ≤ tol :=
  integrate_sound hR (weighted_sound (fun _ _ h => ⟨sq_nonneg _, fdm_sound h⟩)) h

#print axioms fdm_integrate_sound
#print axioms integrate_sound

end SashimiTableWindow
