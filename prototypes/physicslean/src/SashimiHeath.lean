import SashimiNumerics
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.Analysis.SpecialFunctions.Pow.Deriv
import Mathlib.MeasureTheory.Integral.IntervalIntegral.IntegrationByParts

/-! Normalized linear growth for the flat dust and cosmological-constant Heath
model. Write x = √a and y = (1−m)x⁶/(m+(1−m)x⁶). The derivative of the truncated
primitive is exactly the positive integrand times (1−qₙ(y)). Monotonicity of
the residual then sandwiches the integral by two computable rationals.
The empirical collapse threshold 1.686 is a separate input convention. -/
namespace SashimiHeath
open Real Set MeasureTheory FloatLib.Numerics

def coefficient : Nat → ℚ
  | 0 => 1
  | n+1 => coefficient n * ((6*n+9:Nat) : ℚ)/(6*n+11:Nat)

def poly {K : Type*} [Field K] [CharZero K] : Nat → K → K
  | 0, _ => 1
  | n+1, y => poly n y + (coefficient (n+1):K) * y^(n+1)

private noncomputable def slope : Nat → ℝ → ℝ
  | 0, _ => 0
  | n+1, y => slope n y + (coefficient (n+1):ℝ) * ((n:ℝ)+1) * y^n

private theorem coefficient_step (n : ℕ) :
    (coefficient (n+1):ℝ) * (6*(n:ℝ)+11) = (coefficient n:ℝ) * (6*(n:ℝ)+9) := by
  simp [coefficient]; field_simp

private theorem poly_deriv (n : Nat) (y : ℝ) : HasDerivAt (poly n) (slope n y) y := by
  induction n with
  | zero => simpa only [poly, slope] using hasDerivAt_const y (1:ℝ)
  | succ n ih =>
    convert!
      ih.add (((hasDerivAt_id' y).pow (n+1)).const_mul (coefficient (n+1):ℝ)) using 1
    simp [slope]; ring

private theorem poly_identity (n : Nat) (y : ℝ) :
    (5-9*y)*poly n y + 6*y*(1-y)*slope n y =
      5 - (6*(n:ℝ)+9)*(coefficient n:ℝ)*y^(n+1) := by
  induction n with
  | zero => simp [poly,slope,coefficient]
  | succ n ih =>
    simp only [poly,slope, Nat.cast_add, Nat.cast_one, pow_succ]
    linear_combination ih + (y^n*y) * coefficient_step n

private theorem coefficient_pos (n : Nat) : 0 < coefficient n := by
  induction n with
  | zero => decide
  | succ n ih => unfold coefficient; positivity

private theorem poly_cast (n : Nat) (y : ℚ) : ((poly n y : ℚ) : ℝ) = poly n (y:ℝ) := by
  induction n <;> simp [poly, *]

def rad {K : Type*} [Field K] (m x : K) := m + (1-m)*x^6
def fraction {K : Type*} [Field K] (m x : K) := (1-m)*x^6 / rad m x
noncomputable def integrand (m x : ℝ) := 2*x^4/(rad m x * sqrt (rad m x))
noncomputable def primitive (n : Nat) (m x : ℝ) := (2/5)* x^5/(rad m x*sqrt (rad m x))*poly n (fraction m x)
def tail {K : Type*} [Field K] [CharZero K] (n : Nat) (y : K) :=
  (6*(n:K)+9)/5*(coefficient n:K)*y^(n+1)

private theorem rad_pos {m : ℝ} (hm : 0 < m) (h1 : m ≤ 1) (x : ℝ) : 0 < rad m x :=
  add_pos_of_pos_of_nonneg hm (mul_nonneg (sub_nonneg.mpr h1) (by positivity))

set_option maxHeartbeats 1800000 in
private theorem prim_deriv (n : Nat) {m : ℝ} (hm : 0 < m) (h1 : m ≤ 1) (x : ℝ) :
    HasDerivAt (primitive n m) (integrand m x * (1-tail n (fraction m x))) x := by
  let R := rad m x
  let L := 1-m
  let y := fraction m x
  let W := sqrt R
  have hR : 0 < R := rad_pos hm h1 x
  have hW : 0 < W := sqrt_pos.2 hR
  have hw2 : W^2 = R := sq_sqrt hR.le
  have hRAD : HasDerivAt (rad m) (6*L*x^5) x := by
    convert (((hasDerivAt_id' x).pow 6).const_mul L).const_add m using 1
    · rfl
    · dsimp; ring
  have hy : HasDerivAt (fraction m) (6*L*x^5/R*(1-y)) x := by
    convert ((((hasDerivAt_id' x).pow 6).const_mul L).div hRAD hR.ne') using 1
    · rfl
    · dsimp [fraction,y,R,L]; field_simp [(rad_pos hm h1 x).ne']
  have hxy : L*x^6=y*R := by dsimp [y,fraction,R,L]; field_simp [(rad_pos hm h1 x).ne']
  have hb : HasDerivAt (fun x => (2/5)*x^5/(rad m x * sqrt (rad m x)))
      (2*x^4/(R*W)*(1-(9/5)*y)) x := by
    have der := (((hasDerivAt_id' x).pow 5).const_mul (2/5:ℝ)).div
      (hRAD.mul (hRAD.sqrt hR.ne')) (mul_ne_zero hR.ne' hW.ne')
    apply der.congr_deriv
    dsimp
    change ((2/5)*(5*x^4*1)*(R*W)-(2/5*x^5)*(6*L*x^5*W+R*(6*L*x^5/(2*W))))/(R*W)^2 = _
    field_simp [hR.ne', hW.ne']
    grind only
  have hd := hb.mul ((poly_deriv n y).comp x hy)
  apply hd.congr_deriv
  change (2*x^4/(R*W)*(1-(9/5)*y))*poly n y +
    (2/5*x^5/(R*W))*(slope n y*(6*L*x^5/R*(1-y))) =
    (2*x^4/(R*W))*(1-tail n y)
  have hs : (2/5*x^5/(R*W))*(slope n y*(6*L*x^5/R*(1-y))) =
      2*x^4/(R*W)*(6/5*y*(1-y)*slope n y) := by
    field_simp [hR.ne', hW.ne']
    grind only
  rw [hs]
  have hp := poly_identity n y
  unfold tail
  linear_combination (2*x^4/(R*W)/5)*hp

noncomputable def heathIntegral (m x : ℝ) := ∫ t in 0..x, integrand m t

private theorem integral_cont (m : ℝ) (hm : 0 < m) (h1 : m ≤ 1) : Continuous (integrand m) := by
  unfold integrand
  have hr : Continuous (rad m : ℝ → ℝ) := by unfold rad; fun_prop
  exact ((continuous_id.pow 4).const_mul 2).div (hr.mul hr.sqrt)
    (fun x => mul_ne_zero (rad_pos hm h1 x).ne' (sqrt_pos.2 (rad_pos hm h1 x)).ne')

private theorem fraction_mono {m u v : ℝ} (hm : 0 < m) (h1 : m ≤ 1) (hu : 0 ≤ u) (huv : u ≤ v) :
    fraction m u ≤ fraction m v := by
  unfold fraction
  rw [div_le_div_iff₀ (rad_pos hm h1 u) (rad_pos hm h1 v)]
  have hh := mul_nonneg (mul_nonneg hm.le (sub_nonneg.mpr h1))
    (sub_nonneg.mpr (pow_le_pow_left₀ hu huv 6))
  unfold rad
  nlinarith only [hh]

private theorem integral_bounds (n : Nat) (m x : ℝ) (hm : 0 < m) (h1 : m ≤ 1) (hx : 0 ≤ x) :
    (1-tail n (fraction m x))*heathIntegral m x ≤ primitive n m x ∧
      primitive n m x ≤ heathIntegral m x := by
  let f := integrand m
  let Q := tail n (fraction m x)
  let g := fun t => f t * (1-tail n (fraction m t))
  have hr : Continuous (rad m : ℝ → ℝ) := by unfold rad; fun_prop
  have hy : Continuous (fraction m : ℝ → ℝ) := ((continuous_id.pow 6).const_mul _).div hr (fun t => (rad_pos hm h1 t).ne')
  have hi := integral_cont m hm h1
  have hg : Continuous g := by
    dsimp [g,f,tail]
    fun_prop
  have he := intervalIntegral.integral_eq_sub_of_hasDerivAt
    (f := primitive n m) (f' := g) (a := 0) (b := x) (fun t _ => prim_deriv n hm h1 t)
    (hg.intervalIntegrable 0 x)
  have hz : primitive n m 0 = 0 := by simp [primitive]
  rw [hz,sub_zero] at he
  have hf (t : ℝ) : 0 ≤ f t := div_nonneg (by positivity) (mul_nonneg (rad_pos hm h1 t).le (sqrt_nonneg _))
  have hy0 (t : ℝ) : 0 ≤ fraction m t :=
    div_nonneg (mul_nonneg (sub_nonneg.mpr h1) (by positivity)) (rad_pos hm h1 t).le
  have hc : 0 ≤ (coefficient n : ℝ) := by exact_mod_cast (coefficient_pos n).le
  have ht (t : ℝ) : 0 ≤ tail n (fraction m t) := by unfold tail; positivity [hy0 t]
  have hbound (t : ℝ) (h : t ∈ Icc 0 x) : tail n (fraction m t) ≤ Q := by
    dsimp [Q,tail]
    gcongr
    · exact hy0 t
    · exact fraction_mono hm h1 h.1 h.2
  rw [← he]
  constructor
  · rw [heathIntegral, ← intervalIntegral.integral_const_mul]
    apply intervalIntegral.integral_mono_on hx ((hi.const_mul _).intervalIntegrable _ _) (hg.intervalIntegrable _ _)
    intro t h
    dsimp [g]
    rw [mul_comm (1-_)]
    exact mul_le_mul_of_nonneg_left (sub_le_sub_left (hbound t h) 1) (hf t)
  · apply intervalIntegral.integral_mono_on hx (hg.intervalIntegrable _ _) (hi.intervalIntegrable _ _)
    intro t _
    exact mul_le_of_le_one_right (hf t) (by linarith [ht t])

/-- a is the scale factor in units where a₀ = 1. The change of variable a'=x²
removes the endpoint singularity of E(a) inside the traditional Heath integral. -/
noncomputable def growingX (m x : ℝ) := sqrt (rad m x) / x^3 * heathIntegral m x
noncomputable def growing (m a : ℝ) := growingX m (sqrt a)
def density {K : Type*} [Field K] (m a : K) := m + (1-m)*a^3
def yAt {K : Type*} [Field K] (m a : K) := (1-m)*a^3 / density m a
def approx {K : Type*} [Field K] [CharZero K] (n : Nat) (m a : K) := (2/5)*a/density m a*poly n (yAt m a)

private theorem growth_bounds (n : Nat) {m a : ℝ} (hm : 0 < m) (h1 : m ≤ 1) (ha : 0 < a) :
    (1-tail n (yAt m a))*growing m a ≤ approx n m a ∧ approx n m a ≤ growing m a := by
  let x := sqrt a
  let E := sqrt (rad m x)/x^3
  have hx : 0 < x := sqrt_pos.2 ha
  have hr := rad_pos hm h1 x
  have he : 0 ≤ E := div_nonneg (sqrt_nonneg _) (pow_nonneg hx.le _)
  have hx2 : x^2 = a := sq_sqrt ha.le
  have hx6 : x^6 = a^3 := by rw [← hx2]; ring
  have hR : rad m x = density m a := by simp [rad, density, hx6]
  have hy : fraction m x = yAt m a := by simp [fraction,yAt,hR,hx6]
  have heq : E * primitive n m x = approx n m a := by
    dsimp [E,primitive,approx]; rw [← hR, ← hy, ← hx2]
    field_simp [hx.ne',hr.ne',(sqrt_pos.2 hr).ne']
  have hh := integral_bounds n m x hm h1 hx.le
  change _*(E*heathIntegral m x) ≤ _ ∧ _ ≤ E*heathIntegral m x
  rw [← heq, ← hy]
  refine ⟨?_, mul_le_mul_of_nonneg_left hh.2 he⟩
  simpa only [mul_left_comm] using mul_le_mul_of_nonneg_left hh.1 he
/-- Dimensionless E=H/H₀ for flat dust + Λ, with a normalized to one today. -/
noncomputable def expansion (m a : ℝ) := sqrt (density m a)/(a*sqrt a)

private noncomputable def scaleIntegrand (m a : ℝ) :=
  a*sqrt a / (density m a * sqrt (density m a))

set_option maxHeartbeats 800000 in
/-- Recover the usual scale-factor version of the Heath integral from its
regular parametrization. Its constant normalization cancels in D(z). -/
theorem growing_eq_heath {m a : ℝ} (hm : 0<m) (h1 : m≤1) (ha : 0<a) :
    growing m a = expansion m a * ∫ s in 0..a, 1/(s*expansion m s)^3 := by
  let x := sqrt a
  have hx : 0<x := sqrt_pos.2 ha
  have hx2 : x^2 = a := sq_sqrt ha.le
  have hr (s : ℝ) (hs : 0≤s) : 0<density m s := by unfold density; positivity [sub_nonneg.mpr h1]
  have hrad (t : ℝ) : density m (t^2) = rad m t := by unfold density rad; ring
  have hcp : Continuous (fun s : ℝ => s * sqrt s) := continuous_id.mul continuous_sqrt
  have hcd : Continuous (density m : ℝ → ℝ) := by unfold density; fun_prop
  have subst := intervalIntegral.integral_comp_mul_deriv' (a := 0) (b := x) (f := fun t : ℝ => t^2)
    (f' := fun t : ℝ => 2*t) (g := scaleIntegrand m)
    (fun t _ => by convert! (hasDerivAt_id' t).pow 2 using 1; simp)
    (by fun_prop) (by
      unfold scaleIntegrand
      apply hcp.continuousOn.div (hcd.mul hcd.sqrt).continuousOn
      rintro _ ⟨t,_,rfl⟩
      have h := hr (t^2) (sq_nonneg _)
      exact mul_ne_zero h.ne' (sqrt_pos.2 h).ne')
  simp only [zero_pow (by decide : 2 ≠ 0), hx2] at subst
  have hl : heathIntegral m x = ∫ t in (0:ℝ)..x, (scaleIntegrand m ∘ (fun v : ℝ => v^2)) t*(2*t) := by
    apply intervalIntegral.integral_congr
    intro t ht
    have ht' : 0 ≤ t := ((uIcc_of_le hx.le ▸ ht : t ∈ Icc 0 x)).1
    simp only [Function.comp_apply,scaleIntegrand, integrand, sqrt_sq ht', hrad]
    ring
  have hs : (∫ s in (0:ℝ)..a, scaleIntegrand m s) = ∫ s in 0..a, 1/(s*expansion m s)^3 := by
    apply intervalIntegral.integral_congr
    intro s hs
    have hs' : 0 ≤ s := ((uIcc_of_le ha.le ▸ hs : s ∈ Icc 0 a)).1
    rcases hs'.eq_or_lt with h|h
    · subst s; simp [scaleIntegrand,expansion]
    · have hW := sqrt_pos.2 (hr s h.le)
      have hS := sqrt_pos.2 h
      have hw2 := sq_sqrt (hr s h.le).le
      have hs2 := sq_sqrt h.le
      unfold scaleIntegrand expansion
      field_simp [h.ne',hS.ne',hW.ne',(hr s h.le).ne']
      grind only
  change sqrt (rad m x)/x^3 * heathIntegral m x = _
  rw [hl, subst, hs, expansion, ← hx2, hrad, sqrt_sq hx.le]
  ring

noncomputable def D (m z : ℝ) := growing m (1/(1+z)) / growing m 1
noncomputable def rate (m z : ℝ) :=
  let a := 1/(1+z)
  (3/2)*m*a/density m a * D m z - (a^2/density m a)/growing m 1
noncomputable def barrier (m z : ℝ) := (1686/1000) / D m z

private theorem growing_pos {m a : ℝ} (hm : 0<m) (h1 : m ≤ 1) (ha : 0<a) : 0 < growing m a := by
  apply lt_of_lt_of_le _ (growth_bounds 0 hm h1 ha).2
  simp only [approx, poly, mul_one, density]
  positivity [sub_nonneg.mpr h1]

private theorem growingX_deriv {m x : ℝ} (hm : 0<m) (h1 : m ≤ 1) (hx : 0<x) :
    HasDerivAt (growingX m) ((-3*m/(rad m x*x))*growingX m x + 2*x/rad m x) x := by
  let R := rad m x
  let W := sqrt R
  have hR : 0 < R := rad_pos hm h1 x
  have hW : 0 < W := sqrt_pos.2 hR
  have hw2 : W^2 = R := sq_sqrt hR.le
  have hr : HasDerivAt (rad m) (6*(1-m)*x^5) x := by
    convert! (((hasDerivAt_id' x).pow 6).const_mul (1-m)).const_add m using 1
    dsimp; ring
  have he : HasDerivAt (fun t => sqrt (rad m t) / t^3)
      ((-3*m/(R*x))*(W/x^3)) x := by
    apply ((hr.sqrt hR.ne').div ((hasDerivAt_id' x).pow 3) (pow_ne_zero 3 hx.ne')).congr_deriv
    change ((6*(1-m)*x^5/(2*W))*x^3-W*(3*x^2*1))/(x^3)^2=_
    field_simp [hR.ne', hW.ne', hx.ne']
    have hdef : R = m + (1-m)*x^6 := rfl
    grind only
  have hi := integral_cont m hm h1
  have hd : HasDerivAt (heathIntegral m) (integrand m x) x :=
    intervalIntegral.integral_hasDerivAt_right (hi.intervalIntegrable _ _)
      hi.aestronglyMeasurable.stronglyMeasurableAtFilter hi.continuousAt
  apply (he.mul hd).congr_deriv
  change ((-3*m/(R*x))*(W/x^3))*heathIntegral m x + W/x^3*(2*x^4/(R*W)) =
      -3*m/(R*x)*((W/x^3)*heathIntegral m x) + 2*x/R
  field_simp [hR.ne', hW.ne', hx.ne']

/-- Differentiate the normalized Heath growing mode at fixed cosmological density.
The term containing G(1) comes from differentiating the integral's moving endpoint. -/
theorem D_deriv {m z : ℝ} (hm : 0<m) (h1 : m ≤ 1) (hz : 0 ≤ z) :
    HasDerivAt (D m) (rate m z) z := by
  let a := 1/(1+z)
  let x := sqrt a
  have hZ : 0 < 1+z := by linarith
  have ha : 0 < a := one_div_pos.mpr hZ
  have hx : 0 < x := sqrt_pos.2 ha
  have hx2 : x^2 = a := sq_sqrt ha.le
  have hR := rad_pos hm h1 x
  have hG := growing_pos hm h1 (by norm_num : (0:ℝ)<1)
  have hd : density m a = rad m x := by unfold density rad; rw [← hx2]; ring
  have hdA : HasDerivAt (fun t : ℝ => 1/(1+t)) (-a^2) z := by
    convert! ((hasDerivAt_const z (1:ℝ)).div ((hasDerivAt_id' z).const_add 1) hZ.ne') using 1
    simp [a,one_div,inv_pow,neg_div]
  have hdX := hdA.sqrt ha.ne'
  apply (((growingX_deriv hm h1 hx).comp z hdX).div_const (growing m 1)).congr_deriv
  change ((-3*m/(rad m x*x)*growingX m x + 2*x/rad m x)*(-a^2/(2*x)))/(growing m 1) =
    3/2*m*a/density m a * (growingX m x / growing m 1) - (a^2/density m a)/growing m 1
  rw [hd, ← hx2]
  field_simp [hx.ne', hR.ne', hG.ne']; ring

/-- Fused recurrence (last term, partial sum), with the same rational arithmetic
as the polynomial and remainder proven against the integral. -/
def terms (y : ℚ) : Nat → ℚ × ℚ
  | 0 => (1,1)
  | n+1 => let (last,total) := terms y n
           let next := (last * ((6*n+9 : Nat):ℚ)/(6*n+11:Nat))*y
           (next,total+next)

private theorem terms_eq (y : ℚ) (n : ℕ) :
    terms y n = (coefficient n*y^n, poly n y) := by
  induction n with
  | zero => simp [terms,coefficient,poly]
  | succ n ih =>
    simp only [terms, ih, coefficient, poly, pow_succ, Rat.cast_id]
    congr 1 <;> ring

def factor (n : ℕ) (m a : ℚ) : RationalInterval :=
  let y := yAt m a
  let (last,total) := terms y n
  let t := (2/5)*a/density m a * total
  let q := ((6*n+9:Nat):ℚ)/5*last*y
  ⟨t, t/(1-q)⟩

def residual (n : ℕ) (m a : ℚ) : ℚ :=
  ((6*n+9:Nat):ℚ)/5*(terms (yAt m a) n).1*yAt m a

theorem factor_sound {n : ℕ} {m a : ℚ} (hm : 0<m) (h1 : m≤1) (ha : 0<a)
    (hq : residual n m a < 1) : (factor n m a).Contains (growing m a) := by
  have hm' : (0:ℝ)<m := by exact_mod_cast hm
  have h1' : (m:ℝ)≤1 := by exact_mod_cast h1
  have ha' : (0:ℝ)<a := by exact_mod_cast ha
  have har : ((approx n m a : ℚ):ℝ) = approx n (m:ℝ) (a:ℝ) := by
    simp [approx, poly_cast, yAt, density]
  have hr : ((residual n m a:ℚ):ℝ) = tail n (yAt (m:ℝ) (a:ℝ)) := by
    simp [residual, terms_eq, tail, yAt, density, pow_succ]; ring
  have hf : factor n m a = ⟨approx n m a, approx n m a/(1-residual n m a)⟩ := by
    simp [factor, approx, residual, terms_eq]
  obtain ⟨hl,hu⟩ := growth_bounds n hm' h1' ha'
  have hq' : (residual n m a:ℝ)<1 := by exact_mod_cast hq
  rw [← har, ← hr] at hl
  rw [← har] at hu
  rw [hf]
  refine ⟨hu, ?_⟩
  push_cast
  exact (le_div_iff₀ (sub_pos.mpr hq')).2 (by linarith only [hl, mul_comm (1-(residual n m a:ℝ)) (growing m a)])

structure Row where
  D : RationalInterval
  dD_dz : RationalInterval
  delta_c : RationalInterval
deriving DecidableEq, Repr

/-- An external generator may round endpoints outward by any amount. This check
verifies each endpoint using only certified rational arithmetic. -/
def check (n : Nat) (m z : ℚ) (out : Row) : Bool := Id.run do
  let a := 1/(1+z)
  let r := density m a
  let g := factor n m a
  let u := factor n m 1
  let d : RationalInterval := if z=0 then ⟨1,1⟩ else ⟨g.lo/u.hi,g.hi/u.lo⟩
  let alpha := (3/2)*m*a/r
  let beta := a^2/r
  return decide (0<m ∧ m≤1 ∧ 0≤z ∧ residual n m a < 1 ∧ residual n m 1 < 1 ∧
    0<u.lo ∧ 0<u.hi ∧ 0<d.lo ∧ 0<d.hi ∧ out.D.lo ≤ d.lo ∧ d.hi ≤ out.D.hi ∧
    out.dD_dz.lo ≤ alpha*d.lo-beta/u.lo ∧ alpha*d.hi-beta/u.hi ≤ out.dD_dz.hi ∧
    out.delta_c.lo ≤ (1686/1000)/d.hi ∧ (1686/1000)/d.lo ≤ out.delta_c.hi)

private theorem ratio_bounds {u g : RationalInterval} {x y : ℝ} (hx : g.Contains x)
    (hy : u.Contains y) (hu : 0<u.lo) (hU : 0<u.hi) (hX : 0≤x) :
    (⟨g.lo/u.hi, g.hi/u.lo⟩ : RationalInterval).Contains (x/y) := by
  have hul : (0:ℝ) < u.lo := by exact_mod_cast hu
  have huh : (0:ℝ) < u.hi := by exact_mod_cast hU
  have hY : 0<y := lt_of_lt_of_le hul hy.1
  constructor
  · push_cast
    apply (div_le_div_iff₀ huh hY).2
    exact (mul_le_mul_of_nonneg_right hx.1 hY.le).trans (mul_le_mul_of_nonneg_left hy.2 hX)
  · push_cast
    apply (div_le_div_iff₀ hY hul).2
    exact (mul_le_mul_of_nonneg_left hy.1 hX).trans (mul_le_mul_of_nonneg_right hx.2 hY.le)

theorem check_sound {n : Nat} {m z : ℚ} {out : Row} (h : check n m z out = true) :
    out.D.Contains (D m z) ∧ out.dD_dz.Contains (deriv (D m) z) ∧ out.delta_c.Contains (barrier m z) := by
  let a : ℚ := 1/(1+z)
  let r := density m a
  let g := factor n m a
  let u := factor n m 1
  let d : RationalInterval := if z=0 then ⟨1,1⟩ else ⟨g.lo/u.hi,g.hi/u.lo⟩
  let alpha : ℚ := (3/2)*m*a/r
  let beta : ℚ := a^2/r
  change decide (_ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _ ∧ _) = true at h
  obtain ⟨hm,h1,hz,hq,hq0,hu,hU,hdlo,hdhi,hl,hh,hrl,hrh,hbl,hbh⟩ := of_decide_eq_true h
  change out.dD_dz.lo ≤ alpha*d.lo-beta/u.lo at hrl
  change alpha*d.hi-beta/u.hi ≤ out.dD_dz.hi at hrh
  change out.delta_c.lo ≤ (1686/1000)/d.hi at hbl
  change (1686/1000)/d.lo ≤ out.delta_c.hi at hbh
  change out.D.lo ≤ d.lo at hl
  change d.hi ≤ out.D.hi at hh
  have ha : 0<a := by dsimp [a]; positivity
  have hm' : (0:ℝ)<m := by exact_mod_cast hm
  have h1' : (m:ℝ)≤1 := by exact_mod_cast h1
  have hz' : (0:ℝ)≤z := by exact_mod_cast hz
  have ha' : (0:ℝ)<a := by exact_mod_cast ha
  have hac : (a:ℝ) = 1/(1+(z:ℝ)) := by simp [a]
  have hg : g.Contains (growing m a) := factor_sound hm h1 ha hq
  have hh0 : u.Contains (growing m 1) := by simpa [u] using factor_sound hm h1 (by decide : (0:ℚ)<1) hq0
  have hpos := growing_pos hm' h1' (by norm_num : (0:ℝ)<1)
  have hd : d.Contains (D m z) := by
    dsimp only [d]
    split
    · rename_i he
      simp [he, D, div_self hpos.ne', RationalInterval.Contains]
    · rw [D, ← hac]
      exact ratio_bounds hg hh0 hu hU (growing_pos hm' h1' ha').le
  have hdpos : (0:ℝ)<d.lo := by exact_mod_cast hdlo
  have hdpos' : (0:ℝ)<d.hi := by exact_mod_cast hdhi
  have hul : (0:ℝ)<u.lo := by exact_mod_cast hu
  have huh : (0:ℝ)<u.hi := by exact_mod_cast hU
  have hden : (0:ℚ)<r := by dsimp [r,density]; positivity [sub_nonneg.mpr h1]
  have hα : (0:ℝ)≤alpha := by exact_mod_cast (show (0:ℚ)≤alpha by dsimp [alpha]; positivity)
  have hβ : (0:ℝ)≤beta := by exact_mod_cast (show (0:ℚ)≤beta by dsimp [beta]; positivity)
  have heq : deriv (D m) z = (alpha:ℝ)*D m z-(beta:ℝ)/growing m 1 := by
    rw [(D_deriv hm' h1' hz').deriv]
    simp [rate,alpha,beta,r,density,hac]
  refine ⟨⟨(show (out.D.lo:ℝ) ≤ d.lo by exact_mod_cast hl).trans hd.1,
    hd.2.trans (by exact_mod_cast hh)⟩, ?_, ?_⟩
  · rw [heq]
    constructor
    · refine (show (out.dD_dz.lo:ℝ) ≤ (alpha:ℝ)*d.lo-(beta:ℝ)/u.lo by exact_mod_cast hrl).trans ?_
      exact sub_le_sub (mul_le_mul_of_nonneg_left hd.1 hα) (div_le_div_of_nonneg_left hβ hul hh0.1)
    · refine le_trans ?_ (show (alpha:ℝ)*d.hi-(beta:ℝ)/u.hi ≤ out.dD_dz.hi by exact_mod_cast hrh)
      exact sub_le_sub (mul_le_mul_of_nonneg_left hd.2 hα) (div_le_div_of_nonneg_left hβ hpos hh0.2)
  · have hp := lt_of_lt_of_le hdpos hd.1
    constructor
    · exact (show (out.delta_c.lo:ℝ)≤ (1686/1000)/d.hi by simpa only [Rat.cast_div,Rat.cast_ofNat] using (Rat.cast_le (K := ℝ)).2 hbl).trans
        (div_le_div_of_nonneg_left (by norm_num) hp hd.2)
    · exact (div_le_div_of_nonneg_left (by norm_num : (0:ℝ)≤1686/1000) hdpos hd.1).trans
        (show (1686/1000)/(d.lo:ℝ)≤out.delta_c.hi by simpa only [Rat.cast_div,Rat.cast_ofNat] using (Rat.cast_le (K := ℝ)).2 hbh)

#print axioms growing_eq_heath
#print axioms D_deriv
#print axioms check_sound

end SashimiHeath
