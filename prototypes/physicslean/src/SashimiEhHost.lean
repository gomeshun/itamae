import SashimiEhFullGrid
import SashimiHalo
import SashimiHeath

/-! Matched host assembly at initial redshift zero. The Correa fitted progenitor
ratio is a real expression; the Fourier integral uses its exact mass radius. -/
namespace SashimiEhHost
open SashimiEisensteinIntegral SashimiEhFull SashimiEhFullGrid SashimiSpectrum
open FloatLib.Numerics FloatLib.Numerics.Interval

def host : ℚ := 10^12
def qExpr := SashimiHalo.formationMassRatio (.const host)
def childExpr := Expr.binary .div (.const host) qExpr
noncomputable def child := childExpr.eval (fun _ => 0)
noncomputable def childRadius := massRadius ((matched.criticalDensity*matched.om:ℚ):ℝ)
  (child*(matched.h:ℝ))

def childRadiusExpr := positivePowExpr (.binary .div
  (.const (3*matched.h) * childExpr) (.const (4*(matched.criticalDensity*matched.om)) * piExpr)) (.const (1/3))

theorem childRadiusExpr_eq : childRadiusExpr.eval (fun _ => 0) = childRadius := by
  have hx : 0 < child := by
    show 0 < (host:ℝ)/(SashimiHalo.formationMassRatio (.const host)).eval (fun _ => 0)
    change 0 < (host:ℝ)/((4.137:ℚ) / Real.exp _)
    have hh : (0:ℝ)<host := by norm_num [host]
    positivity
  have hp : 0 < ((Expr.binary .div (.const (3*matched.h) * childExpr)
    (.const (4*(matched.criticalDensity*matched.om)) * piExpr)).eval (fun _ => 0)) := by
    change 0 < ((3*matched.h:ℚ):ℝ)*child/(((4*(matched.criticalDensity*matched.om):ℚ):ℝ)*piExpr.eval _)
    rw [eval_piExpr]
    norm_num [matched]; positivity
  rw [childRadiusExpr,eval_positivePowExpr hp]
  change (((3*matched.h:ℚ):ℝ)*child/(((4*(matched.criticalDensity*matched.om):ℚ):ℝ)*piExpr.eval _))^(((1/3:ℚ):ℝ)) = _
  rw [eval_piExpr]
  norm_num [childRadius, massRadius, matched]
  ring_nf

noncomputable def deltaS : ℝ := normalize8 matched.sigma8 (fullRaw matched 8)
  (fullRaw matched childRadius - fullRaw matched (radius matched host))

/-- Enclose a difference after adding each tail to its own side. -/
def fromPaired (cr pr : ℚ) (difference anchor : RationalInterval) :=
  let fullDiff : RationalInterval := ⟨difference.lo-allowance pr,difference.hi+allowance cr⟩
  (fullDiff.scaleNonnegative (matched.sigma8^2)).divPositive (pad anchor (allowance 8))

theorem paired_sound {cr pr : ℚ} {I A : RationalInterval} (hcr : (cr:ℝ) ≤ childRadius)
    (hpr : (pr:ℝ) ≤ radius matched host)
    (hp : 0<cr ∧ 1≤highEdge*cr ∧ 0<pr ∧ 1≤highEdge*pr ∧ 0<A.lo)
    (hI : I.Contains (raw matched lowEdge highEdge childRadius - raw matched lowEdge highEdge (radius matched host)))
    (hA : A.Contains (raw matched lowEdge highEdge 8)) :
    (fromPaired cr pr I A).Contains deltaS := by
  obtain ⟨cp,ck,pp,pk,ap⟩ := hp
  have hc := fullRaw_tail hcr cp ck
  have hr := fullRaw_tail hpr pp pk
  have ha := fullRaw_sound (rmin:=8) (by norm_num) (by norm_num) (by norm_num [highEdge]) hA
  have hd : (⟨I.lo-allowance pr,I.hi+allowance cr⟩:RationalInterval).Contains
      (fullRaw matched childRadius-fullRaw matched (radius matched host)) := by
    constructor
    · change ((I.lo - allowance pr:ℚ):ℝ) ≤ _
      push_cast
      linarith [hI.1]
    · change _ ≤ ((I.hi + allowance cr:ℚ):ℝ)
      push_cast
      linarith [hI.2]
  simpa only [fromPaired,deltaS,normalize8,Rat.cast_pow] using
    RationalInterval.contains_divPositive (RationalInterval.contains_scaleNonnegative hd (sq_nonneg matched.sigma8)) ha ap

/-- The two inputs are D'(0) and S(M₀/q)−S(M₀). Both are dimensionless. -/
def expression (out : SashimiHalo.AssemblyOutput) (z : ℚ) : Expr :=
  let m := SashimiHalo.massHistory (.const host) (.const z) (.const 1) (.var 0) (.var 1)
  let matter := matched.om * (1+z)^3
  let delta := SashimiHalo.virialOverdensity (.const (matter/(1-matched.om+matter)))
  match out with
  | .m200 => m
  | .virialRate => SashimiHalo.massRate (SashimiHalo.virialMassFit m (SashimiHalo.correa200 m z) delta)
      (.const z) (.const 1) (.var 0) (.var 1)

noncomputable def target (out : SashimiHalo.AssemblyOutput) (z : ℚ) :=
  (expression out z).eval (fun | 0 => deriv (SashimiHeath.D matched.om) 0 | _ => deltaS)

theorem eval_sound {out z slope ds bits degree I} (hs : slope.Contains (deriv (SashimiHeath.D matched.om) 0))
    (hd : ds.Contains deltaS) (h : (expression out z).gridEnclosureWith? [slope,ds] bits degree = some I) :
    I.Contains (target out z) := by
  apply Expr.contains_gridEnclosureWith? (h:=h)
  intro i box hi
  cases i with
  | zero => simp at hi; subst box; exact hs
  | succ i => cases i with
    | zero => simp at hi; subst box; exact hd
    | succ i => simp at hi

def checkRadius (R : RationalInterval) (bits degree : Nat) : Bool :=
  (childRadiusExpr.gridEnclosure? bits degree).any fun J => decide (R.lo ≤ J.lo ∧ J.hi ≤ R.hi)

theorem radius_check_sound {R bits degree} (h : checkRadius R bits degree) :
    R.Contains childRadius := by
  unfold checkRadius at h
  cases he : childRadiusExpr.gridEnclosure? bits degree with
  | none => simp [he] at h
  | some J =>
    have hsubset : R.lo ≤ J.lo ∧ J.hi ≤ R.hi := by simpa [he] using h
    have hJ := (show J.Contains (childRadiusExpr.eval (fun _ => 0)) from Expr.contains_gridEnclosure? he)
    rw [childRadiusExpr_eq] at hJ
    exact ⟨(by exact_mod_cast hsubset.1 : (R.lo:ℝ) ≤ J.lo).trans hJ.1,
      hJ.2.trans (by exact_mod_cast hsubset.2)⟩

#print axioms childRadiusExpr_eq
#print axioms paired_sound
#print axioms eval_sound
end SashimiEhHost
