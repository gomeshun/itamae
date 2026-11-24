import SashimiGrowth
import SashimiSpectrumScalar
import FloatLib.Numerics.Enclosure.Rational.Converges

/-! SASHIMI-C non-cusp CDM mass variance. Mass is in solar masses (not
solar masses/h), and the native empirical pivot is 10¹⁰/h solar masses. -/
namespace EpsLudlow
open SashimiSpectrum Real FloatLib.Numerics FloatLib.Numerics.Interval
open scoped Topology

def pivot : ℚ := 10^10 / (0.674:ℚ)

structure Row where
  sigma : RationalInterval
  variance : RationalInterval
  derivative : RationalInterval
deriving DecidableEq

def growth? (z : ℚ) (bits degree : Nat) : Option RationalInterval :=
  if z = 0 then some (.point 1)
  else SashimiGrowth.enclose .growth 0.315 0.685 0.674 z 1 bits degree

theorem growth_sound {z bits degree G} (hg : growth? z bits degree = some G) :
    G.Contains (SashimiGrowth.growth 0.315 0.685 z) := by
  unfold growth? at hg
  split at hg
  · rename_i hz
    subst z; cases hg
    have hz := SashimiGrowth.growth_zero (m := (0.315:ℝ)) (by norm_num) (by norm_num)
    norm_num at hz ⊢
    rw [hz]
    exact RationalInterval.contains_point 1
  · have hh := (SashimiGrowth.enclose_sound hg).1
    simpa [SashimiGrowth.value] using hh

noncomputable def sigmaAt (z mass : ℝ) : ℝ :=
  SashimiGrowth.growth 0.315 0.685 z * ludlowSigma0 pivot mass

private theorem sigma0_positive (p m : ℝ) : 0 < ludlowSigma0 p m := by
  unfold ludlowSigma0 ludlowLogShape ludlowDen
  positivity

/-- The fit gives a logarithmic slope; transform it back to the derivative
of S with respect to physical mass, including the fixed-redshift D². -/
theorem variance_derivative (z m : ℝ) (hm : 0 < m) :
    deriv (fun x => (sigmaAt z x)^2) m =
    (sigmaAt z m)^2 * ludlowLogSlope (log (pivot:ℝ) - log m) * (2/m) := by
  let p : ℝ := pivot
  have hs := (((ludlow_log_slope p (log m)).exp).comp m (hasDerivAt_log hm.ne'))
  have ha : HasDerivAt (ludlowSigma0 p)
      (ludlowSigma0 p m * ludlowLogSlope (log p - log m) * m⁻¹) m := by
    convert hs.congr_of_eventuallyEq (f₁ := ludlowSigma0 p) (by
      filter_upwards [Ioi_mem_nhds hm] with x hx
      simp [exp_log (show 0 < x from hx), exp_log (sigma0_positive p x)]) using 1
    simp [exp_log hm, exp_log (sigma0_positive p m)]
  have hd := ((ha.const_mul (SashimiGrowth.growth 0.315 0.685 z)).pow 2).deriv
  change deriv (fun t => (SashimiGrowth.growth 0.315 0.685 z * ludlowSigma0 p t)^2) m = _
  refine hd.trans ?_
  unfold sigmaAt; dsimp [p]; ring

/-- Exact rational arithmetic propagates the two FloatLib enclosures and
the growth factor into σ, S and the signed derivative. -/
def evaluate? (mass z : ℚ) (bits degree : Nat) : Option Row := do
  guard (0 < mass ∧ 0 ≤ z)
  let D ← growth? z bits degree
  let s ← (ludlowExpr (.const pivot) (.const mass)).gridEnclosure? bits degree
  let l ← (ludlowSlopeExpr (.const pivot) (.const mass)).gridEnclosure? bits degree
  let sigma := D.mul s
  let S := sigma.mul sigma
  return ⟨sigma, S, (S.mul l).scaleNonnegative (2/mass)⟩

def Refines (mass z : ℝ) (r : Row) : Prop :=
  r.sigma.Contains (sigmaAt z mass) ∧ r.variance.Contains ((sigmaAt z mass)^2) ∧
    r.derivative.Contains (deriv (fun m => (sigmaAt z m)^2) mass)

theorem sound {mass z bits degree r} (h : evaluate? mass z bits degree = some r) :
    Refines mass z r := by
  unfold evaluate? at h
  obtain ⟨_,hg,h⟩ := Option.bind_eq_some_iff.mp h
  have hdom : 0 < mass ∧ 0 ≤ z := by by_contra hn; simp [guard, hn] at hg
  obtain ⟨D,hD,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨s,hs,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨l,hl,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  have h1 := Expr.contains_gridEnclosure? hs
  have h2 := Expr.contains_gridEnclosure? hl
  rw [eval_ludlowExpr] at h1
  rw [eval_ludlowSlopeExpr] at h2
  have g : (D.mul s).Contains (sigmaAt z mass) := RationalInterval.contains_mul (growth_sound hD) h1
  have hS := RationalInterval.contains_mul g g
  rw [← sq] at hS
  refine ⟨g,hS,?_⟩
  rw [variance_derivative _ _ (by exact_mod_cast hdom.1)]
  have hslope := RationalInterval.contains_mul hS h2
  have hscale := RationalInterval.contains_scaleNonnegative hslope (factor := 2/mass) (div_nonneg (by norm_num) hdom.1.le)
  simpa only [Expr.eval, Rat.cast_div, Rat.cast_ofNat, mul_comm ((2:ℝ)/(mass:ℝ))] using hscale

end EpsLudlow
