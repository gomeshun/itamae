import SashimiMassCalculus
import SashimiEisensteinHu

/-! Native binary64 quadrature for the baryonic Eisenstein–Hu comparison.
The caller supplies the signed transfer amplitude in h/Mpc. Sigma8 fixes
one common amplitude for the variance and its physical-mass derivative. -/
namespace SashimiEhVariance

structure Sample where
  k : Float
  weight : Float

structure Moment where
  raw : Float
  logSlopeNumerator : Float

/-- Fourier top hat and x*dW/dx. Their regular Taylor expansions avoid
cancellation about x=0. For x≥0.1 use x W'(x) = 3 sin(x)/x - 3 W(x). -/
def window (x : Float) : Float × Float :=
  if x < 0.1 then
    let u := x*x
    let w := 1 + u*(-0.1 + u*(1/280 + u*(-1/15120 + u*(1/1330560 - u/172972800))))
    let dw := u*(-0.2 + u*(1/70 + u*(-1/2520 + u*(1/166320 - u/17297280))))
    (w,dw)
  else
    let s := x.sin / x
    let w := 3 * (s - x.cos) / (x*x)
    (w, 3*s - 3*w)

/-- Composite Simpson in ln k; include dk/k in each spectral weight. -/
def nodes (transfer : Float → Float) (ns : Float) (logLo logHi : Float)
    (panels : Nat) : IO (Array Sample) := do
  unless 0 < panels && panels % 2 == 0 && logLo < logHi do
    throw (IO.userError "invalid Simpson mesh")
  let step := (logHi - logLo) / panels.toFloat
  let pi := Float.acos (-1)
  let coeff := step / (6*pi*pi)
  let mut data := Array.emptyWithCapacity (panels + 1)
  for i in [0:panels+1] do
    let logK := logLo + i.toFloat*step
    let k := logK.exp
    let t := transfer k
    let power := ((ns+3)*logK).exp * t*t
    let w := if i == 0 || i == panels then 1 else if i%2 == 0 then 2 else 4
    let weight := w * coeff * power
    unless weight.isFinite && weight ≥ 0 do throw (IO.userError "nonfinite spectral weight")
    data := data.push ⟨k,weight⟩
  return data

/-- Compensated sums for S and its signed derivative with respect to ln M.
dW(kR)/d ln M = x W'(x)/3, with x=kR and R ∝ M^(1/3). -/
def integrate (samples : Array Sample) (R : Float) : Moment := Id.run do
  let mut sum := 0.0
  let mut slope := 0.0
  let mut comp := 0.0
  let mut gradComp := 0.0
  for p in samples do
    let (w,dw) := window (p.k*R)
    let base := p.weight * w
    let y := base*w - comp
    let next := sum+y
    comp := (next-sum)-y
    sum := next
    let dy := (2/3 : Float)*base*dw - gradComp
    let nextSlope := slope+dy
    gradComp := (nextSlope-slope)-dy
    slope := nextSlope
  return ⟨sum,slope⟩

/-- R in Mpc/h for a physical mass in solar masses. The density constant
is Colossus 1.4.0's tabulated critical density with its h² factored out. -/
def radius (mass h omegaM : Float) : Float :=
  (3*mass*h/(4*Float.acos (-1)*277536627245.708*omegaM)).pow (1/3)

/-- Parent-history mass at which S enters the native SASHIMI-C fit. -/
def progenitor (host : Float) : Float × Float × Float :=
  let l := Float.log10 host
  let zf := -0.0064*l*l + 0.0237*l + 1.8837
  let q := 4.137/(zf.pow 0.9476)
  (host/q,q,zf)

/-- Exact real spectrum used in the normalization and signed derivative
below. The amplitude is fixed only after integrating the R=8 window. -/
noncomputable def realPower (k : ℝ) : ℝ :=
  k ^ (0.965 : ℝ) *
    (SashimiEisensteinHu.transferReal
      (SashimiEisensteinHu.prepareReal 0.674 0.315 0.049 2.7255) k)^2

/-- ρcrit Ωm / h compensates for input mass in physical solar masses
when the radius and wavenumber use Mpc/h and h/Mpc. -/
noncomputable def massDensity : ℝ := 277536627245.708 * 0.315 / 0.674

theorem radius_identity (M : ℝ) :
    SashimiSpectrum.massRadius massDensity M =
      SashimiSpectrum.massRadius (277536627245.708*0.315) (0.674*M) := by
  unfold SashimiSpectrum.massRadius massDensity
  congr 1; ring

noncomputable def rawAtRadius (r a b : ℝ) : ℝ :=
  ∫ k in a..b, k^2 * realPower k * (SashimiSpectrum.topHat (k*r))^2 / (2*Real.pi^2)

noncomputable def normalizedS (M a b : ℝ) : ℝ :=
  (0.811^2 / rawAtRadius 8 a b) * rawAtRadius (SashimiSpectrum.massRadius massDensity M) a b

noncomputable def normalizedDerivative (M a b : ℝ) : ℝ :=
  (0.811^2 / rawAtRadius 8 a b) *
    ∫ k in a..b, SashimiMassVariance.kernel .topHatDerivative realPower massDensity M 1 k

private theorem rawAtMass (M a b : ℝ) :
    rawAtRadius (SashimiSpectrum.massRadius massDensity M) a b =
      ∫ k in a..b, SashimiMassVariance.kernel .topHat realPower massDensity M 1 k := by
  unfold rawAtRadius
  apply intervalIntegral.integral_congr
  intro k _
  unfold SashimiMassVariance.kernel
  ring

/-- For the exact model on a positive finite band, differentiating the
window yields the derivative of the same σ₈-normalized S. The proof uses
continuity only on the requested band. -/
theorem massDerivative_spec {M a b : ℝ} (hm : 0 < M) (ha : 0 < a) (hab : a ≤ b)
    (hp : ContinuousOn realPower (Set.Icc a b)) :
    HasDerivAt (fun t => normalizedS t a b) (normalizedDerivative M a b) M := by
  have hr : 0 < massDensity := by norm_num [massDensity]
  have hh := (SashimiMassVariance.topHatCell_of_continuousOn realPower hr hm ha hab hp (1:ℝ)).const_mul
    (0.811^2 / rawAtRadius 8 a b)
  have heq : (fun t => normalizedS t a b) = fun t =>
      (0.811^2 / rawAtRadius 8 a b) * ∫ k in a..b, SashimiMassVariance.kernel .topHat realPower massDensity t 1 k := by
    funext t
    exact congrArg (fun y => (0.811^2 / rawAtRadius 8 a b) * y) (rawAtMass t a b)
  rw [heq]
  exact hh

/-- Return exact bit patterns of computed binary64 quantities to the reader. -/
def runMesh (transfer : Float → Float) (masses : Array Float) (low high : Int)
    (perDecade : Nat) : IO Unit := do
  let ns := 0.965
  let sig8 := 0.811
  let version := fun (x : Float) => x.toBits.toNat
  let mid := (high-low).toNat * perDecade
  let ln10 := Float.log 10
  let samples ← nodes transfer ns (low.toFloat*ln10) (high.toFloat*ln10) mid
  let anchor := integrate samples 8
  unless anchor.raw > 0 && anchor.raw.isFinite do throw (IO.userError "invalid R8 integral")
  let norm := sig8*sig8/anchor.raw
  let (extra,q,zf) := progenitor 1e12
  IO.println s!"anchor {low} {high} {perDecade} {version anchor.raw} {version norm} {version q} {version zf}"
  for mass in masses.push extra do
    let R := radius mass 0.674 0.315
    let i := integrate samples R
    let S := norm*i.raw
    let grad := norm*i.logSlopeNumerator/mass
    unless mass > 0 && S > 0 && S.isFinite && grad.isFinite do
      throw (IO.userError "invalid variance or gradient")
    IO.println s!"row {version mass} {version R} {version S.sqrt} {version S} {version grad} {version (i.logSlopeNumerator / i.raw)}"

end SashimiEhVariance
