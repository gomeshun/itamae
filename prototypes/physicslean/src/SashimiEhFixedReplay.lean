import SashimiEhFixedCells
import CertifiedIntegralRatio

/-! Fresh block evaluation directly from rational edges. Independent kernel
witnesses for consecutive blocks compose into one continuous band integral. -/
namespace SashimiEhFixedReplay
open FloatLib.Numerics SashimiEhFixedCells SashimiEisensteinIntegral SashimiEhFullGrid
open MeasureTheory Set SashimiSpectrum

def block? (R : RationalInterval) (bits degree : Nat) (a : ℚ) (ends : List ℚ) :
    Option RationalInterval := do
  let boxes ← boxes? matched bits degree
  let cells ← build? matched boxes bits degree a ends
  return variance R bits degree a cells

theorem block_sound {R J : RationalInterval} {bits degree a ends} {x : ℝ}
    (hr : R.Contains x) (hp : 0 < R.lo)
    (h : block? R bits degree a ends = some J) : J.Contains (raw matched a (ends.getLastD a) x) := by
  obtain ⟨C,hC,h⟩ := Option.bind_eq_some_iff.mp h
  obtain ⟨qs,hqs,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  obtain ⟨hG,he⟩ := build_sound hC hqs
  simpa only [← he] using variance_sound (bits:=bits) (degree:=degree) hG hr hp

private theorem insideBand {r : ℝ} {a b : ℚ} (hr : 0 ≤ r)
    (ha : lowEdge ≤ a ∧ a ≤ highEdge) (hb : lowEdge ≤ b ∧ b ≤ highEdge) :
    IntervalIntegrable (fun k => SashimiEhSweep.weightReal matched k * topHat (k*r)^2)
      volume a b := by
  apply (matched_integrable r hr).mono_set
  exact uIcc_subset_uIcc (Icc_subset_uIcc ⟨by exact_mod_cast ha.1,by exact_mod_cast ha.2⟩)
    (Icc_subset_uIcc ⟨by exact_mod_cast hb.1,by exact_mod_cast hb.2⟩)

theorem join {r : ℝ} {a b c : ℚ} {I J : RationalInterval} (hr : 0 ≤ r)
    (ha : lowEdge ≤ a ∧ a ≤ highEdge) (hb : lowEdge ≤ b ∧ b ≤ highEdge)
    (hc : lowEdge ≤ c ∧ c ≤ highEdge)
    (hI : I.Contains (raw matched a b r)) (hJ : J.Contains (raw matched b c r)) :
    (I.add J).Contains (raw matched a c r) := by
  have h := RationalInterval.contains_add hI hJ
  have h1 := insideBand hr ha hb
  have h2 := insideBand hr hb hc
  change (I.add J).Contains (∫ k in (a:ℝ)..c, SashimiEhSweep.weightReal matched k * topHat (k*r)^2)
  rw [← intervalIntegral.integral_add_adjacent_intervals h1 h2]
  exact h

/-- Signed integral samples with two window functions and one positive weight.
Rounding the three input ranges once allows exact integer bisection on dyadics. -/
def shared (N A : ℚ → ℚ → RationalInterval) (bits : Nat) (q a : ℚ) :
    List Cell → RationalInterval
  | [] => .point 0
  | b :: bs => (CertifiedIntegralRatio.contribution q
      ⟨SashimiEhSweep.rounded bits b.integral, SashimiEhSweep.rounded bits (N a b.right),
        SashimiEhSweep.rounded bits (A a b.right)⟩).add (shared N A bits q b.right bs)

theorem shared_sound {w n d : ℝ → ℝ} {N A : ℚ → ℚ → RationalInterval} {a qs bits} {q : ℚ}
    (hw : Measurable w) (hf : Measurable (fun x => w x * (n x - q*d x)))
    (hnn : ∀ x, 0 ≤ w x) (hg : Good w a qs) (hq : 0 ≤ q)
    (hN : ∀ (l u : ℚ), 0 < l → l ≤ u → ∀ x ∈ Icc (l:ℝ) u, (N l u).Contains (n x))
    (hA : ∀ (l u : ℚ), 0 < l → l ≤ u → ∀ x ∈ Icc (l:ℝ) u, (A l u).Contains (d x)) :
    IntervalIntegrable (fun x => w x*(n x-q*d x)) volume a (lastEdge a qs) ∧
      (shared N A bits q a qs).Contains (∫ k in (a:ℝ)..(lastEdge a qs:ℝ), w k*(n k-q*d k)) := by
  induction qs generalizing a with
  | nil => simp [shared,lastEdge,RationalInterval.point,RationalInterval.Contains]
  | cons b bs ih =>
    obtain ⟨hv,ht⟩ := hg
    let D := (SashimiEhSweep.rounded bits (N a b.right)).sub
      ((SashimiEhSweep.rounded bits (A a b.right)).scaleNonnegative q)
    have hD (k : ℝ) (hk : k ∈ Icc (a:ℝ) b.right) : D.Contains (n k-q*d k) :=
      RationalInterval.contains_sub (SashimiEhSweep.round_sound _ (hN _ _ hv.pos hv.ordered k hk))
        (RationalInterval.contains_scaleNonnegative (SashimiEhSweep.round_sound _ (hA _ _ hv.pos hv.ordered k hk)) hq)
    have hwI := (SashimiEhSweep.cellIntegral_sound hw hv.ordered hv.range).1
    have hfI := (SashimiEhSweep.cellIntegral_sound (I := b.pointwise.mul D) hf hv.ordered
      (fun k hk => RationalInterval.contains_mul (hv.range k hk) (hD k hk))).1
    have hJ := CertifiedWeightedIntegral.weighted_sound (by exact_mod_cast hv.ordered) hwI hfI
      (fun x _ => hnn x) (SashimiEhSweep.round_sound bits hv.mass) hD
    obtain ⟨hi,hr⟩ := ih ht
    exact ⟨hfI.trans hi,(intervalIntegral.integral_add_adjacent_intervals hfI hi) ▸
      RationalInterval.contains_add hJ hr⟩

def signed (R : RationalInterval) (q : ℚ) (bits degree : Nat) :=
  shared (fun l u => SashimiTighterWindow.windowSq (l*R.lo) (u*R.hi) bits degree)
    (fun l u => SashimiTighterWindow.windowSq (l*8) (u*8) bits degree) bits q

theorem signed_sound {R : RationalInterval} {bits degree q a b qs} {x : ℝ}
    (hG : Good (SashimiEhSweep.weightReal matched) a qs) (he : lastEdge a qs = b)
    (ha : lowEdge ≤ a ∧ a ≤ highEdge) (hb : lowEdge ≤ b ∧ b ≤ highEdge)
    (hr : R.Contains x) (hp : 0 ≤ R.lo) (hq : 0 ≤ q) :
    (signed R q bits degree a qs).Contains (raw matched a b x - q*raw matched a b 8) := by
  let w := SashimiEhSweep.weightReal matched
  let n := fun k : ℝ => topHat (k*x)^2
  let d := fun k : ℝ => topHat (k*8)^2
  let f := fun k : ℝ => w k*(n k-q*d k)
  have hf : f = fun k => w k*n k-q*(w k*d k) := by funext k; dsimp [f]; ring
  have meas (r : ℝ) : Measurable (fun k => w k*topHat (k*r)^2) :=
    SashimiAnalyticWindow.measurableIntegrand (weighted matched)
      (fun i => ((coefficients matched).toList[i]?).getD 0) r
  have hmf : Measurable f := by rw [hf]; exact (meas x).sub (measurable_const.mul (meas 8))
  have hh := (shared_sound (n:=n) (d:=d) (q:=q) (bits:=bits)
    (N:=fun l u => SashimiTighterWindow.windowSq (l*R.lo) (u*R.hi) bits degree)
    (A:=fun l u => SashimiTighterWindow.windowSq (l*8) (u*8) bits degree)
    (SashimiEhDerivativeBounds.measurable_weight matched) hmf
    (fun k => by unfold SashimiEhSweep.weightReal; positivity) hG hq (by
      intro l u hl _ k hk
      have hl' : (0:ℝ) ≤ l := by exact_mod_cast hl.le
      have hr' : (0:ℝ) ≤ R.lo := by exact_mod_cast hp
      apply SashimiTighterWindow.windowSq_sound (mul_nonneg hl.le hp)
      · push_cast; exact mul_le_mul hk.1 hr.1 hr' (hl'.trans hk.1)
      · push_cast; exact mul_le_mul hk.2 hr.2 (hr'.trans hr.1) ((hl'.trans hk.1).trans hk.2)) (by
      intro l u hl _ k hk
      apply SashimiTighterWindow.windowSq_sound (mul_nonneg hl.le (by norm_num))
      · push_cast; gcongr; exact hk.1
      · push_cast; gcongr; exact hk.2)).2
  change (signed R q bits degree a qs).Contains (∫ k in (a:ℝ)..(lastEdge a qs : ℝ), f k) at hh
  rw [he,hf,intervalIntegral.integral_sub] at hh
  · rwa [intervalIntegral.integral_const_mul] at hh
  · exact insideBand ((by exact_mod_cast hp : (0:ℝ) ≤ R.lo).trans hr.1) ha hb
  · exact (insideBand (r:=8) (by norm_num) ha hb).const_mul _

/-- Closed physical target for rational sign tests on any certified fixed-cell
mesh spanning the matched band. The anchor needs only a positive lower bound. -/
theorem full_ratio {m l u : ℚ} {R A : RationalInterval} {bits degree qs}
    (hg : Good (SashimiEhSweep.weightReal matched) lowEdge qs) (he : lastEdge lowEdge qs = highEdge)
    (hr : R.Contains (radius matched m)) (hp : 0 < R.lo) (hk : 1 ≤ highEdge*R.lo)
    (hA : A.Contains (raw matched lowEdge highEdge 8)) (ha : 0 < A.lo)
    (hl : 0 ≤ l) (hu : 0 ≤ u)
    (hL : l*allowance 8 ≤ (signed R l bits degree lowEdge qs).lo)
    (hU : (signed R u bits degree lowEdge qs).hi + allowance R.lo ≤ 0) :
    ((⟨l,u⟩:RationalInterval).scaleNonnegative (matched.sigma8^2)).Contains (S m) := by
  have hband : lowEdge ≤ highEdge := by decide +kernel
  apply CertifiedIntegralRatio.full_sound hr.1 hp hk hA ha
    (signed_sound (bits:=bits) (degree:=degree) hg he ⟨le_rfl,hband⟩ ⟨hband,le_rfl⟩ hr hp.le hl)
    (signed_sound (bits:=bits) (degree:=degree) hg he ⟨le_rfl,hband⟩ ⟨hband,le_rfl⟩ hr hp.le hu)
  exact ⟨hl,hu,hL,hU⟩

end SashimiEhFixedReplay
