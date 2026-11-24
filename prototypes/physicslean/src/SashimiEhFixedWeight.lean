import CertifiedWeightedIntegral
import SashimiEhDerivativeBounds
import SashimiEhTighterPower
import EhCenteredWeight
import SashimiEhMoment

/-! Integrate the common positive spectral weight once per original cell.
The mesh records pointwise and integrated bounds separately. -/
namespace SashimiEhFixedCells
open FloatLib.Numerics MeasureTheory Set SashimiEhSweep SashimiEisensteinIntegral

structure Cell where
  right : ℚ
  pointwise : RationalInterval
  integral : RationalInterval
deriving Repr, DecidableEq

def meet (I J : RationalInterval) : RationalInterval := ⟨max I.lo J.lo, min I.hi J.hi⟩

theorem meet_sound {I J : RationalInterval} {x : ℝ} (hI : I.Contains x) (hJ : J.Contains x) :
    (meet I J).Contains x :=
  ⟨by simpa [meet] using max_le hI.1 hJ.1, by simpa [meet] using le_min hI.2 hJ.2⟩

def buildCell? (c : Cosmology) (boxes : List RationalInterval) (a b : ℚ) (bits degree : Nat) :
    Option Cell := do
  let P ← SashimiEhTighterPower.weight? c boxes a b bits degree
  let rectangle := P.scaleNonnegative (b-a)
  let fallback := fun () => match SashimiEhMoment.integral? c boxes a b bits degree with
    | some J => meet rectangle J
    | none => rectangle
  let A := match EhCenteredWeight.center? c boxes a b bits degree with
    | some (C,D) => meet rectangle (CertifiedCenteredIntegral.cell a b C D)
    | none => fallback ()
  return ⟨b,P,A⟩

/-- Whichever integration rule succeeds, its result refines the pointwise
rectangle on the very same cell, including cases where the bounds coincide. -/
theorem buildCell_refines {c boxes a b bits degree q}
    (h : buildCell? c boxes a b bits degree = some q) :
    (q.pointwise.scaleNonnegative (b-a)).lo ≤ q.integral.lo ∧
    q.integral.hi ≤ (q.pointwise.scaleNonnegative (b-a)).hi := by
  unfold buildCell? at h
  obtain ⟨P,_,h⟩ := Option.bind_eq_some_iff.mp h
  cases h
  dsimp only
  cases EhCenteredWeight.center? c boxes a b bits degree with
  | some cd => cases cd; exact ⟨le_max_left .., min_le_left ..⟩
  | none =>
    cases SashimiEhMoment.integral? c boxes a b bits degree with
    | some _ => exact ⟨le_max_left .., min_le_left ..⟩
    | none => exact ⟨le_rfl,le_rfl⟩

structure Valid (w : ℝ → ℝ) (a : ℚ) (q : Cell) : Prop where
  pos : 0 < a
  ordered : a ≤ q.right
  range : ∀ k ∈ Icc (a:ℝ) q.right, q.pointwise.Contains (w k)
  mass : q.integral.Contains (∫ k in (a:ℝ)..q.right, w k)

theorem buildCell_sound {c boxes cb cd a b bits degree q} (hb : boxes? c cb cd = some boxes)
    (hq : buildCell? c boxes a b bits degree = some q) :
    q.right = b ∧ Valid (weightReal c) a q := by
  unfold buildCell? at hq
  obtain ⟨P,hP,hq⟩ := Option.bind_eq_some_iff.mp hq
  cases hq
  obtain ⟨ha,hab,hr⟩ := SashimiEhTighterPower.weight_sound hb hP
  have hm := SashimiEhDerivativeBounds.measurable_weight c
  have hi := (cellIntegral_sound hm hab hr).2
  refine ⟨rfl,ha,hab,hr,?_⟩
  dsimp only
  generalize hx : EhCenteredWeight.center? c boxes a b bits degree = cd
  cases cd with
  | some cd =>
    obtain ⟨C,D⟩ := cd
    apply meet_sound hi
    exact CertifiedWeightedIntegral.affine_band_sound (f := weightReal c) (a := a) (b := b)
      (C := C) (D := D) hab hm
      (EhCenteredWeight.center_sound (c := c) (boxes := boxes) (cb := cb) (cd := cd)
        (bits := bits) (degree := degree) (l := a) (u := b) (CD := (C,D)) hb hx)
  | none =>
    generalize hmom : SashimiEhMoment.integral? c boxes a b bits degree = mom
    cases mom with
    | some J => exact meet_sound hi (SashimiEhMoment.integral_sound hb hmom)
    | none => exact hi

def Good (w : ℝ → ℝ) (a : ℚ) : List Cell → Prop
  | [] => True
  | q::qs => Valid w a q ∧ Good w q.right qs

def lastEdge (a : ℚ) : List Cell → ℚ
  | [] => a
  | q::qs => lastEdge q.right qs

def build? (c : Cosmology) (boxes : List RationalInterval) (bits degree : Nat) (a : ℚ) :
    List ℚ → Option (List Cell)
  | [] => some []
  | b::bs => do
    let q ← buildCell? c boxes a b bits degree
    return q :: (← build? c boxes bits degree b bs)

theorem build_sound {c boxes cb cd bits degree a ends qs} (hb : boxes? c cb cd = some boxes)
    (hq : build? c boxes bits degree a ends = some qs) :
    Good (weightReal c) a qs ∧ lastEdge a qs = ends.getLastD a := by
  induction ends generalizing a qs with
  | nil => cases hq; exact ⟨trivial,rfl⟩
  | cons b bs ih =>
    obtain ⟨q,h,h'⟩ := Option.bind_eq_some_iff.mp hq
    obtain ⟨rest,hh,hh'⟩ := Option.bind_eq_some_iff.mp h'
    cases hh'
    obtain ⟨he,hv⟩ := buildCell_sound hb h
    obtain ⟨hg,hl⟩ := ih hh
    exact ⟨⟨hv,he ▸ hg⟩,by change lastEdge q.right rest = _; rw [he,List.getLastD_cons]; exact hl⟩

#print axioms buildCell_refines
#print axioms buildCell_sound
end SashimiEhFixedCells
