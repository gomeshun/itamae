import SashimiEhSweep

/-! Composition of independently kernel-checked analytic EH transfer and moment blocks. -/
namespace SashimiEhCertificate

open FloatLib.Numerics SashimiEisensteinIntegral SashimiEhSweep SashimiSpectrum

theorem good_append {f a xs ys} (hx : Good f a xs) (hy : Good f (last a xs) ys) :
    Good f a (xs ++ ys) := by
  induction xs generalizing a with
  | nil => exact hy
  | cons x xs ih => exact ⟨hx.1,hx.2.1,hx.2.2.1,ih hx.2.2.2 hy⟩

theorem last_append (a : ℚ) (xs ys : List Bin) :
    last a (xs ++ ys) = last (last a xs) ys := by
  induction xs generalizing a with
  | nil => rfl
  | cons x xs ih => exact ih _

theorem total_append (R : RationalInterval) (bits degree : Nat) (a : ℚ) (xs ys : List Bin) :
    total R bits degree a (xs ++ ys) =
      (total R bits degree a xs).add (total R bits degree (last a xs) ys) := by
  induction xs generalizing a with
  | nil =>
    simp [total, totalBy, last, RationalInterval.point, RationalInterval.add_def]
  | cons x xs ih =>
    change (panel R bits degree a x).add (total R bits degree x.right (xs ++ ys)) = _
    rw [ih]
    change RationalInterval.add _ (RationalInterval.add _ _) = RationalInterval.add (RationalInterval.add _ _) _
    simp [RationalInterval.add_def, add_assoc, last, total]

/-- The exact finite-band numerator, the independent radius-8 denominator, and the
normalized variance all use the same certified continuum power on each cell. -/
theorem report {c m a b bins R N A bits degree}
    (hg : Good (weightReal c) a bins) (hEnd : last a bins = b)
    (hR : radiusBox? c m bits degree = some R) (hR0 : 0 ≤ R.lo)
    (hN : total R bits degree a bins = N)
    (hA : total (.point 8) bits degree a bins = A) (hpos : 0 < A.lo) :
    N.Contains (raw c a b (radius c m)) ∧ A.Contains (raw c a b 8) ∧
    ((N.scaleNonnegative (c.sigma8^2)).divPositive A).Contains
      (normalize8 c.sigma8 (raw c a b 8) (raw c a b (radius c m))) := by
  have hn := total_sound (bits := bits) (degree := degree) hg (radiusBox_sound hR) hR0
  have ha := total_sound (bits := bits) (degree := degree) hg
    (show (RationalInterval.point 8).Contains (8:ℝ) by norm_num [RationalInterval.Contains, RationalInterval.point])
    (by norm_num [RationalInterval.point] : 0 ≤ (RationalInterval.point 8).lo)
  rw [hEnd, hN] at hn
  rw [hEnd, hA] at ha
  refine ⟨hn,ha,?_⟩
  simpa only [normalize8, Rat.cast_pow] using RationalInterval.contains_divPositive
    (RationalInterval.contains_scaleNonnegative hn (sq_nonneg c.sigma8)) ha hpos

#print axioms report

end SashimiEhCertificate
