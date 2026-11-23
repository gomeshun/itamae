import SashimiTable

/-! Join independently checked portions of one tabulation. Consecutive portions
share their whole boundary knot, including its power ordinate. -/
namespace SashimiTable

private theorem cells_cons_cons (x y : Knot) (xs : List Knot) :
    cells (x :: y :: xs) = ⟨x,y⟩ :: cells (y :: xs) := rfl

private theorem cells_split (left : List Knot) (mid : Knot) (right : List Knot) :
    cells (left ++ mid :: right) = cells (left ++ [mid]) ++ cells (mid :: right) := by
  induction left with
  | nil => rfl
  | cons x xs ih =>
    cases xs with
    | nil => simp [cells]
    | cons y ys =>
      simp only [List.cons_append, cells_cons_cons] at ih ⊢
      rw [ih]

/-- Removing the repeated first knot of the right block preserves every
adjacent cell, including the first cell after the block boundary. -/
theorem cells_join {xs ys : List Knot} {mid : Knot}
    (hx : xs.getLast? = some mid) (hy : ys.head? = some mid) :
    cells (xs ++ ys.tail) = cells xs ++ cells ys := by
  obtain ⟨left,rfl⟩ := List.getLast?_eq_some_iff.mp hx
  obtain ⟨right,rfl⟩ := List.head?_eq_some_iff.mp hy
  simpa only [List.tail_cons, List.append_assoc, List.cons_append, List.nil_append]
    using cells_split left mid right

/-- A verified portion of the original knot list and its exact list of
adjacent cells. Endpoints are explicit so joining large blocks is symbolic. -/
structure Block (mode : Interpolation) (rows : List Knot) (cs : List Cell) (first last : Knot) : Prop where
  head : rows.head? = some first
  lastKnot : rows.getLast? = some last
  size : 2 ≤ rows.length
  fields : rows.all (fun r => decide (0 < r.k ∧ 0 ≤ r.power ∧ (mode = .logLog → 0 < r.power))) = true
  ordered : cs.all (fun c => decide (c.left.k < c.right.k)) = true
  adjacent : cells rows = cs

theorem Block.join {mode xs ys cs ds first mid last} (x : Block mode xs cs first mid)
    (y : Block mode ys ds mid last) : Block mode (xs ++ ys.tail) (cs ++ ds) first last := by
  constructor
  · rw [List.head?_append, x.head]; rfl
  · have hsize := y.size
    rw [List.getLast?_append, List.getLast?_tail, ite_eq_right (by omega), y.lastKnot]; rfl
  · rw [List.length_append]; exact x.size.trans (Nat.le_add_right _ _)
  · rw [List.all_append, x.fields]
    simp only [Bool.true_and, List.all_eq_true]
    intro r hr
    exact List.all_eq_true.mp y.fields r (List.mem_of_mem_tail hr)
  · rw [List.all_append, x.ordered, y.ordered]; rfl
  · rw [cells_join x.lastKnot y.head, x.adjacent, y.adjacent]

theorem Block.accepts {mode rows cs first last} (x : Block mode rows cs first last) {a b : ℚ}
    (hbound : first.k ≤ a ∧ a ≤ b ∧ b ≤ last.k) : SashimiTable.accepts mode rows a b = true := by
  unfold SashimiTable.accepts
  rw [x.head,x.lastKnot]
  simp only [x.fields, x.adjacent, x.ordered, Bool.and_true, decide_eq_true_eq]
  exact ⟨hbound.1,hbound.2.1,hbound.2.2,x.size⟩

#print axioms Block.join
#print axioms Block.accepts
end SashimiTable
