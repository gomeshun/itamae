import SashimiNormalization

/-! A single pure request interpreter for the executable and concrete kernel replays. -/
namespace SashimiTableRequest

open FloatLib.Numerics FloatLib.Numerics.Interval SashimiTable

/-- Columns in file order. WMAP7 puts P before k; CAMB puts k before P.
Their numerical units use h/Mpc and (Mpc/h)^3. Canonical input uses Mpc units. -/
def canonicalRows (layout : String) (h : ℚ) (rows : List Knot) : Option (List Knot) :=
  match layout with
  | "canonical" => if h = 1 then some rows else none
  | "camb" => if 0 < h then some (rows.map (fromHUnits h)) else none
  | "wmap7" => if 0 < h then some (rows.map fun r => fromHUnits h ⟨r.power,r.k⟩) else none
  | _ => none

/- Parse ASCII bases using reducible character operations. The native string
parsers in this Lean release do not reduce under `decide +kernel`. -/
private def digits (cs : List Char) : Option Nat := do
  guard (!cs.isEmpty)
  cs.foldlM (init := 0) fun n c => do
    guard ('0' ≤ c && c ≤ '9')
    return 10 * n + (c.toNat - '0'.toNat)

private def integer : List Char → Option Int
  | '-' :: cs => return -Int.ofNat (← digits cs)
  | cs => return Int.ofNat (← digits cs)

def rational (s : String) : Except String ℚ := do
  match s.toList.splitOn '/' with
  | [n, d] =>
    let some n := integer n | throw "invalid numerator"
    let some d := digits d | throw "invalid denominator"
    if d = 0 then throw "zero denominator"
    return mkRat n d
  | [n] =>
    let some n := integer n | throw "invalid integer"
    return n
  | _ => throw "invalid rational"

private def nat (s : String) : Except String Nat :=
  (digits s.toList).elim (throw "invalid natural number") pure

/-- Convert the W fitting scale in Mpc/h to Mpc before applying it to k in 1/Mpc. -/
def physicalAlphaExpr (m Ω h : ℚ) : Expr :=
  .binary .div (SashimiSpectrum.wdmAlphaExpr (.const m) (.const Ω) (.const h)) (.const h)

theorem eval_physicalAlphaExpr {m Ω h : ℚ} (hm : 0 < m) (hΩ : 0 < Ω) (hh : 0 < h)
    (env : Nat → ℝ) :
    (physicalAlphaExpr m Ω h).eval env = SashimiSpectrum.wdmAlpha m Ω h / h := by
  change (SashimiSpectrum.wdmAlphaExpr _ _ _).eval env / _ = _
  rw [SashimiSpectrum.eval_wdmAlphaExpr] <;>
    first | rfl | (change (0 : ℝ) < (_ : ℚ); exact_mod_cast ‹_›)

/-- A single physical scale, or the exact particle mass, abundance and h for a W fit. -/
private def parameter (model text : String) : Except String Expr := do
  match (text.toList.splitOn ',').map String.ofList with
  | [m, Ω, h] =>
    if model != "wdm-q5" && model != "wdm-q10" then throw "particle inputs require WDM"
    let m ← rational m
    let Ω ← rational Ω
    let h ← rational h
    if m ≤ 0 || Ω ≤ 0 || h ≤ 0 then throw "particle inputs must be positive"
    return physicalAlphaExpr m Ω h
  | [x] =>
    let x ← rational x
    if model != "cdm" && x ≤ 0 then throw "expected a positive transfer parameter"
    return .const x
  | _ => throw "invalid transfer parameter"

private structure Anchor where
  sigma : ℚ
  radius : ℚ
  low : ℚ
  high : ℚ
  depth : Nat

/-- Parse and evaluate all CLI arguments, returning an error for malformed requests.
A valid request can return none when its enclosure misses the requested precision. -/
def run (args : List String) (rows : List Knot) : Except String (Option RationalInterval) := do
  let [action, coords, place, model, param, radius, low, high, j, depth, bits, degree, budget] := args.take 13
    | throw "expected a spectrum-table request"
  let mode ← match coords with
    | "linear" => pure Interpolation.linear
    | "log-log" => pure Interpolation.logLog
    | _ => throw "unknown interpolation"
  let placement ← match place with
    | "continuous" => pure Placement.continuous
    | "at-knots" => pure Placement.atKnots
    | _ => throw "unknown transfer placement"
  let p ← parameter model param
  let radius ← rational radius
  let a ← rational low
  let b ← rational high
  let tol ← rational budget
  let depth ← nat depth
  let bits ← nat bits
  let degree ← nat degree
  let j ← nat j
  let ratio : Expr → Expr ← match model with
    | "cdm" => pure (fun _ => .const 1)
    | "wdm-q5" => pure (SashimiSpectrum.wdmPowerExpr .publishedQ5 p)
    | "wdm-q10" => pure (SashimiSpectrum.wdmPowerExpr .standardT2Q10 p)
    | "fdm" => pure (SashimiSpectrum.fdmPowerExpr p)
    | "gaussian" => pure (SashimiSpectrum.cutoffExpr p)
    | _ => throw "unknown power ratio"
  if radius < 0 || tol < 0 then throw "invalid physical inputs or tolerance"
  let anchor : Option Anchor ← match args.drop 13 with
    | [] => pure none
    | [s, r, l, u, d] => do
      let σ ← rational s
      let r ← rational r
      if σ ≤ 0 || r < 0 then throw "invalid sigma8 anchor"
      pure (some ⟨σ,r,← rational l,← rational u,← nat d⟩)
    | _ => throw "invalid sigma8 anchor arguments"
  let w : Expr → Expr := if radius = 0 then fun _ => .const 1
    else fun x => SashimiSpectrum.topHatExpr (x * .const radius)
  let f := moment j mode placement ratio w
  match action, anchor with
  | "power", none => return lookup? mode placement ratio rows a tol bits degree
  | "moment", none => return integrate? f mode rows a b tol depth bits degree
  | "power", some t =>
    return (SashimiNormalization.point? mode placement ratio rows a bits degree).bind fun I =>
      SashimiNormalization.normalized? t.sigma tol I mode placement ratio rows
        t.radius t.low t.high t.depth bits degree
  | "moment", some t =>
    return (SashimiNormalization.band? f mode rows a b depth bits degree).bind fun I =>
      SashimiNormalization.normalized? t.sigma tol I mode placement ratio rows
        t.radius t.low t.high t.depth bits degree
  | _, _ => throw "unknown table quantity"

end SashimiTableRequest
