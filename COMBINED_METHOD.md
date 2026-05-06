# Combined Method — implementation spec

Source: `Диссертация.doc` by А.Ю. Дёмин (МАИ, 2002), supervisor Г.Ф. Хахулин.
Title: «Комбинированный метод решения линейных задач целочисленного программирования
и его применение для оптимизации топологии локальных вычислительных сетей».

This document is a working translation/condensation of Chapter 2 (and the relevant
software-architecture bits of Chapter 3) sufficient to re-implement the method in
code. Section numbers below match the dissertation. Russian terms are kept in
parentheses where they are load-bearing.

> **Note on formulas.** The original `.doc` stores every equation as an embedded
> Microsoft Equation 3.0 OLE object (898 streams, MTEF binary). `textutil`,
> `pandoc` and Quick Look all drop them, leaving `EMBED Equation.3` placeholders
> in the converted text. The descriptive prose around each equation in the
> dissertation is, however, complete and unambiguous, so the formulas below have
> been reconstructed from prose + standard OR conventions. Where reconstruction
> is judgement-call rather than mechanical, it is flagged with **[recon]**.

---

## 1. Problem statement (§1.1, §2.2.1)

### 1.1 Original (user-facing) form — equation (2.1)

```
min   c · x
s.t.  A x  {≤,=,≥}  b           (2.1.1)  main constraints
      0 ≤ x ≤ d                  (2.1.2)  two-sided bounds
      x ∈ Z^n_+                  (2.1.3)  integrality
```

* `x ∈ Z^n_+` — `n`-vector of decision variables, non-negative integers.
* `c ∈ R^n` — row of objective coefficients, **arbitrary real** in the user form
  (sign is unrestricted; minimisation is forced internally — see §3).
* `A ∈ R^{m×n}`, `b ∈ R^m` — main constraints, mixed `≤`, `=`, `≥` allowed.
* `d ∈ Z^n_+` — per-variable upper bound. Every `xⱼ` has at least an integer upper
  bound (in network-topology problems it falls out of the formulation; in
  general it is supplied as an a-priori envelope).

The application context that drives the dissertation is the LAN-topology problem
(§1.1.1, §1.1.2): selecting a hierarchy of hubs/concentrators with hard limits on
ports per hub, total length of cabling, etc. That problem is purely
combinatorial after discretisation and motivates the size class the method
targets.

### 1.2 Internal canonical form — equation (2.17)

The Combined Method requires the problem in this shape:

```
min   c̃ · y
s.t.  Ã y ≤ b̃                    (2.17)
      0 ≤ y ≤ d̃
      y ∈ Z^n_+
      c̃ⱼ ≥ 0  for all j
      c̃₁ ≤ c̃₂ ≤ … ≤ c̃ₙ          (sorted ascending)
```

Reaching (2.17) from (2.1):

1. **Force minimisation.** If the original is `max c·x`, replace `c ← -c`.
2. **Force non-negative cost coefficients (2.4).** For every `j` with `c̃ⱼ < 0`,
   substitute `yⱼ := dⱼ − xⱼ` so the new coefficient is `−cⱼ ≥ 0`. The constant
   term that appears in the objective is dropped (it does not affect the
   argmin); the same substitution is performed in *every* main constraint.
   Equation (2.18): `xⱼ ↦ dⱼ − yⱼ`.
3. **Re-number variables ascending by `c̃ⱼ`** so `c̃₁ ≤ c̃₂ ≤ … ≤ c̃ₙ`. Required
   both by the lattice enumerator and by Stage 2 of the method (§2.2.3).
4. **Force `≤` constraints.**
   * `=`-rows are split into a `≤` and a `≥` pair.
   * `≥`-rows are multiplied by `−1` to become `≤`.
5. The resulting matrix and RHS are denoted `Ã ∈ R^{m̃×n}`, `b̃ ∈ R^{m̃}`
   (`m̃ ≥ m` after splitting equalities). The optimum of (2.17), if it exists,
   is denoted `q*`.

### 1.3 LP-companion form — equation (2.19)

For the simplex passes the same problem is taken to canonical LP form, but
**without** the `≥`→`≤` flip and **without** equality-splitting (those were only
needed for the lattice search). One slack `xₙ₊ᵢ` is added per main row and
minimisation is converted to maximisation:

```
max  −c̃' · x'
s.t.  Â x' = b̂             (2.19)
      0 ≤ x' ≤ d̂
      x' ∈ R^{n+m}_+        (continuous LP relaxation)
```

* `x' = (y, s)` where `y` are the structural vars (numbered identically to (2.17))
  and `s` are slacks introduced to turn rows into equalities.
* `c̃'` is `c̃` padded with zeros for slacks. Sign-flipped so simplex maximises.
  Hence `obj_value(2.17) = −obj_value(2.19)`.
* `Â` is built from the *original* `A` after the substitution (2.18) plus slack
  columns; the `≤`-only rewrite of step 4 is **not** applied here.
* `b̂` is `b` after substitution (2.18); `d̂` extends `d` with slack upper bounds
  (effectively `+∞` unless the slack maps to a `=` or two-sided row).
* The optimum of (2.19), when it exists, is denoted `q'`.

This dual representation — `(2.17)` for the lattice enumerator, `(2.19)` for the
simplex — is held in step throughout the algorithm.

---

## 2. Vector-lattice implicit enumeration (§2.2.2)

This is the integer-search primitive that the method calls in Stage 2 and
Stage 4. It is a specialised B&B over the integer lattice
`L = {y ∈ Z^n : 0 ≤ yⱼ ≤ dⱼ}`, exploring it as a tree by ±1 steps.

### 2.1 Lattice traversal

A *route* (маршрут) is a path through `L` formed by **forward** and **backward**
steps. A forward step on variable `j` (denoted "step forward by jⱼ") increments
`yⱼ` by 1 and moves to a deeper level. A backward step decrements `yⱼ` and
returns one level. The traversal starts and finishes at the origin `0`. It
finishes when, at the origin, every legal forward step has been explored.

To prevent revisiting the same point along different orderings, at every node
`y` the algorithm restricts forward steps to a subset `J(y)` of variable
indices:

```
J(y) = { j ∈ {1..n} :
          j ≥ j_last(y)          (lex/level guard, eq. 2.3)
          AND yⱼ < dⱼ }          (upper bound not hit)
```

`j_last(y)` is the index of the variable last incremented to reach `y` (the
"level number" — каждый шаг вперёд связан с `j`, see §2.2.2).
A point with `J(y) = ∅` is **terminal** (конечная); after entering a terminal
point on a forward step a backward step must follow.

### 2.2 First reduction rule — non-negative `c̃` (eq. 2.4)

Already enforced as part of the canonical-form step 1.2. With every
`c̃ⱼ ≥ 0`, every forward step strictly worsens the objective. Therefore:

> **Rule 1.** As soon as a forward step lands in a feasible point, do not take
> any further forward step from it — back up.

Corollary: the algorithm spends almost all its time in *infeasible* points;
feasible leaves end branches. If `0` itself is feasible, it is the global
optimum and the search is over.

### 2.3 Slacks for constraint tracking (eq. 2.5–2.6)

After the (2.17) rewrite all main constraints are `≤`. Introduce a slack `sᵢ`
per row, requiring `sᵢ ≥ 0`:

```
Σⱼ ãᵢⱼ yⱼ + sᵢ = b̃ᵢ                      (2.6)
sᵢ(y) := b̃ᵢ − Σⱼ ãᵢⱼ yⱼ
```

A point `y` is feasible for (2.17) iff `sᵢ(y) ≥ 0 ∀ i`. Each forward step on
`yⱼ` updates every slack: `sᵢ ← sᵢ − ãᵢⱼ`. For a `≤` row this means slack
*decreases* whenever `ãᵢⱼ > 0` and *increases* (helps feasibility) whenever
`ãᵢⱼ < 0`.

### 2.4 Infeasibility criterion **КН** (eqs. 2.7–2.9)

Applied at every newly entered point.

Two index sets at point `y`:

```
I⁻(y) = { i : sᵢ(y) < 0 }                         (2.7)  unsatisfied rows
J⁻ᵢ(y) = { j ∈ J(y) : ãᵢⱼ < 0 }                   (2.8)  forward steps that
                                                          increase row i's slack
```

A forward step on `j ∈ J⁻ᵢ(y)` increases `sᵢ` by `−ãᵢⱼ > 0` and only those
steps can ever cure row `i`. The maximum future increase achievable in row `i`
from `y` is bounded by:

```
Δᵢ_max(y) = Σ_{j ∈ J⁻ᵢ(y)}  (−ãᵢⱼ) · (dⱼ − yⱼ)         [recon, the
                                                     dissertation states the
                                                     argument verbatim]
```

**КН fires** at `y` iff there exists `i ∈ I⁻(y)` with

```
sᵢ(y) + Δᵢ_max(y) < 0          (2.9)
```

i.e. even taking all helping steps to their bounds cannot make row `i`
feasible. When КН fires, prune: take a backward step.

### 2.5 Plan-aware alternative-elimination criterion **КПИА** (eqs. 2.10–2.12)

Switches on as soon as the first feasible solution `y*_best` (with objective
`f*_best = c̃ · y*_best`) is known. Adds a **filter constraint** to the working
problem:

```
c̃ · y < f*_best            (2.10)  strict — used inside the lattice search
```

In slack form (with a dedicated filter slack `s_f`):

```
Σⱼ c̃ⱼ yⱼ + s_f = f*_best,   s_f ≥ 0       (2.11)  [recon, dissertation text
                                                  says “фильтрующее ограничение
                                                  преобразуется к (2.11)”]
```

**КПИА fires** on a candidate forward step on `yⱼ` from `y` iff

```
s_f(y) − c̃ⱼ ≤ 0           (2.12)
```

i.e. taking that step would non-strictly violate the filter, so no descendant
can beat `f*_best`. Prune that step (do **not** descend; try the next `j`).

Whenever a strictly better feasible `y` is found, replace `y*_best := y`,
`f*_best := c̃ · y`. The filter tightens monotonically.

### 2.6 Preferred-variable criterion **КПП** (eq. 2.13) — heuristic

Applied at every freshly entered infeasible point that survived КН: it
*reorders* `J(y)` (does not prune). Sort candidate forward steps `j ∈ J(y)`
**ascending** by

```
score(j) = Σ_{i ∈ I⁻(y) : ãᵢⱼ < 0}  ãᵢⱼ          (2.13)   [recon]
```

Intuition: `−ãᵢⱼ` is the gain in row `i`'s slack from one step on `j`; summing
the (negative) `ãᵢⱼ` over currently-violated rows means a *more-negative* score
== a *bigger* total help to violated constraints, so it should be tried first.
The dissertation explicitly calls КПП **heuristic**: it does not guarantee
faster discovery of a feasible, only that feasibility tends to be reached
earlier and КПИА kicks in earlier.

### 2.7 Two extra capabilities of the lattice enumerator

These are used directly by the Combined Method:

**(a) Search restricted to a partial hyperparallelepiped.** The enumerator can
run inside any axis-aligned box

```
P = { y : lⱼ ≤ yⱼ ≤ dⱼ }       partial hyperparallelepiped
```

with `l ≥ 0` instead of the full `[0, d]`. Because most search time is spent in
infeasible points, shrinking from `[0, d]` to `[l, d]` reduces the number of
visited points by far more than the volume reduction would suggest.

**(b) Variable fixing for slice-by-slice expansion.** The enumerator can fix
`yⱼ = vⱼ` and search only the slice. Two slices on different fixed values of the
same variable do not overlap, so the algorithm can grow the search region
incrementally without redoing work. Concretely, having searched `P =
[l, d]`, choosing some `j` with `lⱼ > 0`, fixing

```
vⱼ = lⱼ − 1                   (2.14)
```

and searching the slice

```
{ y : yⱼ = vⱼ;  lₖ ≤ yₖ ≤ dₖ for k ≠ j }     (2.15)
```

extends the cumulative search to

```
P' = { y : l'ⱼ ≤ yⱼ ≤ dⱼ;  lₖ ≤ yₖ ≤ dₖ for k ≠ j }   with l'ⱼ = vⱼ   (2.16)
```

This is the mechanism Stage 2 uses to grow the search box outward when the
initial box is empty of feasible points.

---

## 3. Combined Method — four stages

Stages, mirroring §2.2.4, §2.2.5, §2.3.

### Stage 1 — LP relaxation

1. Form (2.19).
2. Solve with the **special bounded-variable simplex** (§3.3.1) — a primal
   simplex variant that carries the `0 ≤ x' ≤ d̂` bounds *implicitly* in the
   tableau (variables can be at lower or upper bound; the basis algebra is the
   same as classical simplex). Keep the final tableau, basis `B` and the
   "variables-at-upper-bound" set `H` — the dissertation calls them out by name
   ("`H` is an attribute of the special bounded simplex used here").
3. Branch on the LP outcome:
   * **`Ω' = ∅`** (LP infeasible) ⇒ `Ω = ∅` (so (2.17) is infeasible too) ⇒
     `q* = ∅`, `f(q*) = ∞`, jump to Stage 5 (output).
   * **`q' is integer`** ⇒ `q* := q'`, jump to Stage 5.
   * **`q' fractional`** ⇒ compute `q*'` from `q'` by component-wise floor
     (eq. 2.20):

     ```
     q*'ⱼ = ⌊q'ⱼ⌋  for all j         (2.20)
     ```

     `q*'` is in general infeasible for both (2.17) and (2.19); it serves as
     the upper-right corner of the initial Stage-2 search box.

### Stage 2 — find a sub-optimal feasible

The Stage-2 box is

```
P₀ = { y : 0 ≤ yⱼ ≤ q*'ⱼ }              (2.21)
```

and the Stage-2 problem is

```
min c̃·y  s.t.  y ∈ P₀ ∩ Ω             (2.22)
```

solved by the lattice enumerator of §2 with `d ← q*'`. Because `q*'` sits next
to the LP optimum it is geometrically close to `Ω`'s boundary, so the
enumeration in (2.22) is cheap.

#### 2.1.1 Initial enumeration

Run the lattice enumerator on (2.22). If it finds a feasible — call it
`y_init` — go to 2.2. Else go to 2.1.2.

#### 2.1.2 Box expansion (eqs. 2.23, 2.24, 2.25)

When `P₀ ∩ Ω = ∅` (or whatever current `P` ∩ `Ω` is empty), grow `P` by slicing
along one variable at a time using §2.7(b). Choose the variable to expand by:

> **Expansion order rule.** Expand first along the variable with the smallest
> objective coefficient `c̃ⱼ`, i.e. **the smallest index** (since variables are
> already sorted by `c̃` ascending in §1.2 step 3).

For the chosen `j`:

```
vⱼ_new = lⱼ_current − 1                   (2.24)
solve   min c̃·y on slice { yⱼ = vⱼ_new, l_other ≤ y_other ≤ d_other }   (2.23)
```

Update `lⱼ ← vⱼ_new` after the slice search completes (using §2.7(b)
incrementally, no point is revisited). Skip variables for which:

* `c̃ⱼ = 0` is invalid for expansion in this direction (text on p.2.2.4 — those
  are non-improving anyway and the dissertation lists them as "невозможно").
* the КН criterion fired at the very start of the slice search — that variable
  is permanently excluded from further expansion (no slice will ever intersect
  `Ω`).

If every legal expansion has been tried and no feasible found ⇒
`Ω = ∅`, set `q* := ∅`, jump to Stage 5. (This is the (2.25) infeasibility
exit — rare, mainly for problems with equality constraints.)

#### 2.2 Cheap improvement of feasible (eq. 2.26)

Given `y_init`, look for "trivial improvements": variables `j` such that
`c̃ⱼ > 0` *and* `yⱼ` can be decreased without leaving `Ω`. From the candidate
set pick the one giving the largest objective decrease (which is
`c̃ⱼ · Δyⱼ_max`), apply the decrement, repeat until no candidate remains. Call
the result `y_sub`. If no improvement was possible, `y_sub := y_init` (2.26).

#### 2.3 Optional early exit by relative tolerance

If the user supplied a relative tolerance `ε` (per cent) on the objective, a
sound upper bound for the relative gap of `y_sub` is

```
gap_upper = (c̃·y_sub − (−obj_value(2.19))) / |c̃·y_sub|
          = (c̃·y_sub − f_LP) / |c̃·y_sub|        (2.27)  [recon — the
                                                        dissertation states
                                                        "оценка сверху
                                                        фактической погрешности"
                                                        based on the LP value]
```

If `gap_upper ≤ ε` (2.28) then accept `y_sub` and jump to Stage 5 with
`q* := y_sub`.

### Stage 3 — filter, then locate the start corner of the optimal-search box

Goal: build a small box that **provably** contains every integer feasible
strictly better than `y_sub` (if any), and start the corner of that box as
close as possible to `Ω`'s boundary so Stage 4 is cheap.

#### 3.1 Add the filter row to the LP (2.30) — *non-strict* form

In the LP setting we use the non-strict version of (2.10) because simplex
cannot handle strict inequalities:

```
c̃' · x'  ≤  f_sub                                 (2.30)

where  f_sub = c̃·y_sub  (objective value of the sub-optimal solution)
```

(Stage 4's lattice search will use the **strict** form (2.29) instead.)

The implementation injects this row directly into the saved final tableau of
(2.19) (which is why we kept it). The mechanical recipe is in §2.3.1 of the
dissertation, summarised:

* `m_LP ← m_LP + 1`, `n_LP ← n_LP + 1` — one new row, one new slack column.
* `b̂` extended with one new RHS `f_sub` (the slack `s_f` enters the basis).
* `Â` extended with a new row whose `j`-th entry is `c̃ⱼ` for `j ≤ n` and `0`
  elsewhere except a `1` in the new slack column. **[recon — the prose
  describes a row of `c̃ⱼ`s and a `1` in the diagonal slack position; the
  formulae layout in §2.3.1 was lost in the OLE export.]**
* `d̂` extended with `+∞` (or a large sentinel) for the new slack.
* New basis `B' = B ∪ {n_LP+1}`. The inverse-basis matrix `B'⁻¹` is grown by
  one row and one column (the explicit update is given in §2.3.1; bordering
  formula). The basic-solution vector gets one new entry equal to
  `f_sub − c̃·q'` (the residual filter slack at the LP optimum).

The filter row does **not** cut `q'` (`q'` is feasible for it because
`f_LP ≤ f_sub`).

#### 3.2 Compute the start corner `q''` of the Stage-4 box (eqs. 2.32–2.35)

We want, for each variable `yⱼ`, the smallest value `yⱼ` can take inside the
filtered region

```
Ω' = { x' ∈ Ω̂ : c̃'·x' ≤ f_sub }                (2.31)
```

so that the starting corner is

```
q''ⱼ = ⌈ min_{x' ∈ Ω'} x'ⱼ ⌉                   (2.32)
                                                (2.33)
```

with the rounding (2.34): `⌈z⌉ := smallest integer ≥ z`. Stage 4's search box
runs from `q''` up to `d̃`.

##### Which variables actually need an LP solve

Per (2.35) — and this is the key efficiency claim of Stage 3 — we only need
to solve (2.33) for variables in

```
S = B(q') ∪ H(q')                              (2.35)
```

* `B(q')` = basic variables of the LP optimum.
* `H(q')` = non-basic variables that sat at their **upper** bound in `q'`
  (the bounded-simplex's "upper-bound set"). Sitting at upper means
  `q'ⱼ = dⱼ`, not `0`; the LP is at an extreme point of `Ω̂`, so for any
  variable that is non-basic at zero, zero already *is* the minimum and there
  is nothing to compute — set `q''ⱼ = 0`.

For each `j ∈ S`:

##### LP variant (2.33)

Solve

```
min  x'ⱼ                                       (2.33)
s.t. x' ∈ Ω̂  AND  c̃'·x' ≤ f_sub
```

starting **from the saved tableau** of the filtered (2.19). Because the
objective is "isolate one variable" and the filter row does not cut `q'`,
typical iteration count per variable is 0 (when `q'ⱼ` is already minimal) or 1
(when the minimum is across one edge from `q'`). Procedural recipe in §2.3.2
(also lost from the OLE export — see the gap log at the end). Then

```
q''ⱼ = ⌈ x'ⱼ_min ⌉                              (2.34)
```

For `j ∉ S`: `q''ⱼ = 0`.

### Stage 4 — provably-optimal lattice search

Run the lattice enumerator (§2) on

```
min  c̃·y
s.t. y ∈ Ω
     q''ⱼ ≤ yⱼ ≤ dⱼ                            (2.36)
     y ∈ Z^n_+
with the filter row in **strict** form
     c̃·y < f_sub                              (2.29)
active from the very first step (so КПИА fires from step 0 with
y*_best := y_sub, f*_best := f_sub).
```

Two short-circuit cases:

* **Empty box.** If `q'' ≥ d̃` componentwise breaks down, or if the strict
  filter and the box together imply emptiness on inspection (e.g. the Stage-4
  region is provably a subset of the Stage-2 region already exhausted), declare
  `q* := y_sub` (the sub-optimal *was* optimal) and skip the lattice run.
* **Lattice search yields a feasible better than `y_sub`.** Then `q* :=`
  that result.
* **Otherwise** (lattice search exhausts without improvement) `q* := y_sub`.

### Stage 5 — output

Map `q*` back to the original variable set by inverting (2.18) (`xⱼ = dⱼ − yⱼ`
on whichever indices were substituted) and the re-numbering of step 1.2(3).
Apply the same inverse if the algorithm exited from Stage 2.3 on tolerance.

---

## 4. Heuristic plug-ins (§2.4)

> "Применение эвристических алгоритмов … не влияет на возможность нахождения
> строго оптимального решения, а лишь обеспечивает более быстрый его поиск."

The architectural commitment is that heuristics are confined to **Stage 2**:
they may *replace* the lattice enumerator's job of producing the first feasible
`y_init`, but they may not influence Stages 1, 3 or 4. Concretely:

* The heuristic is exposed as a module that takes the input problem (in form
  (2.17)) plus the LP solution `q'` from Stage 1 and returns either an integer
  feasible point or "no result".
* If it returns a point, that point is treated exactly as if the lattice
  enumerator had returned it: cheap improvement (§3 stage 2.2) and tolerance
  check (§2.3) are applied.
* If it returns "no result", Stage 2 falls back to the lattice enumerator.

For the LAN-topology application the dissertation supplies a domain-specific
heuristic (Chapter 4); it is not part of the generic method.

КПП (§2.6) is the *only* heuristic embedded inside the lattice enumerator
itself, and it only affects ordering, never pruning — so optimality survives
even with КПП disabled.

---

## 5. Properties (§2.5)

Claims the dissertation makes about the method, in implementation-relevant
terms:

* **Optimality guarantee.** If `Ω ≠ ∅` and the objective is bounded over the
  integers (granted by the upper bounds `d`), Stages 1–4 return a strict
  optimum. Stage 4's lattice search inside `[q'', d̃]` with the strict filter
  is the certificate: either it finds a strictly better solution, or the
  exhausted search proves `y_sub` was optimal.
* **Numerical robustness vs. cutting-plane / B&B.** The filter row is a
  *single* extra LP row added once; there is no cumulative cut family that
  would cause coefficient blow-up the way Gomory cuts do. Branching is
  replaced by lattice enumeration in a small box. The dissertation flags
  numerical safety as a major selling point against МПО (правильные
  отсечения = cutting planes) and МВГ (ветвей и границ = B&B).
* **Scaling.** Author's a-priori claim: gap to МПО/МВГ widens with `n`. The
  small-test results in §3.2 use eight crafted problems exercising all
  branches; Chapter 4 uses a real LAN-topology instance solved at РКК
  «Энергия».
* **Heuristic openness.** A user can plug in a Stage-2 heuristic without
  risking optimality.

---

## 6. Software architecture (§3.1)

For reference; useful mainly for choosing module boundaries.

* Language: **C++** (Borland C++ 5.5 in the dissertation, also tested with
  GCC under UNIX). Reimplementation in any modern language is fine.
* Module map (§3.1.3, "Схема межмодульных связей"):
  * `IN`     — read problem in the dissertation's input format and echo-print.
  * `OUT`    — print result (`q*` and `c · q*`).
  * `MAIN`   — orchestrates the four stages.
  * `TLP`    — bounded-variable simplex *and* the dual bounded-variable
    simplex (the latter is needed both inside the regular B&B comparator and
    by Stage 3's per-variable LPs reusing the saved tableau).
  * `CONV`   — transformations to/from canonical forms (1.2, 1.3 above).
  * `SVL`    — lattice enumerator (§2).
  * Heuristic modules slot under Stage 2 with the contract from §4.
* The control flow inside `MAIN` is linear (Stage 1 → 2 → 3 → 4 → 5) but may
  short-circuit (LP infeasible, integer LP optimum, tolerance hit, Stage 4
  empty).
* The bounded-variable simplex code is shared across the comparison
  experiments (B&B, cutting-plane, combined) — same data structure, three
  drivers.

---

## 7. Reference pseudocode

```
function combined_method(c, A, b, sense, d, eps=None):
    # ----- 1.2  Build canonical (2.17) -----
    if sense == "max":
        c = -c
    flip = []                       # vars where x = d - y substitution happened
    for j where c[j] < 0:
        c[j], A[:,j], b = -c[j], -A[:,j], b - A[:,j]*d[j]
        flip.append(j)
    perm = argsort(c)               # ascending by c
    apply perm to (c, A, d)
    Atil, btil = split_equalities_and_flip_geq(A, b)   # only used by SVL

    # ----- 1.3  Build LP form (2.19) -----
    Ahat, bhat, dhat, chat = build_LP_canonical(c, A, b, sense, d)

    # ----- Stage 1 -----
    lp_status, q', tableau, B, H = bounded_simplex(Ahat, bhat, dhat, chat, maximise=True)
    if lp_status == INFEASIBLE:                      return None
    if is_integer(q'):                               return decode(q', perm, flip, d)

    # ----- Stage 2 -----
    qstar' = floor(q')                    # (2.20)
    box   = HyperBox(low=0, high=qstar')  # (2.21)

    y_init = svl_search(Atil, btil, c, box, filter=None)
    while y_init is None:
        j = pick_smallest_c_index_with_low_gt_0(box)
        if j is None:                       return None    # Ω = ∅  (2.25)
        if KH_at_start_of_slice(box, j):    box.exclude_var(j); continue
        v = box.low[j] - 1                  # (2.24)
        y_slice = svl_search_slice(Atil, btil, c, box, j, v)
        if y_slice is not None:             y_init = y_slice
        box.low[j] = v                      # (2.16)

    y_sub = trivial_improve(y_init, c, Atil, btil, d)   # (2.26)

    if eps is not None:
        f_LP  = -tableau.objective          # back to min sign
        f_sub = c · y_sub
        if (f_sub - f_LP) / abs(f_sub) <= eps:
            return decode(y_sub, perm, flip, d)

    # ----- Stage 3 -----
    f_sub = c · y_sub
    add_filter_row_to_tableau(tableau, c, f_sub)        # §2.3.1

    qpp = zeros(n)                                       # q''
    S   = basis(B) ∪ at_upper_bound(H)                   # (2.35)
    for j in S:
        z_min = solve_min_xj(tableau, j)                 # (2.33), warm-started
        qpp[j] = ceil(z_min)                             # (2.34)

    # ----- Stage 4 -----
    box4 = HyperBox(low=qpp, high=d)
    if box4.is_empty() or box4 ⊆ already_searched_in_stage_2:
        q_star = y_sub
    else:
        better = svl_search(Atil, btil, c, box4,
                            filter=StrictFilter(c, f_sub, y_sub))   # (2.36) + (2.29)
        q_star = better if better is not None else y_sub

    # ----- Stage 5 -----
    return decode(q_star, perm, flip, d)
```

```
function svl_search(A, b, c, box, filter):
    # Pre: c[j] ≥ 0, sorted asc, integer.
    s_init  = b - A @ box.low                            # slack at the corner
    j_last  = -1
    y       = box.low.copy()
    s       = s_init
    best    = filter.best if filter else None
    f_best  = filter.f_best if filter else +inf
    stack   = []           # for backtracking

    enter(y, j_last):
        if all(s ≥ 0):                                   # feasible
            if c·y < f_best:
                best, f_best = y.copy(), c·y
                filter.update(best, f_best)
            backtrack()
            return
        if KH_fires(s, A, box, y, j_last):              # (2.7)–(2.9)
            backtrack()
            return
        cands = [ j ≥ j_last : y[j] < box.high[j]
                              and (filter is None or s_filter - c[j] > 0) ]   # (2.12)
        order cands ascending by score = Σ_{i ∈ I⁻(y) : A[i,j] < 0} A[i,j]    # (2.13)
        for j in cands:
            push(j_last); y[j]++; s -= A[:,j]; enter(y, j); y[j]--; s += A[:,j]
        backtrack()

    enter(y, -1)
    return best
```

---

## 8. Parameters and tunables

* `eps` — relative tolerance for early exit on Stage 2.3 (default: not set).
* `expand_order` — overridable rule for which variable Stage 2.1.2 expands on
  next; default is "smallest `c̃ⱼ` index" (which after sorting is just "smallest
  index").
* `heuristic` — pluggable Stage-2 producer of `y_init`.
* `enable_KPP` — toggle the КПП ordering inside the lattice enumerator (does
  not change correctness, only running time).
* `simplex_pivot_rule` — affects warm-starts in Stage 3. Bland's rule keeps
  the per-variable LP solves bounded; default is what the bounded-variable
  simplex used in `TLP`.

---

## 9. Gaps in the source / things to verify against the original

These are the bits where the OLE-equation export left the formula blank and
where the prose alone is not 100% mechanical. Anything you implement that
depends on these should be tested against a hand example before trusting it.

1. **(2.9) КН exact form.** The prose justifies it but the `Δᵢ_max` summation
   is reconstructed. It is the natural one and matches every standard
   reference for vector-lattice enumeration; verify on §3.2 example 1.
2. **(2.13) КПП score.** Same situation — reconstructed sign convention;
   ascending order should pick "most-helpful" steps first.
3. **§2.3.1 tableau bordering formulae.** The dissertation gives explicit
   bordering formulae for the inverse basis and the basic-solution vector
   when the filter row is added; the formula entries did not survive export
   (rows of `EMBED Equation.3`). The prose is enough to write correct code,
   but cross-check with example 5 (the "all stages exercised" test).
4. **§2.3.2 warm-start formulae** for solving (2.33) per variable. Same
   situation as 3.
5. **(2.27) gap_upper.** The dissertation says "upper bound on the actual
   error", which for a min problem with LP lower bound `f_LP` is
   `(f_sub − f_LP) / |f_sub|`. Verify wording before shipping a tolerance
   gate.
6. **Filter slack `+∞` upper bound.** Implementation detail: the new slack
   column for the filter row needs a sentinel upper bound; choose something
   like `1e18` (or model-aware) and assert it never binds.
7. **Variable expansion exclusion in Stage 2.1.2.** The dissertation says
   "exclude variables where КН fired at the start of the slice". Make sure
   "start of the slice" means the freshly-fixed corner, not the parent box
   corner.

If any of (1)–(4) need to be exact for compliance, the cleanest fix is to
re-render the `.doc` to PDF on a machine with Word/LibreOffice (the install
on this machine did not complete) and read those formula images directly.
The 898 `Equation Native` OLE streams are intact in `source.doc`, so the
information is recoverable; it just was not rendered by `textutil`/`pandoc`.
