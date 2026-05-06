# Combined Method — implementation spec

Source: `Диссертация.doc` by А.Ю. Дёмин (МАИ, 2002), supervisor Г.Ф. Хахулин.
Title: «Комбинированный метод решения линейных задач целочисленного программирования
и его применение для оптимизации топологии локальных вычислительных сетей».

This is a working translation/condensation of Chapter 2 (and the relevant
software-architecture bits of Chapter 3) sufficient to re-implement the method
in code. Section numbers below match the dissertation. Notation follows the
dissertation verbatim where possible — Russian/Cyrillic terms in parentheses
are load-bearing names you'll see in the original.

> **Provenance.** The original `.doc` stores every equation as an embedded
> Microsoft Equation 3.0 OLE object (898 streams). After installing
> LibreOffice and converting `.doc → .docx → .md` via pandoc, all formulas
> render correctly. The text below uses the dissertation's exact symbols.
> Equation numbers cite the dissertation directly. Where the rendered
> markdown for a formula was layout-broken (a few cases in §2.3.1 / §2.3.2),
> the formula has been reconstructed from prose — flagged **[recon]**.

---

## 0. Notation cheat-sheet

| Symbol                   | Meaning                                                                                                               |
| ------------------------ | --------------------------------------------------------------------------------------------------------------------- |
| `x`                      | `n`-vector of decision variables (integer in (2.1), (2.17); continuous in (2.19)). Entries `x[j]`                     |
| `c`                      | row-vector of objective coefficients, length `n`                                                                      |
| `A`                      | constraint matrix, `m × n`                                                                                            |
| `b`                      | RHS vector, length `m`                                                                                                |
| `h`                      | per-variable upper-bound vector, length `n`, integer, finite                                                          |
| `D`                      | feasible set of (2.1) (the user form)                                                                                 |
| `^p` superscript         | "problem-form" — variant of objects after the canonical rewrite (2.17). Lattice search consumes these                 |
| `^l` superscript         | "LP-form" — variant for the bounded-variable simplex (2.19)                                                           |
| `x^l = (x_0^l, x_d^l)^T` | LP variables: `x_0^l` are the structural vars (matching `x^p`); `x_d^l` are slacks introduced to make rows equalities |
| `D^p`, `D^l`             | feasible sets of (2.17) and (2.19)                                                                                    |
| `x_opt^l`, `Z_opt^l`     | LP optimum point and objective                                                                                        |
| `x_opt^p`, `Z_opt^p`     | integer optimum of (2.17)                                                                                             |
| `x_min^p`                | floor of `x_opt^l` — Stage-2 search-box upper corner (eq. 2.20)                                                       |
| `x_min`                  | start corner of the Stage-4 search box (eq. 2.34)                                                                     |
| `x_D^p`, `z_D^p`         | first feasible found in Stage 2, and its objective                                                                    |
| `x̃_D^p`, `z̃_D^p`       | improved feasible after Stage-2 sub-step 2.2                                                                          |
| `B_opt^l`                | optimal basis of (2.19)                                                                                               |
| `N(B_opt^l)`             | index set of basic variables at LP optimum                                                                            |
| `N_opt^l`                | index set of variables sitting at their **upper** bound at LP optimum (specific to bounded-variable simplex)          |
| `J_x`                    | indices of forward steps available from lattice point `x` (eq. 2.2)                                                   |
| `j_x`                    | "level" of `x`: the largest index `j` with `x[j] ≠ 0` (eq. 2.3)                                                       |
| `y[i]`                   | slack of the `i`-th main constraint, used inside the lattice search (eq. 2.6)                                         |
| `y[0]`                   | slack of the filter constraint (eq. 2.11)                                                                             |
| `I_x`                    | indices of currently-violated rows at `x` (eq. 2.7)                                                                   |
| `J^x_i`                  | forward steps from `x` whose coefficient in row `i` is negative (eq. 2.8) — only these can repair row `i`             |
| `[a]`                    | floor (целая часть) — used in (2.20)                                                                                  |
| `]a[`                    | smallest integer ≥ a (ceiling) — used in (2.34)                                                                       |
| КН, КПИА, КПП            | the three lattice-search criteria (infeasibility / plan-aware / preferred-variable)                                   |

---

## 1. Problem statement (§2.2.1, eq. 2.1)

```
min   Z = c · x
s.t.  D = { x ∈ N^n : A x {≤,=,≥} b              (2.1.1)
                    0 ≤ x ≤ h                    (2.1.2) }
```

* `x ∈ N^n` — non-negative integers, `n`-vector.
* `c ∈ R^n` — **arbitrary real** in user form. The dissertation chooses to write
  `min` — that direction is forced internally by the canonical rewrite below.
* `A`, `b` — main constraints; rows are mixed `≤`, `=`, `≥`.
* `h ∈ N^n` — every variable has a finite integer upper bound. In LAN-topology
  problems it falls out of the formulation; in the general case it is supplied
  as an a-priori envelope.

The application that drives Chapter 4 is the LAN-topology problem (§1.1):
selecting a hub hierarchy with hard limits on ports per hub, total cabling
length, etc. For the generic library implementation that's irrelevant — we
just solve (2.1).

---

## 2. Vector-lattice implicit enumeration (§2.2.2)

This is the integer-search primitive used in Stages 2 and 4.

### 2.1 Lattice traversal

A *route* (маршрут) is a path through the integer lattice `[0..h] ⊂ N^n` made
of **forward** and **backward** steps. A forward step on variable `j`
("step forward by `j`"): increment `x[j]` by 1, descend one level. A backward
step: decrement `x[j]`, return one level. The traversal starts and finishes at
the origin `0`. It finishes when, at the origin, every legal forward step has
been explored.

To prevent revisiting a point along different orderings, at every node `x`
the algorithm restricts forward steps to a subset `J_x` of variable indices:

```
j_x = max { j ∈ 1..n : x[j] ≠ 0 }              (2.3)

J_x = { j ∈ 1..n : j ≥ j_x }                    (2.2)
       — and j = j_x is included in J_x only if x[j_x] < h[j_x]
```

A point with `x[n] = h[n]` is **terminal** (конечная); after entering a
terminal point on a forward step a backward step must follow.

### 2.2 First reduction rule — non-negative `c` (eq. 2.4)

Required for the lattice search to terminate. For any `j` with `c[j] < 0`,
substitute

```
x[j]  =  h[j] − x'[j]                          (2.4)
```

throughout the objective and **all** constraints. This makes
`c'[j] = −c[j] > 0`. The constant term that appears in the objective is
dropped (does not affect the argmin).

After step (2.4) on every offending variable, every `c[j] ≥ 0`. Forward steps
strictly worsen the objective, so:

> **Rule 1.** Once a forward step lands in a feasible point, do not take any
> further forward step from it — back up.

Corollary: the algorithm spends almost all its time in *infeasible* points;
feasible leaves end branches. If `0` itself is feasible, it is the global
optimum and the search is over.

### 2.3 Slacks for constraint tracking (eqs. 2.5–2.6)

After the (2.17) rewrite (§3) all main constraints are `≤`. Introduce a slack
`y[i]` per row, requiring `y[i] ≥ 0`:

```
Σⱼ a[i,j] x[j] + y[i] = b[i],   y[i] ≥ 0,   i = 1..m           (2.6)
y_x[i] := b[i] − Σⱼ a[i,j] x[j]
```

A point `x` is feasible iff `y_x[i] ≥ 0 ∀ i`. Each forward step on `x[j]`
updates every slack by `y[i] ← y[i] − a[i,j]`.

### 2.4 Infeasibility criterion **КН** (eqs. 2.7–2.9)

Applied at every newly entered point.

```
I_x   = { i : y_x[i] < 0 }                                     (2.7)
J^x_i = { j ∈ J_x : a[i,j] < 0 }                                (2.8)
```

`I_x` is the set of currently-violated rows; `J^x_i` is the set of forward
steps that *help* row `i` (one step on `j ∈ J^x_i` increases `y[i]` by
`−a[i,j] > 0`).

**КН fires** at `x` iff there exists `i ∈ I_x` with

```
Σ_{j ∈ J^x_i}  a[i,j] · (h[j] − x[j])   >   y_x[i]              (2.9)
```

Both sides are negative. The LHS is the most negative thing that the
remaining helping-steps can add to `y[i]` (note `a[i,j] < 0`,
`h[j] − x[j] ≥ 0`). The condition says: even taking *all* helping steps to
their bounds cannot lift `y[i]` to non-negative — row `i` cannot be cured. So
prune: take a backward step.

> Equivalent re-statement, easier to read in code:
> ```
> max_increase_in_y[i] = Σ_{j ∈ J^x_i}  (−a[i,j]) · (h[j] − x[j])
> KH fires  ⇔  ∃ i ∈ I_x : y_x[i] + max_increase_in_y[i] < 0
> ```

### 2.5 Plan-aware alternative-elimination criterion **КПИА** (eqs. 2.10–2.12)

Switches on as soon as the first feasible solution `x_∂` (with objective
`z(x_∂) = c · x_∂`) is known. Adds a **filter constraint** to the working
problem:

```
Σⱼ c[j] x[j]   <   z(x_∂)                                       (2.10)
```

In slack form (with a dedicated filter slack `y[0]`):

```
Σⱼ c[j] x[j] + y[0]  =  z(x_∂),   y[0] > 0                      (2.11)
y_x[0] = z(x_∂) − Σⱼ c[j] x[j]
```

Each forward step on `x[j]` updates the filter slack by `y[0] ← y[0] − c[j]`.

**КПИА fires** on a candidate forward step on `x[j]` from `x` iff

```
c[j]  ≥  y_x[0]                                                 (2.12)
```

i.e. taking that step would non-strictly violate the filter, so no descendant
can beat the incumbent. Prune that step (do **not** descend; try the next `j`).

Whenever a strictly better feasible `x` is found, replace `x_∂ := x`,
`z(x_∂) := c · x`. The filter tightens monotonically.

### 2.6 Preferred-variable criterion **КПП** (eq. 2.13) — heuristic

Applied at every freshly entered infeasible point that survived КН: it
*reorders* `J_x` (does not prune). Sort candidate forward steps `j ∈ J_x`
**ascending** by

```
score(j)  =  Σ_{i ∈ I_x}  a[i,j] · (h[j] − x[j])                 (2.13)
```

The sum runs over **all** currently-violated rows (not only those where
`a[i,j] < 0`). For `i ∈ I_x` with `a[i,j] < 0`, the term is non-positive (a
forward step on `j` helps row `i`); with `a[i,j] > 0`, the term is
non-negative (a forward step hurts row `i`). Ascending order picks the
biggest *net help* first.

The dissertation calls КПП **heuristic**: it does not guarantee faster
discovery of feasibility, only that feasibility tends to be reached earlier
and so КПИА kicks in earlier. Disabling it does not affect correctness.

### 2.7 Two extra capabilities of the lattice enumerator

Both used directly by the Combined Method.

**(a) Search restricted to a partial hyperparallelepiped.** The enumerator can
run inside any axis-aligned box

```
P = { x : x_min[j] ≤ x[j] ≤ h[j],  j = 1..n }
```

with `x_min ≥ 0` instead of the full `[0, h]`. Because most search time is
spent in infeasible points, shrinking from `[0, h]` to `[x_min, h]` reduces
the visited count by far more than the volume reduction would suggest.

**(b) Variable fixing for slice-by-slice expansion.** The enumerator can fix
`x[k] = v` and search only the slice. Two slices on different fixed values of
the same variable do not overlap, so the algorithm grows the search region
incrementally without redoing work. Having searched `P = {x_min[j], h[j]}`
for all `j`, choosing some `k` with `x_min[k] > 0` and setting

```
x_min[k] := x_min[k] − 1                                         (2.14)
```

then searching the slice

```
{ x : x[k] = x_min[k];   x_min[j] ≤ x[j] ≤ h[j],  j ≠ k }        (2.15)
```

extends the cumulative search to the new (larger) box

```
{ x_min[j] ≤ x[j] ≤ h[j],  j = 1..n }                            (2.16)
```

with the updated `x_min[k]`. This is the mechanism Stage 2 uses to grow the
search box outward when the initial box is empty of feasibles.

---

## 3. Canonical forms (§2.2.3)

### 3.1 (2.17) — for the lattice enumerator

The Combined Method requires the problem in the form

```
min   Z^p  =  c^p · x^p
s.t.  D^p = { x^p ∈ N^n :  A^p x^p ≤ b^p,   0 ≤ x^p ≤ h^p }      (2.17)
      c^p[j] ≥ 0  ∀ j
      c^p[1] ≤ c^p[2] ≤ … ≤ c^p[n]
```

Reaching (2.17) from (2.1):

1. **Force minimisation.** If the original is `max c·x`, replace `c ← −c`.
2. **Force non-negative cost coefficients.** For every `j` with `c[j] < 0`,
   apply (2.4): `x[j] := h[j] − x^p[j]`. Substitute everywhere; drop the
   constant term that appears in the objective. After this, all `c^p[j] ≥ 0`.
3. **Re-number variables ascending by `c^p[j]`** so `c^p[1] ≤ … ≤ c^p[n]`.
   `h^p` is `h` permuted accordingly. Required both by the lattice enumerator
   (it relies on the level/index correspondence) and by Stage 2.1.2 (which
   expands variables in this order).
4. **Force `≤` constraints.** Equality rows split into `≤` and `≥`; each `≥`
   row multiplied by `−1`. Result: matrix `A^p ∈ R^{m^p × n}` with
   `m^p ≥ m`, RHS `b^p ∈ R^{m^p}`.

The optimum of (2.17), if it exists, is `(Z_opt^p, x_opt^p)`.

### 3.2 (2.19) — for the bounded-variable simplex

The LP companion is built **without** the `≥ → ≤` flip and **without**
splitting equalities; it uses only the substitution (2.18).

```
max   Z^l  =  c^l · x^l
s.t.  D^l = { x^l ∈ R^{n^l} :  A^l x^l = b^l,   0 ≤ x^l ≤ h^l }  (2.19)
```

with

```
x^l = (x_0^l, x_d^l)^T,   length n^l = n + m
c^l = (−c^p, 0, …, 0)
h^l = (h^p, ∞, …, ∞)^T
```

* `x_0^l` are the structural vars, numbered identically to `x^p` (same vector
  after substitution and re-numbering).
* `x_d^l` are slacks introduced to turn each main row into an equality.
* `c^l = (−c^p, 0…)` because canonical-form simplex maximises. Hence
  `Z^p = −Z^l` for the LP relaxation: the values differ only in sign.
* `A^l` is built from the *original* `A` after substitution (2.18), with
  slack columns appended. The flips and splittings done for `A^p` are **not**
  applied here.
* `b^l` is `b` after substitution (2.18).
* `h^l` extends `h^p` with `+∞` slack upper bounds.

The optimum of (2.19), if it exists, is `(Z_opt^l, x_opt^l)`.

> The dissertation insists on the **special bounded-variable simplex** (with
> algorithmic two-sided bounds) — variables can sit at lower bound, basic, or
> at upper bound; the basis algebra is classical but the LP solver tracks an
> extra index set `N_opt^l` of vars currently at their upper bound. This is
> needed by Stage 3 (the start-corner derivation references both `N(B_opt^l)`
> and `N_opt^l` — see (2.35)).

---

## 4. The four stages

Stages mirror §2.2.4–§2.2.5 and the algorithm in §2.3 verbatim.

### Stage 1 — LP relaxation (algorithm step 2)

1. Form (2.19).
2. Solve with the special bounded-variable simplex. Keep the **final tableau**
   (will be reused in Stage 3): basis `B_opt^l`, basis index set
   `N(B_opt^l)`, upper-bound index set `N_opt^l`, basic-solution vector
   `b(B_opt^l, N_opt^l)`, inverse basis `B^{-1}`.
3. Branch on the LP outcome:
   * **`D^l = ∅`** ⇒ `D^p = ∅` and `D = ∅`, jump to Stage 5 (output
     "infeasible").
   * **`x_opt^l ∈ D^p`** (i.e. integer feasible for (2.17)) ⇒
     `x_opt^p := x_opt^l`, jump to Stage 5.
   * **`x_opt^l ∉ D^p`** ⇒ compute `x_min^p` by component-wise floor:

     ```
     x_min^p[j]  =  [ x_opt^l[j] ]   for j = 1..n                 (2.20)
     ```

     `x_min^p` is in general infeasible for both (2.17) and (2.19); it
     defines the upper-right corner of the initial Stage-2 search box.

### Stage 2 — sub-optimal feasible (algorithm step 3)

The Stage-2 search box is

```
P = { x^p : x_min^p[j] ≤ x^p[j] ≤ h^p[j], j = 1..n }              (2.21)
```

> ⚠ **Convention note.** `x_min^p` is the *upper-right* corner of `P` here
> (the box runs from somewhere to `h^p` and `x_min^p` ≤ `h^p`), but the
> dissertation later "moves `x_min^p` toward the origin" in Stage 2.1.2.
> The cleanest reading: the box always runs from `[x_min^p, h^p]` (lower-left
> corner is `x_min^p`, upper-right is `h^p`), and "shifting `x_min^p` toward
> 0" *enlarges* the box. The (2.20) floor of the LP optimum gives the
> *largest* `x_min^p` that still touches the LP corner, and Stage 2.1.2
> shrinks `x_min^p` if needed.

The Stage-2 problem is

```
min  c^p · x^p   s.t.  x^p ∈ P ∩ D^p                              (2.22)
```

solved by the lattice enumerator of §2 with `box = [x_min^p, h^p]`. Because
`x_min^p` sits next to the LP optimum, the enumeration is cheap.

#### Step 3.1.1 — initial enumeration

Run the lattice enumerator on (2.22). If a feasible is found — call it
`x_D^p`, with objective `z_D^p` — go to step 3.2. Else go to 3.1.4.

#### Step 3.1.4 — choose expansion variable; 3.1.5 — solve sliced sub-problem

When `P ∩ D^p = ∅` (current box has no feasible), grow `P` by slicing along
one variable at a time (§2.7(b)). Variable choice rule:

> **Expansion order.** Expand first along the variable with the smallest
> objective coefficient `c^p[j]` — i.e. **the smallest index** (variables are
> already sorted ascending by `c^p` per §3.1).

Skip variables for which:

* `x_min^p[j] = 0` — cannot decrement.
* КН fired at the very start (initial corner) of the slice search on `j` —
  permanently exclude `j`, no slice on `j` will ever intersect `D^p`.

For the chosen `k`:

```
x_min^p[k] := x_min^p[k] − 1                                      (2.24)
```

then solve the sliced (2.23):

```
min c^p · x^p
s.t. x^p ∈ N^n, A^p x^p ≤ b^p, x[k] = x_min^p[k]_new,
     x_min^p[j] ≤ x[j] ≤ h^p[j] for j ≠ k                          (2.23)
```

with the lattice enumerator. (The variable-fix capability of §2.7(b) lets
this slice be searched without overlapping any previously-searched region.)
Update the running `x_min^p[k]` after the slice search; loop to 3.1.1 to
either find feasibility or pick the next expansion variable.

If every legal expansion has been tried and no feasible found:
`D^p = ∅` and `D = ∅` (note this implies (2.25): `D^l ≠ ∅ ∩ D^p = ∅`). Jump
to Stage 5 (rare; mainly for problems with equality constraints).

#### Step 3.2 — cheap improvement of feasible

Given `(z_D^p, x_D^p)`, look for "trivial improvements" — variables `j` such
that `c^p[j] ≥ 0` (every variable, after (2.17)) and `x_D^p[j]` can be
decreased without leaving `D^p`. Improvement size:

```
Δz[j]  =  c^p[j] · ( x_D^p[j] − x̃_D^p[j] )
```

where `x̃_D^p[j]` is the smallest feasible value `x[j]` can take with all
other coordinates held. Pick the `j` with biggest `Δz[j]`, apply, repeat
until no candidate remains. Call the result `(z̃_D^p, x̃_D^p)`. If no
improvement was possible:

```
( z̃_D^p, x̃_D^p )  ≡  ( z_D^p, x_D^p )                            (2.26)
```

#### Step 3.3 — optional early exit by relative tolerance

If the user supplies a relative tolerance `δ_Z` (per cent) on the objective:

```
δ̂_Z  =  ( z̃_D^p − Z_opt^l ) / z̃_D^p  ·  100%                    (2.27)
```

(`Z_opt^l` here is the LP-relaxation optimum *in (2.17) sign convention* —
i.e. `−`(value returned by (2.19)'s simplex). It is a lower bound on
`Z_opt^p`, so `δ̂_Z` upper-bounds the true relative gap.) If

```
δ̂_Z  ≤  δ_Z                                                       (2.28)
```

accept `x̃_D^p` and jump to Stage 5 with `x_opt^p := x̃_D^p`.

### Stage 3 — filter, then locate `x_min` (algorithm step 4)

Goal: build a small box that **provably** contains every integer feasible
strictly better than `x̃_D^p` (if any), with its lower-left corner `x_min`
as close as possible to `D^l`'s boundary so Stage 4 is cheap.

#### Step 4.1 — add filter row to the LP (2.30)

LP cannot handle strict inequalities, so:

```
c^p · x_0^l   ≤   z̃_D^p                                          (2.30)
```

(non-strict; Stage 4 will use the strict form (2.29) inside the lattice
search.) Inject this row into the saved final tableau of (2.19). Mechanical
recipe (§2.3.1):

```
m̃ = m + 1                  one new row
ñ = n^l + 1                 one new variable (filter slack y[0])

b̃^l   =  ( b^l, z̃_D^p )^T                          (RHS)
h̃^l   =  ( h^l, ∞ )^T                              (upper bounds)

Ã^l  =  ┌  A^l                          | 0 ┐
        │  c^p[1] … c^p[n]  | 0 … 0     | 1 │
        └                                       ┘
       (top: A^l with a 0 column for the new slack;
        bottom row: c^p in the n structural-var positions,
        zeros in the m existing-slack positions, 1 in the filter-slack column)

Ñ(B)  =  N(B_opt^l) ∪ {ñ}        (extended basis index set)

B̃^{-1}  =  ┌  B^{-1}      | 0 ┐
            │  d_m̃        | 1 │
            └                  ┘
        where  d_m̃[i] = − Σ_{i'=1..m}  c^p[ j_{i'} ] · (B^{-1})[i', i]
                       — equivalently, d_m̃ = − c_B(B_opt^l) · B^{-1}
                       — j_{i'} are the basis indices of B_opt^l.

basic-solution vector  b̃(B,N)
        =  ( b(B_opt^l, N_opt^l),   z̃_D^p − c^p · x_opt^l[1..n] )^T
                                  └────────── filter slack at the LP optimum
```

The filter row does **not** cut `x_opt^l` (because `c^p · x_opt^l[1..n] =
−Z_opt^l ≤ z̃_D^p`), so the augmented tableau is still primal-feasible at
`x_opt^l`. It can serve as a warm start for the per-variable LPs that follow.

#### Step 4.2 — compute the start corner `x_min` (eqs. 2.32–2.35)

For each variable `j` we want the smallest value `x[j]` can take inside the
filtered LP region

```
D_f^l = D^l ∩ { c^p · x_0^l ≤ z̃_D^p }                            (2.31)
```

so that the start corner is

```
x_min[j]  =  ] z_jopt^l [                                         (2.34)

where  z_jopt^l = min_{x^l ∈ D_f^l}  x^l[j]                       (2.32)
                = max_{x^l ∈ D_f^l}  −x^l[j]                      (2.33)
```

`]a[` = smallest integer ≥ a (ceiling, non-standard notation in Russian OR).

##### Which variables actually need an LP solve

```
J  =  N(B_opt^l) ∪ N_opt^l                                        (2.35)
```

* `j ∉ J` ⇒ `x_opt^l[j] = 0` and is non-basic at zero. Zero already is its
  minimum across `D_f^l` (the filter row does not cut `x_opt^l`, the LP is at
  an extreme point, no edge to follow). Set `x_min[j] := 0` directly.
* `j ∈ J` ⇒ solve (2.33) for `x[j]`.

##### Per-variable LP recipe (§2.3.2)

For each `j ∈ J`, solve a fresh LP (2.33) **warm-started from the augmented
tableau of step 4.1**:

```
N(B_1)   = Ñ(B)                          initial basis = augmented basis
N_1      = N_opt^l                       initial upper-bound set
m_1      = m̃                             rows
n_1      = ñ                             cols
A̅        = Ã^l                            constraint matrix (filter included)
b̅        = b̃^l                            RHS
h̅        = h̃^l                            upper bounds
c̅        = (0, …, 0, −1 at position j, 0, …, 0)         ← objective: −x[j]
                                                          (so simplex maximises
                                                          and the optimum is
                                                          −x[j]_min)
```

The bordering of `B^{-1}` for the per-variable LP gives a one-row update for
the cost row of the simplex tableau:

```
d̄_{m̃+1}  =  − c_B(B_opt^l) · B̃^{-1}
        =  ( 0, …, 0 )       if   j ∉ Ñ(B)
        =  (specific entries) if  j ∈ Ñ(B)
```

(when the variable being minimised is not in the basis, the reduced costs of
the basic vars are unchanged; when it *is* in the basis, the reduced-cost row
needs the row of `B̃^{-1}` corresponding to where `j` sits in the basis).
The entry of `b̃(B,N)` corresponding to the new row obeys

```
b̃(B,N)[m̃+1]  =  b(B_opt^l, N_opt^l)[ j_i ]   if  j = j_i ∈ N(B_opt^l)
              =  (other case)                  if  j ∈ N_opt^l
```

Run primal/dual simplex on this LP starting from the warm state. Per the
dissertation, typical iteration count is **0** (when `x_opt^l[j]` is already
the minimum) or **1** (when the minimum is across one edge from `x_opt^l`).
Then

```
x_min[j]  =  ] (negated optimum value) [                          (2.34)
```

For `j ∉ J`, set `x_min[j] = 0`.

After processing all `j`, `x_min` is the integer lower-left corner of the
Stage-4 box.

### Stage 4 — provably-optimal lattice search (algorithm step 5)

Run the lattice enumerator (§2) on

```
min  c^p · x^p
s.t. x^p ∈ D^p
     x_min[j]  ≤  x^p[j]  ≤  h^p[j]                                (2.36)
     x^p ∈ N^n
with the strict filter
     c^p · x^p  <  z̃_D^p                                          (2.29)
active from the very first step (КПИА fires from step 0, with
x_∂ := x̃_D^p, z(x_∂) := z̃_D^p).
```

Algorithm step 5 short-circuits:

* **Step 5.1** — if `x_min ≤ x_min^p` componentwise, the Stage-4 region is a
  subset of the (already exhausted) Stage-2 region, so the lattice search
  cannot find anything new. `x_opt^p := x̃_D^p`, jump to step 6.
* **Step 5.2** — otherwise run the lattice enumerator. If `D_min = ∅` (no
  feasible found), `x_opt^p := x̃_D^p`. If a feasible better than `x̃_D^p` is
  found, `x_opt^p := that`.

### Stage 5 — output (algorithm step 6)

Map `x_opt^p` back to the original variables by inverting:

* the re-numbering from §3.1 step 3, and
* the substitution (2.18) `x[j] = h[j] − x^p[j]` on whichever indices it was
  applied to.

Same inversion if the algorithm exited from step 3.3 on tolerance.

---

## 5. Heuristic plug-ins (§2.4)

Architectural commitment:

> "Применение эвристических подходов в комбинированном методе не влияет на
> возможность нахождения строго оптимального решения, а лишь обеспечивает
> более быстрый его поиск."

Heuristics are confined to **Stage 2**: they may *replace* the lattice
enumerator's job of producing the first feasible `x_D^p`, but they may not
influence Stages 1, 3 or 4. The contract:

* The heuristic module receives the input problem (in canonical form (2.17))
  plus the LP optimum `x_opt^l` from Stage 1, and returns either an integer
  feasible `x_D^p` or "no result".
* If it returns a point: treat it exactly as if the lattice enumerator had
  returned it — apply the cheap-improvement step 3.2 and the tolerance check
  3.3, then continue normally.
* If it returns "no result": fall back to the lattice enumerator on (2.22).

For the LAN-topology application Chapter 4 supplies a domain-specific
heuristic; it is not part of the generic library.

КПП (§2.6) is the only heuristic embedded inside the lattice enumerator
itself, and it only reorders forward-step candidates — never prunes — so
optimality survives even with КПП disabled.

---

## 6. Properties (§2.5)

* **Optimality guarantee.** If `D ≠ ∅` and the objective is bounded over the
  integers (granted by the upper bounds `h`), Stages 1–4 return a strict
  optimum. Stage 4's lattice search inside `[x_min, h^p]` with the strict
  filter is the certificate: either it finds a strictly better solution, or
  the exhausted search proves `x̃_D^p` was optimal.
* **Numerical robustness vs. cutting-plane / B&B.** The filter row is a
  *single* extra LP row added once; there is no cumulative cut family that
  would cause coefficient blow-up the way Gomory cuts do. Branching is
  replaced by lattice enumeration in a small box. The dissertation flags
  numerical safety as a major selling point against МПО (правильные
  отсечения = cutting planes) and МВГ (ветвей и границ = B&B).
* **Scaling.** Author's a-priori claim: gap to МПО/МВГ widens with `n`. §3.2
  uses eight crafted small problems exercising every algorithm branch;
  Chapter 4 uses a real LAN-topology instance solved at РКК «Энергия»
  (n = 153, m = 108, full enumeration volume ~2.3·10^143).
* **Heuristic openness.** A user can plug in a Stage-2 heuristic without
  risking optimality.

---

## 7. Software architecture (§3.1)

For reference; useful mainly for choosing module boundaries.

* Language in the dissertation: **C++** (Borland C++ 5.5; also tested with
  GCC under UNIX). Reimplementation in any modern language is fine.
* Module map (§3.1.3, "Схема межмодульных связей"):
  * `IN`     — read problem in the dissertation's input format and echo-print.
  * `OUT`    — print result (`x_opt^p` and `c · x_opt^p`).
  * `MAIN`   — orchestrates the four stages.
  * `TLP`    — bounded-variable simplex *and* the dual bounded-variable
    simplex (the dual variant is needed both inside the regular B&B
    comparator and by Stage 3's per-variable LPs reusing the saved tableau).
  * `CONV`   — transformations to/from canonical forms (§3 above).
  * `SVL`    — lattice enumerator (§2).
  * Heuristic modules slot under Stage 2 with the contract from §5.
* Control flow inside `MAIN` is linear (Stage 1 → 2 → 3 → 4 → 5) but may
  short-circuit (LP infeasible, integer LP optimum, tolerance hit, Stage 4
  empty).
* The bounded-variable simplex code is shared across the comparison
  experiments (B&B, cutting-plane, combined) — same data structure, three
  drivers.

---

## 8. Reference pseudocode

```
function combined_method(c, A, b, sense, h, eps=None, heuristic=None):
    # ---- §3.1   Build (2.17) ----
    if sense == "max":
        c = -c
    flip = []                                       # vars where x = h - x^p substitution happened
    for j where c[j] < 0:
        c[j] = -c[j]
        A[:,j] = -A[:,j]
        b      = b - A[:,j] * h[j]                  # was the OLD A[:,j], careful about order
        flip.append(j)
    perm     = argsort(c)                           # ascending by c
    apply perm to (c, A, h)
    Ap, bp   = split_eq_and_flip_geq(A, b)          # for SVL only

    # ---- §3.2   Build (2.19) ----
    Al, bl, hl, cl = build_LP_canonical(c, A, b, h)

    # ---- Stage 1 ----
    lp_status, x_opt_l, tableau, B_opt_l, N_opt_l = bounded_simplex(Al, bl, hl, cl, maximise=True)
    if lp_status == INFEASIBLE:                     return None
    if is_integer(x_opt_l) and inside(x_opt_l, Ap, bp):
        return decode(x_opt_l, perm, flip, h)

    # ---- Stage 2 ----
    x_min_p = floor(x_opt_l)                        # (2.20)
    box     = HyperBox(low=x_min_p, high=hl_xp_part(hl))   # (2.21)

    if heuristic is not None:
        x_D_p = heuristic(c, A, b, h, x_opt_l)
    else:
        x_D_p = svl_search(Ap, bp, c, box, filter=None)

    while x_D_p is None:
        k = pick_smallest_index_with_box_low_gt_0_and_not_excluded(box)
        if k is None:
            return None                              # D = ∅  (2.25)
        # try slicing on k
        if KH_at_initial_corner_of_slice(box, k):
            box.exclude(k); continue
        box.low[k] -= 1                              # (2.24)
        x_D_p = svl_search_slice(Ap, bp, c, box, k, fixed_value=box.low[k])

    z_D_p = c · x_D_p

    # cheap improvement (step 3.2)
    x_tilde, z_tilde = trivial_improve(x_D_p, z_D_p, c, Ap, bp, hl_xp_part(hl))   # (2.26)

    # tolerance early exit (step 3.3)
    if eps is not None:
        Z_opt_l_in_p_sign = -tableau.objective       # back to (2.17) sign
        gap_upper = (z_tilde - Z_opt_l_in_p_sign) / abs(z_tilde) * 100
        if gap_upper <= eps:
            return decode(x_tilde, perm, flip, h)

    # ---- Stage 3 ----
    add_filter_row_to_tableau(tableau, c, z_tilde,           # (2.30) §2.3.1
                              B_opt_l, N_opt_l)

    x_min = zeros(n)                                          # x_min
    J     = set(N(B_opt_l)) | set(N_opt_l)                    # (2.35)
    for j in J:
        z_j_opt_l = solve_min_xj_warmstart(tableau, j)        # (2.33), §2.3.2
        x_min[j]  = ceil(z_j_opt_l)                           # (2.34)
    # j ∉ J: x_min[j] is already 0

    # ---- Stage 4 ----
    if all(x_min[j] <= x_min_p[j] for j in range(n)):         # step 5.1
        x_opt_p = x_tilde
    else:
        box4 = HyperBox(low=x_min, high=hl_xp_part(hl))
        better = svl_search(Ap, bp, c, box4,
                            filter=StrictFilter(c, z_tilde))   # (2.36) + (2.29)
        x_opt_p = better if better is not None else x_tilde

    # ---- Stage 5 ----
    return decode(x_opt_p, perm, flip, h)
```

```
function svl_search(A, b, c, box, filter):
    # Pre: c[j] ≥ 0, sorted ascending, integer.
    # Searches { x ∈ Z^n : box.low ≤ x ≤ box.high, A x ≤ b },
    # minimising c·x, optionally subject to a filter c·x < f_best.
    n  = len(c); m = len(b)
    x  = box.low.copy()
    y  = b - A @ x                                  # main slacks       (2.6)
    if filter is not None:
        y0 = filter.f_best - c @ x                  # filter slack      (2.11)
    j_x = max((j for j in range(n) if x[j] != 0), default=-1)   # (2.3)

    best = None; f_best = filter.f_best if filter else +inf

    def visit(x, y, y0, j_x):
        nonlocal best, f_best
        if x[n-1] == box.high[n-1]:                 # terminal — caller backs up
            consider_feasibility_then_return()
            return
        if all(yi >= 0 for yi in y):                # feasible
            f = c @ x
            if f < f_best:
                best, f_best = x.copy(), f
                if filter is not None: filter.update(best, f_best)
            return                                  # back up (Rule 1)

        # КН (2.7-2.9)
        I_x = [i for i, yi in enumerate(y) if yi < 0]
        for i in I_x:
            J_xi = [j for j in candidate_J_x(x, j_x, box) if A[i,j] < 0]
            if sum(A[i,j] * (box.high[j] - x[j]) for j in J_xi) > y[i]:
                return                              # KH fires, back up

        # candidate set J_x
        cands = [j for j in candidate_J_x(x, j_x, box)]

        # КПИА (2.10-2.12)
        if filter is not None:
            cands = [j for j in cands if c[j] < y0]

        # КПП (2.13)  — heuristic ordering, optional
        cands.sort(key=lambda j: sum(A[i,j] * (box.high[j] - x[j]) for i in I_x))

        for j in cands:
            x[j] += 1
            y    -= A[:,j]
            if filter is not None: y0 -= c[j]
            visit(x, y, y0, j)
            x[j] -= 1
            y    += A[:,j]
            if filter is not None: y0 += c[j]

    visit(x, y, y0 if filter else None, j_x)
    return best


def candidate_J_x(x, j_x, box):
    # (2.2): j ∈ {1..n : j ≥ j_x}; if j == j_x, require x[j_x] < box.high[j_x]
    n = len(x)
    out = []
    for j in range(max(j_x, 0), n):
        if j == j_x and x[j] >= box.high[j]:
            continue
        if x[j] >= box.high[j]:
            continue
        out.append(j)
    return out
```

---

## 9. Parameters and tunables

* `eps` (`δ_Z`) — relative tolerance for early exit on Stage 2.3 (default:
  not set; method runs to optimum).
* `expand_order` — pluggable rule for which variable Stage 2.1.2 expands on
  next; default is "smallest `c^p[j]` index" (= "smallest index" after sort).
* `heuristic` — pluggable Stage-2 producer of `x_D^p`.
* `enable_KPP` — toggle the КПП ordering inside the lattice enumerator (does
  not change correctness, only running time).
* `simplex_pivot_rule` — affects warm-starts in Stage 3. Bland's rule keeps
  the per-variable LP solves bounded; default is whatever the bounded-variable
  simplex used in `TLP`.

---

## 10. Cross-checks against the eight test cases (§3.2)

The dissertation supplies eight crafted test problems exercising every branch
of the algorithm. When implementing, run all eight and compare:

1. Case `D^p ≠ ∅` while LP-floor `x_min^p` is already feasible → Stage 2.1.1
   succeeds first try.
2. Stage 2.1.2 expansion is needed once.
3. Stage 2.1.2 multi-step expansion: initial slice triggers КН on one
   direction, switch to another variable.
4. Stage 2.2 trivial-improvement is exercised.
5. End-to-end exercise of all four stages.
6. LP optimum is integer → algorithm exits at end of Stage 1.
7. Stage 4 lattice search is short-circuited (`x_min ≤ x_min^p`).
8. No Stage-2 improvement possible (`x̃_D^p ≡ x_D^p`).

§3.2 confirms hand-solved results match machine results for all eight.
Reproducing this table is a strong sanity check before scaling up.

---

## 11. Remaining gaps

After the LibreOffice → docx → pandoc pipeline almost every formula is now
clean. Two tableau cells in §2.3.1 / §2.3.2 had layout-broken renders and
their content is reconstructed from prose:

1. The `b̃(B,N)[m̃+1]` cell after adding the filter row (§2.3.1, end). The
   dissertation's pandoc render is truncated. Reading the bordering algebra
   for bounded-variable simplex, this cell must equal the filter slack at the
   LP optimum, i.e. `z̃_D^p − c^p · x_opt^l[1..n]`. **Verify against test
   case 5.**
2. The `b̃(B,N)[m̃+1]` cell when forming each per-variable LP (§2.3.2, end).
   The dissertation has a piecewise definition keyed on whether `j` is in
   `N(B_opt^l)` (basic) or in `N_opt^l` (at upper bound); the rendered text
   gives only the first branch. The standard bordering for warm-starts
   resolves this either as the corresponding entry of `b(B_opt^l, N_opt^l)`
   (basic case) or as `h^l[j] − b(B_opt^l, N_opt^l)[?]` (upper-bound case).
   **Verify against test case 5.**

For both, the cleanest way to make this exact in code is to *not* hand-roll
the bordering at all: build the augmented tableau (with the filter row and
the new objective) from scratch, then run one round of phase-1 simplex if
needed to land at a basic feasible point. That gives correct numerics by
construction and only costs at most one extra pivot per variable.

The eight test cases of §3.2 are the right oracle: anything reproducing them
matches the dissertation.
