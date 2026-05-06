# Concrete tasks the dissertation runs the Combined Method on

The dissertation evaluates the Combined Method on **two LAN-topology
optimisation tasks** drawn from the РКК «Энергия» computing centre. Both
sit inside the (2.1) shape that the Combined Method ingests after the
canonical rewrite of `COMBINED_METHOD.md` §3 — that is, they are linear
mixed-integer programs with two-sided variable bounds and a mix of `≤`,
`=`, `≥` rows.

The two variants are:

* **Variant A — different hub types, direct connection** (§1.1.1, §4.1).
  Two-level star: users → concentrators → central commutator. Two hub
  types, with stacking allowed for type 1.
* **Variant B — uniform hubs, indirect connection** (§1.1.2, §4.2).
  Single-type concentrators may be daisy-chained; the path from a user
  to the central commutator can pass through several hubs.

Both are stated below with their numerical data and the transformations
the dissertation applies to feed them to the Combined Method, so each
becomes a self-contained input `(c, A, b, sense, h)` that the spec's
`combined_method(...)` entry point can consume.

---

## Variant A — two hub types, direct connection (§1.1.1)

### A.1 Verbal statement (Содержательная постановка)

Build a two-level star LAN connecting users distributed across `n`
rooms, with `b[j]` users in room `j`. Users connect directly to
concentrators (КЦ); concentrators connect to a single central commutator
(ЦК) located in a fixed room `k`. Each concentrator has `a` ports, each
of which is used either for a user link or for the link to the ЦК. Two
hub types are available:

* **Type 1** can be **stacked** in a single room into groups of up to
  `K_max` units. A whole stack uses **one** port for the link to the ЦК
  and the remaining ports for users (so a stack of `s` type-1 hubs in a
  room provides `s·a − 1` user-facing ports and consumes one cable
  segment from the ЦК regardless of `s`).
* **Type 2** does not stack. Each hub of type 2 in a room provides
  `a − 1` user ports and consumes one cable from the ЦК.

The supply of each type is bounded: `K[1]` and `K[2]`. Distances `d[i,j]`
between rooms are given (`d[j,j] = 0`, `d[i,j] = ∞` if the cable cannot
physically be run between rooms `i` and `j`).

**Goal:** minimise total cable length (user→hub plus hub→ЦК).

### A.2 Numerical data (§1.1.1.2)

```
n        = 18                                    rooms
b        = [1, 1, 0, 2, 1, 4, 1, 11, 4, 2, 1, 0, 1, 2, 1, 2, 3, 0]
k        = 18                                    central commutator location
K_max    = 3                                     max stack height for type 1
K        = [3, 5]                                # of available hubs by type
a        = 8                                     ports per hub
```

Inter-room distance matrix `d[i,j]` (symmetric, in metres). Blank cells
mean **no direct cabling allowed** (`d[i,j] = ∞`); diagonal is `0`. The
dissertation lists rooms in five clusters plus the central commutator
in row/column 18:

```
         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
       ┌
   1:  │ 0  30  10  20  10   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞ 105
   2:  │30   0  20  10  40   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  75
   3:  │10  20   0  10  20   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  95
   4:  │20  10  10   0  30   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  85
   5:  │10  40  20  30   0   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞ 115
   6:  │ ∞   ∞   ∞   ∞   ∞   0  10  30   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  30
   7:  │ ∞   ∞   ∞   ∞   ∞  10   0  20   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  40
   8:  │ ∞   ∞   ∞   ∞   ∞  30  20   0   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  60
   9:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   0  20  40  70  80  90   ∞   ∞   ∞  30
  10:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  20   0  20  50  60  70   ∞   ∞   ∞  50
  11:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  40  20   0  30  40  50   ∞   ∞   ∞  70
  12:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  70  50  30   0  10  20   ∞   ∞   ∞ 100
  13:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  80  60  40  10   0  10   ∞   ∞   ∞ 110
  14:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  90  70  50  20  10   0   ∞   ∞   ∞ 120
  15:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   0   ∞   ∞  30
  16:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   0  10  60
  17:  │ ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞   ∞  10   0  70
  18:  │105 75  95  85 115  30  40  60  30  50  70 100 110 120  30  60  70   0
       └
```

Reported size before reduction: 133 vars, 71 main constraints, brute
hyperparallelepiped volume ≈ 1.756·10^57 (§1.1.1.2).

### A.3 Formal model (§1.1.1.3, eqs. 1.1–1.6)

**Decision variables**

* `x[i,j]`, `i,j = 1..n` — number of user-link ports the hubs in room
  `i` use to serve users in room `j` (equivalently: number of user
  cables run room `i` → room `j`). Integer, ≥ 0.
* `y[t,i]`, `t ∈ {1,2}`, `i = 1..n` — number of type-`t` hubs placed in
  room `i`. Integer, ≥ 0.
* `w[i]`, `i = 1..n` — `1` if at least one type-1 hub is placed in room
  `i`, else `0`. Boolean. `w` exists only to linearise the stacking cost
  (one cable from ЦК per stack, regardless of stack height).

**Objective (1.1)**

```
min  Z =  Σᵢ Σⱼ d[i,j] · x[i,j]            user-cable length
       +  Σᵢ d[i,k] · w[i]                 type-1 stack ↔ ЦК cables
       +  Σᵢ d[i,k] · y[2,i]               type-2 hub ↔ ЦК cables
```

**Constraints**

```
(1.2)  Σⱼ x[i,j]  ≤  a · y[1,i] − w[i] + (a − 1) · y[2,i]    ∀ i = 1..n
       — port balance per room: type-1 stack contributes a·y[1,i] − w[i]
         user ports (one port spent on the ЦК link, only if a stack
         exists, hence "−w[i]"); each type-2 hub contributes a − 1.

(1.3)  Σᵢ x[i,j]  =  b[j]                                    ∀ j = 1..n
       — every user has exactly one cable to a hub.

(1.4)  Σᵢ y[t,i]  ≤  K[t]                                    ∀ t = 1, 2
       — supply limit per hub type.

(1.5)  w[i]  ≤  y[1,i]  ≤  K_max · w[i]                      ∀ i = 1..n
       — link w[i] to y[1,i]:  y[1,i] = 0 ⇒ w[i] = 0;
                                y[1,i] > 0 ⇒ w[i] = 1.

(1.6)  x[i,j] ∈ N,  y[t,i] ∈ N,  w[i] ∈ {0, 1}
```

### A.4 Transformations to feed the Combined Method (§4.1.1)

The dissertation applies these reductions before invoking
`combined_method(...)`:

1. **Flatten to single-index variables.** Both `x[i,j]` and `y[t,i]`
   are re-numbered linearly. The natural row-major order of `x[i,j]`
   gives variable indices, then `y[1,·]`, `y[2,·]`, then `w[·]`.
2. **Drop variables that cannot exist.**
   * Drop every `x[i,j]` with `d[i,j] > D_max` (cable too long /
     impossible to lay).
   * Drop every `x[i,j]` whose target room has `b[j] = 0` (no users to
     serve there) — and drop the corresponding row of (1.3).
   * Optionally drop `x[i,j]` and `y[t,i]` for rooms `i` where it is
     either physically impossible or, after inspecting `d[i,·]`,
     clearly suboptimal to place a hub.
3. **Replace the (1.3) equalities by `≥`-inequalities:**

   ```
   Σᵢ x[i,j]  ≥  b[j]    ∀ j   (replaces 1.3)
   ```

   The combination of (1.2) `≤`-rows, the minimisation direction, and
   bounded variables means this `≥` reformulation is *equivalent*: at
   the optimum the user-port supply gets driven down to the demand
   exactly, so equality is achieved automatically. The dissertation
   prefers this because the lattice-search core (МНПВР) only handles
   inequalities — keeping the equality would force the algorithm to
   split it into two rows, doubling the row count for these
   constraints.
4. **Standardise sides:** move all variable terms to LHS, all
   constants to RHS, and multiply rows so RHS is non-negative (needed
   for the simplex's "minimisation of residuals" feasibility-search
   start).
5. **Split (1.5) into two rows** of single-direction inequalities each:

   ```
   y[1,i] − w[i]              ≥ 0
   K_max · w[i] − y[1,i]      ≥ 0
   ```
6. **Per-variable upper bounds `h[ℓ]`.** For the flattened variable
   index `ℓ` corresponding to:

   * `x[i,j]`:                     `h[ℓ] = b[j]`
   * `y[1,i]` (type-1 hub count):  `h[ℓ] = ⌊(b[j]+1)/a⌋ + 1`
     (literal formula from §4.1.1; here `b[j]` is read in the
     room-`i` aggregation sense — total users in `i`'s catchment)
   * `y[2,i]` (type-2 hub count):  `h[ℓ] = ⌊b[j]/(a−1)⌋ + 1`
     (same caveat)
   * `w[i]`:                        `h[ℓ] = 1`

   These are the `h` vector consumed by `combined_method` per the
   spec's §1.

After steps 1–6 the model is in the (2.1) shape
`(min c·x, A x ≤/≥/= b, 0 ≤ x ≤ h, x ∈ N^n)` and the spec's algorithm
applies directly.

> **Heuristic plug-in (§4.1.2).** Variant A also has a domain-specific
> Stage-2 heuristic (the "ranking" algorithm in §4.1.2 — score each
> room by `KПS[i]` and `SD[i]`, pick the highest-rank room, place
> hubs, connect users, repeat). Per the spec's §5, this slots into the
> `heuristic` parameter of `combined_method`. The dissertation
> measured a 34× speed-up using it; see §A.6 below.

### A.5 Reported result (§4.1.3)

Runs on a Pentium-200 MMX:

| Run                         | Wall time   | Stage-2 sub-opt   | True optimum   |
| --------------------------- | ----------- | ----------------- | -------------- |
| Combined Method, no heur    | 5h 40min    | 705               | 615            |
| Combined Method, w/ heur    | 10 min      | 705               | 615            |

* Stage-2 gap to optimum: **14.634 %**.
* Search volume in Stage 2: **0.08 %** of full enumeration.
* Search volume in Stage 4: **0.04 %**.
* МНПВР alone, B&B alone, cutting-plane alone — all failed to finish
  this instance.

The optimum is the LAN topology shown in dissertation Fig. 4-22.

---

## Variant B — uniform hubs, indirect connection (§1.1.2)

### B.1 Verbal statement

Build a star LAN connecting users distributed across `m` rooms, with
`KP[j]` users in room `j`. Room `S` houses the central commutator
(no users in that room). Single-type concentrators are placed in
user-rooms; each has `a` ports, used either for a user or for a link
to **another concentrator** (so the path user → ЦК can pass through
several hubs in series).

Inputs: cable cost `c[i,j]` from room `i` to room `j` (∞ when no cable
can be laid), and per-hub cost `c_k`.

**Goal:** minimise total cost (cables + hubs).

### B.2 Numerical data (§1.1.2.2)

```
m       = 9                                       rooms (incl. ЦК)
S       = 1                                       ЦК is in room 1
KP      = [0, 15, 16, 3, 9, 10, 20, 15, 10]       users per room
c_k     = 20                                      hub unit cost
a       = 8                                       ports per hub
```

Cost matrix `c[i,j]` (cost of laying cable from room `i` to room `j`).
Off-diagonal `0` here means **no connection possible** (`c = ∞`); the
diagonal is conventionally `0`. The dissertation prints:

```
        1   2   3   4   5   6   7   8   9
      ┌
  1:  │ 0  10  20  30  10  10  70   0  10
  2:  │10   0  20  40   0  20  10  20  40
  3:  │20  20   0  10  30  40   0  30  15
  4:  │30  40  10   0  10  10  10   0   0
  5:  │10   0  30  10   0   0  20   0  20
  6:  │10  20  40  10   0   0  30   0  10
  7:  │70  10   0  10  20  30   0  10  20
  8:  │ 0  20  30   0   0   0  10   0  30
  9:  │10  40  15  20  20  10  20  30   0
      └
```

Read the off-diagonal `0`s as `∞` (no cable possible). The matrix is
not strictly symmetric in the printed form — treat any asymmetry as a
typo and use the symmetrised value if implementing.

Reported size before reduction: **153 variables, 108 main constraints**,
brute hyperparallelepiped volume ≈ 2.3·10^143 (§1.1.2.2).

### B.3 Formal model (§1.1.2.3, eqs. 1.7–1.11)

**Graph.** Build the directed graph `(V, D)` where `V` is the room set
(`|V| = m`), `D` is the set of all *possible* arcs (room pairs `i, j`
with `c[i,j] < ∞`, both directions). For each vertex `v`:

* `D_v^+` = arcs *into* `v`.
* `D_v^-` = arcs *out of* `v`.

**Trick for the ЦК.** Set `KP[S] = −KPS` where
`KPS = Σ_{v∈V} KP[v]` (the total users including the negated ЦК).
This makes the ЦК a sink that "absorbs" `KPS` users via the
flow-balance constraint, so a single conservation law works for every
vertex.

**Decision variables**

* `x[d]`, `d ∈ D` — number of users whose traffic uses arc `d` to
  reach the ЦК. Integer ≥ 0.
* `z[d]`, `d ∈ D` — `1` if a cable is physically laid on arc `d`, else
  `0`. Boolean.
* `y[v]`, `v ∈ V` — number of concentrators stacked in vertex `v`. The
  total port count at `v` is `a · y[v]`. Integer ≥ 0.

**Objective (1.7)**

```
min  Σ_{d∈D}  c[d] · z[d]   +   Σ_{v∈V}  c_k · y[v]
       └ cabling cost ┘            └ hub cost ┘
```

**Constraints**

```
(1.8)  Σ_{d ∈ D_v^+} x[d]  −  Σ_{d ∈ D_v^-} x[d]  =  KP[v]      ∀ v ∈ V
       — flow balance: (incoming users) − (outgoing users) =
         users originating at v. Because KP[S] = −KPS, the ЦК
         vertex has a deficit of −KPS, absorbing all traffic.

(1.9)  x[d]  ≤  KPS · z[d]                                       ∀ d ∈ D
       — link x[d] to z[d]: cable must be laid before any user
         traffic flows on arc d. KPS is the trivial big-M.

(1.10) a · y[v]  ≥  KP[v]  +  Σ_{d ∈ D_v^-} z[d]  −  1           ∀ v ∈ V
       — capacity: hub stack at v must have enough ports to serve
         (local users) + (number of outgoing cable links) − 1.
         The "−1" accounts for the single port the stack uses on
         its own outgoing path toward the ЦК.

(1.11) y[v]  ≤  Σ_{d ∈ D_v^-} x[d]  +  KP[v]  −  1               ∀ v ∈ V
       — anti-stub: do not place a hub at a leaf vertex that has
         only one user (no point putting a hub there).

(integrality)
       x[d] ∈ N         (d ∈ D)
       y[v] ∈ N         (v ∈ V)
       z[d] ∈ {0, 1}    (d ∈ D)

(upper bounds)
       x[d]  ≤  KPS                     (d ∈ D)
       y[v]  ≤  ⌊(KP[v] + |D_v^-| − 1) / (a − 2)⌋ + 1
       y[v]  ≤  ⌊(KP[v] + (n − 1) − 1) / (a − 2)⌋ + 1
              (the dissertation gives both, the second is a looser
              fallback that doesn't depend on the realised D_v^-)
```

### B.4 Transformations to feed the Combined Method (§4.2.1)

1. **Flatten** `x[d]`, `z[d]`, `y[v]` to a single index vector. Order
   in the dissertation: `x` block, then `y` block, then `z` block.
2. **Drop unused variables and constraints.**
   * Drop `y[v]` for vertices with no users (and the rows of (1.10),
     (1.11) referencing them).
   * Optionally drop `y[v]` for vertices where placing a hub is
     impossible or, after inspecting `c`, clearly suboptimal.
3. **Replace the (1.8) flow-balance equalities by `≥`-inequalities**
   in the same way as Variant A:

   ```
   Σ_{d∈D_v^+} x[d]  −  Σ_{d∈D_v^-} x[d]   ≥   KP[v]    ∀ v ∈ V
   ```

   At optimum the inequality is tight (any extra flow incurs cable cost
   without serving extra users), so optimum value is preserved while
   row count drops by half versus splitting the equality.
4. **Standardise sides:** all variable terms on the LHS, constants on
   the RHS, signs flipped if needed so RHS ≥ 0.
5. **Per-variable upper bounds `h[ℓ]`.** For the flattened index `ℓ`
   corresponding to:

   * `x[d]`:                         `h[ℓ] = KPS`
   * `y[v]`:                         `h[ℓ] = ⌊(KP[v] + (m − 1) − 1) / (a − 2)⌋ + 1`
   * `z[d]`:                         `h[ℓ] = 1`

After steps 1–5 the model is in the (2.1) shape and the spec's
`combined_method(...)` applies directly. No domain-specific Stage-2
heuristic is reported for Variant B.

### B.5 Reported result (§4.2.2)

Run on Pentium-200 MMX:

| Stat                                      | Value       |
| ----------------------------------------- | ----------- |
| Wall time                                 | 6 hours     |
| Stage-2 sub-optimum objective             | 360         |
| True optimum objective                    | 280         |
| Stage-2 gap                               | 21.6 %      |
| Search volume in Stage 2                  | 0.006 %     |
| Search volume in Stage 4                  | 0.002 %     |

МНПВР, B&B, and cutting-plane all failed to finish on this instance
within the test budget.

---

## Mapping to the Combined Method's input

Both variants reduce to the same data shape after the §A.4 / §B.4
transformations:

```python
problem = {
    "c":     <flat_objective_vector_n>,
    "A":     <constraint_matrix_m_by_n>,
    "b":     <rhs_vector_m>,
    "sense": ["<=", "<=", ..., ">=", "<=", ...],   # per row
    "h":     <upper_bound_vector_n>,               # all integer
}
result = combined_method(**problem, eps=None)        # exact
# or
result = combined_method(**problem, eps=15.0,
                         heuristic=lan_topology_ranking_heuristic)
```

For the dissertation's runs:

* Variant A pre-reduction: 133 vars / 71 rows. Post-reduction (after
  dropping unreachable cables, zero-user destinations, and split (1.5)
  rows) the counts shift; the dissertation does not separately publish
  the post-reduction sizes for Variant A.
* Variant B pre-reduction: 153 vars / 108 rows.

Either is well outside hand computation but ends up tractable for the
Combined Method because the LP relaxation gives a tight Stage-1 corner
and Stage 4 only needs to sweep < 0.1 % of the box.

## Recovering the exact instance for re-implementation

If you intend to reproduce the dissertation's numbers as a regression
test, the four pieces you need are:

1. The numerical tables (§A.2 and §B.2 above) — present here verbatim.
2. The algebraic model (§A.3 / §B.3) — present here verbatim.
3. The reduction rules (§A.4 / §B.4) — present here verbatim.
4. The Variant A heuristic (§4.1.2 of the dissertation, summarised in
   §A.4 above as "ranking heuristic"). The full pseudocode is in the
   dissertation; if you implement it, the speedup vs. plain МНПВР is
   ≈34× per the times in §A.5.

The Combined Method spec in `COMBINED_METHOD.md` consumes any
`(c, A, b, sense, h)` tuple, so the only project-specific work is
building the flat constraint matrix from §A.3 / §B.3 — a routine
mechanical job once the variable index orders are pinned down.
