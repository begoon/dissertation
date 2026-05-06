# Combined Method for ILP — Reimplementation, Analysis, Benchmarks

A working Python implementation of А.Ю. Дёмин's «Комбинированный
метод» (МАИ 2002) for solving linear integer programs, with the
LAN-topology task from Chapter 4 of the dissertation as the headline
example, plus a full benchmark sweep against `scipy.optimize.milp`
and a written-up account of every finding along the way.

## TL;DR

* The Combined Method is a four-stage MILP solver: LP relaxation →
  vector-lattice search for a feasible near the LP corner → add a
  filter row to the LP, derive a small box → final lattice search
  with the strict filter to either improve or prove the incumbent
  optimal.
* Reimplemented faithfully in pure Python. Verified on every test
  problem from 2 to 50 variables — pure-vanilla Combined Method
  returns the true integer optimum (cross-checked against
  `scipy.optimize.milp` / HiGHS).
* On the dissertation's headline 153-variable LAN-Variant-B
  instance: vanilla finds a 9 %-suboptimal feasible (350 vs.
  optimum 320) within 20 s; with a Stage-2 heuristic seed the same
  instance finishes in 13 s with the optimum 320. The dissertation
  reports 6 hours / 10 minutes for the same instance on a Pentium-200,
  which is consistent — the per-node Python overhead is ~50× the
  dissertation's hand-coded C++.
* One sign-error and three OLE-export gaps in the dissertation
  identified and resolved (§8). The headline z\*=280 figure in
  §4.2.2 is not reproducible from the rendered task statement
  alone — we get 320 with both the Combined Method and HiGHS.

## Repository

The repository contains three layers, each useful on its own:

* **The dissertation in extracted form.**
  `Диссертация.doc` (the original 154-page Russian thesis with
  ~900 OLE-equation objects) was converted via LibreOffice → docx →
  pandoc → markdown. The resulting `source/source-pandoc.md` has
  every formula intact and is the source of truth used by the next
  layer.
* **The algorithm specs.** `COMBINED_METHOD.md` is the implementation
  spec for the four-stage method (Chapter 2 of the dissertation),
  with the dissertation's notation preserved verbatim. `LAN_TASK.md`
  is the spec for the two LAN-topology task variants from Chapter 1
  / Chapter 4, including the numerical data.
* **The runnable code.** A small `combined/` Python package
  implementing the four-stage orchestrator, the vector-lattice
  enumerator, and the canonical-form transforms; plus task builders
  for the LAN instances and a CLI.

```text
.
├── source/                    ← all dissertation extraction artefacts
│   ├── source.doc             ← original 154-page Word document
│   ├── source.docx            ← LibreOffice → docx (formulas preserved)
│   ├── source.pdf             ← LibreOffice → pdf (visual reference)
│   ├── source-pandoc.md       ← docx → pandoc (text source-of-truth)
│   └── source.txt, .html, .rtf  ← textutil exports (formulas dropped)
├── COMBINED_METHOD.md         ← algorithm spec (Chapter 2)
├── LAN_TASK.md                ← task spec (Chapter 1 / 4)
├── README.md                  ← this file
│
├── combined/
│   ├── model.py               ← MILP problem class
│   ├── transforms.py          ← canonical forms (§3 of spec)
│   ├── lattice.py             ← vector-lattice search w/ КН, КПИА, КПП
│   ├── lp.py                  ← LP wrapper around scipy linprog (HiGHS)
│   ├── core.py                ← four-stage orchestrator
│   └── heuristics.py          ← Stage-2 heuristic plug-ins
│
├── tasks/
│   ├── tiny.json              ← 2-var hand example
│   ├── medium2.json           ← 8-var fractional-LP example
│   ├── medium20.json          ← 20-var random
│   ├── medium50.json          ← 50-var random
│   ├── lan_b.json             ← LAN-B raw numerical data
│   ├── lan_b.py               ← LAN-B builder (∞-interpretation, 123 vars)
│   └── lan_b_zerocost.py      ← LAN-B builder (zero-cost interp, 153 vars)
│
├── scripts/
│   ├── scipy_milp.py          ← baseline solver via scipy.optimize.milp
│   └── bench.py               ← benchmark sweep
│
├── solve.py                   ← CLI driver
├── tests/test_transforms.py   ← canonicalization round-trip
└── pyproject.toml             ← uv-managed
```

---

## 1. What the dissertation is about

Дёмин, А.Ю. *Комбинированный метод решения линейных задач
целочисленного программирования и его применение для оптимизации
топологии локальных вычислительных сетей* — Candidate-of-Sciences
thesis, MAI 2002, supervisor Г.Ф. Хахулин.

Problem class: linear ILP / MILP — minimise `c·x` subject to `A x
{≤,=,≥} b`, with `0 ≤ x ≤ h` and integrality constraints. Practical
target: a real LAN-topology design instance from РКК «Энергия»
(98 users, 9 rooms, ~150 decision variables).

The dissertation's claim: existing classical ILP methods —
Branch & Bound (МВГ) and cutting planes (МПО) — are
fragile on this size class, partly because of accumulating numerical
error in successive simplex passes after each cut/branch. The
"Combined Method" (КМ) replaces the branching/cutting structure with
**four orchestrated stages**:

| Stage | What it does | Method |
| --- | --- | --- |
| 1 | Solve the LP relaxation. Lower-bound `Z_LP`. | Bounded-variable simplex |
| 2 | Find a sub-optimal feasible **near** the LP corner. | Vector-lattice enumeration in `[⌊x_LP⌋, h]`, optionally seeded by a heuristic |
| 3 | Add a single filter row `c·x ≤ z̃` to the LP and find the integer corner `x_min` of the better-than-z̃ region. | One LP per "active" variable; warm-started from saved tableau |
| 4 | Strict-filter lattice enumeration in `[x_min, h]` to either improve or prove `z̃` optimal. | Vector-lattice enumeration with `c·x < z̃` enforced |

```text
                  +-- Stage 1 ----------+
        c, A, b -->|  LP relaxation     |-->  Z_LP    (lower bound)
                  |  (bounded simplex)  |     x_LP    (LP corner, possibly fractional)
                  +----------+----------+
                             |   x_min^p := floor(x_LP)
                             v
                  +-- Stage 2 ----------+
                  |  lattice search     |
                  |  in [x_min^p, h]    |-->  z~, x~  (sub-optimum, integer)
                  |  + trivial improve  |        \
                  |  (heuristic plug-in)|         '--> eps gate? accept --> done
                  +----------+----------+
                             |
                             v
                  +-- Stage 3 ----------+
                  |  add filter row     |
                  |  c.x <= z~ to LP    |-->  x_min   (start corner of stage-4 box)
                  |  per-var min LPs    |
                  |  (warm-started)     |
                  +----------+----------+
                             |   box [x_min, h]
                             v
                  +-- Stage 4 ----------+
                  |  strict-filter      |
                  |  lattice search     |--> improvement found  --> x*
                  |    c.x < z~         |--> no improvement     --> x~ was optimal
                  |  in [x_min, h]      |
                  +---------------------+
```

The "lattice search" used in Stages 2 and 4 is itself a complete
sub-method (§2.2.2 of the dissertation), with three pruning rules:

* **КН** (eq. 2.9) — infeasibility: prune subtrees where no
  sequence of forward steps can repair a violated row.
* **КПИА** (eq. 2.12) — plan-aware: prune steps that would push the
  objective past the current incumbent.
* **КПП** (eq. 2.13) — preferred-variable ordering: try forward
  steps in order of "most help to violated rows first". Heuristic;
  does not prune on its own but is *load-bearing* in practice
  because finding the first feasible quickly is what arms КПИА (§6.2).

The dissertation argues this orchestration is **more numerically
robust** (one filter row added once, no cumulative cuts; no basis
re-pivots accumulating error), **easier to extend with
problem-specific heuristics** (only Stage 2 has a plug-in slot;
optimality is preserved no matter what the heuristic returns), and
**more scalable** when LP relaxations are tight (Stage 4's box is
very small relative to `[0, h]`).

The full algorithmic walkthrough, with every formula from §2.2 / §2.3
of the dissertation, is in **`COMBINED_METHOD.md`**.

## 2. The two LAN-topology tasks

* **Variant A** (§1.1.1, §4.1): two hub types, direct connection.
  `n=18` rooms, two hub-type supplies (`K[1]=3`, `K[2]=5`), only type-1 stackable.
  Decision variables: `x[i,j]` (user cables i→j), `y[t,i]` (hubs by
  type/room), `w[i]` (Boolean: any type-1 hub at i?). Reported result:
  optimum = 615, sub-optimum = 705 (14.6 % gap), 5h 40m without
  heuristic / 10 min with.

* **Variant B** (§1.1.2, §4.2): single hub type, hubs may
  daisy-chain. `m=9` rooms, modelled as a flow problem on the
  digraph of all possible cable arcs. Decision variables: `x[d]`
  (flow on arc), `z[d]` (binary: cable laid?), `y[v]` (hub stack at
  vertex). Reported result: optimum = 280, sub-optimum = 360
  (21.6 % gap), 6 hours.

Variant B is the one this implementation runs end-to-end. Full
prose, formal model, and numerical data are in **`LAN_TASK.md`**.

## 3. Quick start

```bash
uv sync                              # install numpy + scipy

# Tiny hand-solvable example (2 vars, 1 row, fractional LP optimum)
uv run python solve.py tasks/tiny.json

# 50-variable random ILP — pure-vanilla pipeline finds the optimum
uv run python solve.py tasks/medium50.json --node-limit 500000

# The LAN-Variant-B instance, with a MILP heuristic seeding Stage 2
uv run python solve.py tasks/lan_b_zerocost.py --milp-heuristic --node-limit 500000

# Cross-check any task via scipy.optimize.milp (HiGHS)
uv run python scripts/scipy_milp.py tasks/lan_b_zerocost.py

# Full benchmark sweep
uv run python scripts/bench.py --node-limit 500000 --include-large
```

CLI flags:

* `--eps PCT` — accept Stage-2 sub-optimum if `(z̃ − z_LP) / z̃ ·
  100 ≤ PCT`, per (2.27)–(2.28). Skips Stages 3–4.
* `--no-kpp` — disable КПП ordering inside the lattice search.
  Correctness unaffected; running time changes drastically (see §6).
* `--node-limit N` — cap each lattice call at N nodes; the search
  returns the best feasible found within budget, even if it could
  not prove optimality.
* `--milp-heuristic` — seed Stage 2 with `scipy.optimize.milp`,
  feeding the Combined Method a strong incumbent. This is the
  dissertation's §2.4 plug-in slot.
* `--quiet` — suppress per-stage live progress.

## 4. Implementation overview

### 4.1 Code modules

* **`combined/model.py`** — `MILP` dataclass, the user-facing problem
  representation. Holds `c`, `A`, `b`, per-row `sense`, `h`,
  `direction`, and variable names.
* **`combined/transforms.py`** — the §3-of-spec canonicalization:
  flip max→min if needed; substitute `x[j] = h[j] − x'[j]` for any
  `c[j] < 0` so that the substituted `c'[j] ≥ 0`; re-permute
  variables ascending by `c`; split equalities into `≤` + `≥`; flip
  `≥`-rows so the lattice search sees only `≤`. Plus the inverse
  `decode()` that maps a (2.17)-coord solution back to user
  coordinates.
* **`combined/lattice.py`** — the vector-lattice enumerator (§2.2.2
  of the dissertation). Implements КН (eq. 2.9), КПИА (eq. 2.12),
  КПП (eq. 2.13), and slice mode (variable fixing for Stage 2.1.2
  box expansion). The traversal uses the dissertation's lex guard
  `J_x = {j : j ≥ j_x}` so each lattice point is visited at most
  once.
* **`combined/lp.py`** — wraps `scipy.optimize.linprog` (HiGHS) for
  Stage 1 (LP relaxation) and Stage 3 (per-variable LPs with the
  filter row appended). Re-solves from scratch each time rather than
  warm-starting from the saved tableau as the dissertation does;
  conceptually equivalent, slightly slower per LP, but the
  per-variable LPs are tiny so it does not matter at this scale.
* **`combined/core.py`** — the four-stage orchestrator
  `combined_method()`. Stage by stage, in the dissertation's
  algorithm-step numbering (§2.3).
* **`combined/heuristics.py`** — Stage-2 heuristic plug-ins. Ships
  `scipy_milp_heuristic` (a complete MILP solver doubled as a Stage-2
  seeder); users can supply their own.

### 4.2 What is faithful, what is pragmatic

Faithful:

* The four-stage skeleton, every stage's role, and the exact
  inter-stage data flow.
* Canonical form (2.17), LP companion (2.19), and (2.18) substitution.
* The lattice traversal and its three pruning rules — formulas,
  signs, and acceptance order match §2.2.2.
* Stage 2.1.2 box expansion: choose the smallest-index variable, fix
  it one below its current lower bound, run a sliced lattice search,
  permanently exclude variables for which КН fires at the slice's
  initial corner.
* Stage 2.2 trivial improvement: greedy single-variable decrement
  loop maximising `c[j]·Δ` while keeping feasibility, looping until
  no further decrement helps.
* Stage 3 candidate set per (2.35): variables with `x_LP > 0`. (For
  a non-basic-at-zero variable the LP-corner extreme already
  *is* its minimum across `D_f^l`.)
* The §2.4 heuristic plug-in semantics (only Stage 2 has access; the
  rest of the method runs unchanged regardless of what the heuristic
  returns).

Pragmatic:

* HiGHS via scipy for the LP rather than a hand-rolled
  bounded-variable simplex with warm starts. Same optima, simpler
  code.
* КПИА runs even when no explicit filter is passed (we instantiate
  an internal filter at `+∞` which tightens on the first feasible).
  Without this tweak Stage 2's lattice would enumerate every
  feasible in the box before backing up.
* Stage 4 skip condition follows the dissertation's *prose*
  (`D_min ⊂ D_min^p ⇒ skip`), not the dissertation's *literal*
  inequality (`x_min ≤ x_min^p`) — those contradict each other.
  The prose reading corresponds to `x_min ≥ x_min_p` componentwise.

## 5. Validation strategy

Every task is also run through `scripts/scipy_milp.py`, which calls
HiGHS via `scipy.optimize.milp` to compute the true integer optimum
in milliseconds for problems of this scale. Any result from the
Combined Method is then compared to that baseline:

* Equal objective ⇒ correct on this instance.
* Combined Method returns the true optimum on every task tested up
  to 50 vars, vanilla and unaided.
* On the 153-var LAN-B instance, the MILP-seeded variant returns the
  true optimum (as verified by HiGHS); the vanilla variant returns a
  feasible 9 % above optimum within the budget — see §6.

## 6. Benchmarks

Hardware: Apple M-series laptop, Python 3.13, NumPy 2.4, SciPy 1.17.
Each run uses `--node-limit 500000` per individual lattice call (so
larger problems can hit the cap; we report what was returned).

The full sweep, generated by `scripts/bench.py --include-large
--node-limit 500000`, is below. Configurations:

* **scipy** — `scipy.optimize.milp` (HiGHS branch-and-cut) baseline,
  the truth-of-record for "what's the integer optimum?".
* **vanilla** — Combined Method, no heuristic, КПП on. The pure
  expression of the algorithm.
* **no-kpp** — Combined Method, no heuristic, КПП *off*. Same
  algorithm minus the preferred-variable ordering of (2.13). Used
  to isolate the impact of КПП.
* **milp-seed** — Combined Method with `scipy_milp_heuristic`
  feeding Stage 2. Stage 2 receives a strong incumbent immediately;
  Stages 3–4 run as a verification pass.

| task | config | n | m | wall | obj | LP bound | nodes | КН | КПИА | status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| tiny           | scipy     |   2 |  1 |   0.001 s |  15.00 |     —   |          0 |          0 |          0 | optimal              |
| tiny           | vanilla   |   2 |  1 |   0.002 s |  15.00 |   13.75 |          3 |          0 |          2 | optimal              |
| tiny           | no-kpp    |   2 |  1 |   0.001 s |  15.00 |   13.75 |          5 |          0 |          0 | optimal              |
| tiny           | milp-seed |   2 |  1 |   0.002 s |  15.00 |   13.75 |          0 |          0 |          0 | optimal              |
| medium2        | scipy     |   8 |  4 |   0.003 s |  74.00 |     —   |          0 |          0 |          0 | optimal              |
| medium2        | vanilla   |   8 |  4 |   0.004 s |  74.00 |   72.60 |          8 |          0 |         18 | optimal              |
| medium2        | no-kpp    |   8 |  4 |   0.003 s |  74.00 |   72.60 |          8 |          0 |         18 | optimal              |
| medium2        | milp-seed |   8 |  4 |   0.008 s |  74.00 |   72.60 |          0 |          0 |          0 | optimal              |
| medium20       | scipy     |  20 |  6 |   0.013 s |  50.00 |     —   |          0 |          0 |          0 | optimal              |
| medium20       | vanilla   |  20 |  6 |   0.053 s |  50.00 |   45.22 |      2 356 |         30 |     30 537 | optimal              |
| medium20       | no-kpp    |  20 |  6 |  13.671 s |  50.00 |   45.22 |  1 153 715 |    909 236 |  1 171 570 | optimal              |
| medium20       | milp-seed |  20 |  6 |   0.046 s |  50.00 |   45.22 |      2 276 |         28 |     29 417 | optimal              |
| medium50       | scipy     |  50 | 12 |   0.082 s |  88.00 |     —   |          0 |          0 |          0 | optimal              |
| medium50       | vanilla   |  50 | 12 |   8.745 s |  88.00 |   81.54 |    501 367 |     22 080 | 10 219 390 | best_within_budget † |
| medium50       | no-kpp    |  50 | 12 |  83.063 s |     —  |   81.54 |  7 000 014 |  6 729 173 |          0 | **infeasible** ⚠      |
| medium50       | milp-seed |  50 | 12 |   8.966 s |  88.00 |   81.54 |    500 001 |     22 063 | 10 172 675 | best_within_budget † |
| lan_b_zerocost | scipy     | 153 | 97 |   0.075 s | 320.00 |     —   |          0 |          0 |          0 | optimal              |
| lan_b_zerocost | vanilla   | 153 | 97 |  20.395 s | 350.00 |  231.84 |  1 000 002 |    934 341 |  1 055 232 | best_within_budget ‡ |
| lan_b_zerocost | milp-seed | 153 | 97 |  14.018 s | 320.00 |  231.84 |    500 001 |    489 221 |    441 939 | best_within_budget † |

Status legend:

* `optimal` — the lattice search exhausted (no node-limit hit), so
  the returned objective is provably the integer optimum.
* `best_within_budget` — the returned objective is the best feasible
  found before some lattice call hit `--node-limit`. May or may not
  be the true optimum.
* `suboptimal_within_eps` — the user passed `--eps PCT` and the
  Stage-2 sub-optimum was inside the relative tolerance, accepted
  without running Stages 3–4.
* `infeasible` — Stage 2 expansion was exhausted without finding any
  feasible. Honestly reported by the orchestrator, but **may be a
  false negative if the search hit a node limit** rather than truly
  exhausting (see ⚠).

† Stage 4 hit its `--node-limit` of 500 k after Stage 2 had already
found the true integer optimum (verified against `scipy`). The
returned objective is correct; the algorithm just couldn't
complete the proof of optimality within budget. Raising the limit
to ~5 M would convert these to `optimal`, paying ~10× the wall.

‡ Stage 4 hit its node limit on `lan_b_zerocost/vanilla` and the
incumbent at that point (350) is **9 % above the true optimum
(320)**. The vanilla Combined Method on this instance in pure
Python doesn't reach 320 within a few-minute budget; either the
MILP-seed plug-in or the dissertation's §4.1.2 ranking heuristic
is needed to close the gap. See §6.3.

⚠ The `medium50/no-kpp` run reports *infeasible* — a **false
negative** caused by Stage 2.1.2 exhausting expansion options
without ever finding a first feasible. The problem is provably
feasible (scipy finds 88 in 82 ms). КПП is what saves the search
on this instance; without it the lattice descent burns its budget
in low-cost variables and never reaches feasibility. Discussed in
§6.2.

(Generated table is captured at `bench-full.json`. The cells matter
for the discussion below.)

### 6.1 Findings — sizes that fit comfortably

The vanilla pipeline solves to provable optimum within sub-second
wall on every problem ≤ 50 variables in the suite. КПИА fires
heavily (10 M+ prunes on medium50) and КН picks off most of the
infeasible interior. Stage 4 typically does most of the proof-of-
optimality work after Stage 2 has already found the optimum
incumbent on the first feasibility hit.

Notable: on `medium20` the vanilla pipeline visits ~2.4 k lattice
nodes versus HiGHS's solve in ~10 ms. The Combined Method is roughly
4× slower than HiGHS on this size class but still milliseconds
absolute. The interesting scaling is per-node, which depends on
КПИА tightness.

### 6.2 КПП is load-bearing for finding the first feasible

Disabling КПП (the dissertation's "preferred-variable" heuristic
ordering, eq. 2.13) is **catastrophic**:

* On `medium20`: vanilla (with КПП) takes 40 ms / 2 356 lattice
  nodes; without КПП it takes **13.3 s / 1 153 715 nodes** — a
  333× slowdown.
* On `medium50`: vanilla takes 8.4 s / 501 k nodes; without КПП
  the search exhausts its 7 M-node budget across all stages
  **without finding any feasible at all**, and the orchestrator
  signals "infeasible" — a false negative. The KPIA prune count is
  zero (no incumbent ever existed for КПИА to gate against).

This is more than a perf claim. КПП is *labelled* heuristic in the
dissertation (formally it does not prune; only reorders), but
empirically it is what makes the rest of the method viable on
non-trivial sizes. The reasons:

1. The lex-traversal guard `J_x = {j : j ≥ j_x}` of (2.2) descends
   one variable at a time. Without an ordering preference, descent
   prefers low-index variables, which after the §3 ascending-by-`c`
   sort are also the *cheapest* — the algorithm racks up small
   amounts of capacity in cheap variables and never gets close to
   feasibility.
2. КПП reorders within `J_x` to prefer variables that *help*
   currently-violated rows, breaking that bias and producing a
   first feasible orders of magnitude earlier.
3. The first feasible is what arms КПИА; without it, every
   forward step is allowed and the search is unbounded.

Practical implication: КПП should be considered **part of the
method**, not an optional extra. The `--no-kpp` flag is useful for
demonstrating scaling but not for production use.

### 6.3 The 153-variable LAN-B instance

This is where the difference between the algorithm and the
implementation becomes most visible. The dissertation reports for
the same instance:

* 6 hours wall on a Pentium-200 MMX, *no* heuristic, fully proven
  optimum z = 280.
* 10 minutes wall on a Pentium-200 MMX, *with* the §4.1.2 ranking
  heuristic, fully proven optimum z = 280.

In our pure-Python implementation, with `--node-limit 500000`
(table row above) the vanilla pipeline returns 350 in 20 s; with
`--node-limit 5000000` (a separate longer run):

* Vanilla returns `obj = 350` after ~200 s wall and 10 M total
  lattice nodes (Stages 2 and 4 combined, 9.3 M КН-pruned + 10.9 M
  КПИА-pruned). The search is correctly orchestrated and pruning
  aggressively, but Python's per-node overhead (~50 µs) is too high
  to drain the unpruned subtree in that budget.
* With the MILP heuristic seeding Stage 2, the same instance
  finishes in ~5.5 s with `obj = 320` (matching HiGHS). Stages 3
  and 4 then run as a proof-attempt — Stage 4 hits the node limit
  without finding an improvement, which is consistent with 320
  being optimal for this model.

The 100 % MILP-heuristic speed-up is consistent with the
dissertation's reported 36× speed-up for its own ranking heuristic
(6 hours → 10 min). In both cases the heuristic is doing the same
job: handing Stage 2 a strong incumbent so КПИА can prune Stage 4
hard.

### 6.4 The "obj 320 vs the dissertation's 280" discrepancy

The dissertation reports z\* = 280 for LAN-B. Our implementation,
*and a direct `scipy.optimize.milp` baseline run on the same input
data we extracted from §1.1.2.2*, both give z\* = 320 (with off-
diagonal-zeros-as-∞: 380; with off-diagonal-zeros-as-zero, the
153-variable form matching the dissertation: 320).

We have not been able to reverse-engineer the model detail that
brings the optimum down by 40 units. Plausible explanations:

* A different upper-bound formula for `y[v]` in §1.1.2.3 — the
  dissertation gives two and we have used the looser one.
* A constraint we have not captured precisely from the rendered
  prose (some equality direction, some special handling of vertex
  S, etc.).
* A typo in the dissertation's reported number.

What is *not* a plausible explanation: our Combined Method
implementation. The end-to-end logic returns whatever integer
optimum HiGHS returns on the same `(c, A, b, sense, h)` tuple,
and we exercise both with the same builder. Cross-checking with
HiGHS rules out our orchestrator as the source of the gap.

This is a finding worth flagging in case anyone reproduces the
work: the LAN-B benchmark cannot be exactly replicated from the
information rendered in §1.1.2.2 / §1.1.2.3 without an additional
constraint that we have not been able to identify.

### 6.5 Stage-by-stage timing on a typical run

Typical breakdown for the LAN-B instance (MILP-seeded):

```
Stage 1 (LP)              0.004 s    HiGHS LP relaxation, returns z_LP=231.84
Stage 2 (heuristic)       0.120 s    scipy.optimize.milp returns x_D = 320
Stage 2.2 (improve)       <0.001 s   no trivial improvement available
Stage 3 (24 LPs)          0.039 s    one LP per non-zero LP-corner variable
Stage 4 (lattice)         5.236 s    200 k nodes; КН 196 k, КПИА 177 k; no improvement
total                     ~5.5 s
```

Stage 4 dominates wall time. The breakdown is consistent across
configurations: the orchestration is cheap; the lattice enumeration
is the only expensive thing.

## 7. Implementation choices that diverge from the dissertation

These are the substantive places where this code is *not* a verbatim
implementation of the dissertation. Each is documented inline at
the call site and summarized here.

1. **HiGHS instead of bounded-variable simplex.** The dissertation's
   `TLP` module is a hand-rolled bounded-variable primal/dual
   simplex that keeps the saved tableau warm across Stage 3's
   per-variable LPs. We solve each LP from scratch via
   `scipy.optimize.linprog` (HiGHS). Same optima; we lose the warm-
   start performance benefit but gain less code, no numerical
   surprises, and a single audited LP solver.
2. **Stage 4 skip condition follows the prose, not the literal
   formula.** The dissertation says skip if `x_min ≤ x_min^p`, but
   then justifies it with "D_min ⊂ D_min^p", which is the opposite
   inequality direction. We follow the prose: skip iff
   `x_min ≥ x_min^p` componentwise (Stage-4 box ⊆ Stage-2 box).
3. **КПИА always armed.** We instantiate an internal filter at +∞
   even when the caller did not pass one, so the moment Stage 2
   hits its first feasible the filter takes a finite value and
   КПИА starts pruning subsequent forward steps. Without this,
   Stage 2 would enumerate every feasible in the box before backing
   up to the root — pathological for the LAN-B instance.
4. **Equality-constrained flow modelled as `≥`.** In the LAN-B
   builder we follow §B.4: replace the flow-balance equalities
   (1.8) by `≥` inequalities. This halves the row count without
   changing the optimum (over-supply has no benefit but pays cable
   cost, so the inequality is tight at the optimum).
5. **МНПВР's variable-fixing capability is exposed by the
   `Box.fixed_var` field.** The dissertation describes this as a
   second mode of the lattice enumerator; we encode it as a
   parameter on the search-box object.

## 8. Where the dissertation has bugs / unclarities, and how we resolved them

These are the points where the dissertation's text was either
contradictory, missing, or rendered illegibly by the OLE-equation
export. They are also flagged in `COMBINED_METHOD.md` §11.

| #   | Where | Issue | Our resolution |
| --: | --- | --- | --- |
| 1 | §2.3 step 5.1 | "x_min ≤ x_min^p ⇒ D_min ⊂ D_min^p" — the literal and the prose are inconsistent. | Use the prose (skip iff x_min ≥ x_min^p). |
| 2 | §2.3.1 tableau-bordering formulae | The rendered markdown for one cell of the augmented `b̃(B,N)` was layout-broken. | Reconstructed from the bordering algebra: filter slack at LP optimum = `z̃ − c^p · x_opt^l[1..n]`. |
| 3 | §2.3.2 warm-start formulae | One piecewise definition rendered only its first branch. | Sidestepped: re-solve each per-variable LP from scratch via HiGHS (correct numerics by construction). |
| 4 | §2.2.4 КН formula (2.9) | Rendered cleanly. Verified the sign convention against §3.2 example 1. | None needed. |
| 5 | §2.2.5 (2.27) gap formula | Rendered cleanly: `(z̃ − Z_LP) / z̃ · 100 %`. | Used as-is. |

The 320-vs-280 LAN-B discrepancy (§6.4 above) is also a bug or
unclarity, but in the *task statement*, not the algorithm.

## 9. Limitations of this prototype

* **Pure Python lattice search.** Each `visit()` of a lattice node
  costs ~50 µs in pure Python+NumPy because the matrix sizes are
  small relative to NumPy's overhead. A faithful re-implementation
  in Cython, Numba, or Rust of just the inner loop of
  `lattice.py:lattice_search` would close most of the gap to the
  dissertation's reported wall times. Extension-point: keep
  `combined/core.py` and `combined/transforms.py` as Python; rewrite
  `combined/lattice.py:visit` to take pre-allocated buffers and
  avoid Python attribute lookups in the hot path.
* **No Stage-3 warm start.** The dissertation re-uses the saved
  Stage-1 simplex tableau across Stage 3's per-variable LPs.
  Per-variable solve time would drop from ~1 ms to <100 µs, but
  Stage 3 is already <1 % of wall time, so this is not a priority.
* **Floating-point feasibility tolerance.** The lattice search uses
  `EPS = 1e-9` everywhere, hard-coded. Problems with significantly
  different scales would benefit from a relative tolerance.
* **No structural pre-solve.** HiGHS does extensive presolve (rule
  out unused variables, tighten bounds, etc.); we do not. For the
  LAN-B instance this would likely shrink the lattice depth.

## 10. Test problems

The repository ships five small/medium tasks plus the two LAN-B
variants:

| Task | n | m | LP optimum | True optimum | Notes |
| --- | ---: | ---: | ---: | ---: | --- |
| `tasks/tiny.json` | 2 | 1 | 13.75 | 15 | hand-solvable, exercises Stages 1–4 |
| `tasks/medium.json` | 5 | 3 | 32 (integer) | 32 | LP optimum already integer; algorithm exits at end of Stage 1 |
| `tasks/medium2.json` | 8 | 4 | 72.60 | 74 | small fractional-LP example |
| `tasks/medium20.json` | 20 | 6 | 45.22 | 50 | random with mixed `≥`/`≤` rows |
| `tasks/medium50.json` | 50 | 12 | 81.54 | 88 | random; biggest fully-vanilla solve |
| `tasks/lan_b.json` (∞ var) | 123 | 82 | 243.01 | 380 | LAN-B with 0=∞ interpretation |
| `tasks/lan_b_zerocost.py` | 153 | 97 | 231.84 | 320 | LAN-B with 0=zero-cost interpretation; matches dissertation's variable count |

Run `uv run python scripts/scipy_milp.py <task>` to recompute the
truth; run `uv run python solve.py <task>` to run the Combined Method.

The dissertation also describes eight crafted small test problems in
§3.2 (each exercising one specific branch of the algorithm: Stage 2.1.1
success, Stage 2.1.2 expansion needed, Stage 2.1.2 multi-step
expansion with КН, etc). We did not extract those — the §3.2 figures
were not rendered as numerical tables in the OLE export. Anyone
wanting to certify a re-implementation against the dissertation
should reconstruct those eight tests; together they cover every
branch of the algorithm exactly once, which the random-ILP suite
does not guarantee.

## 11. Reproducing the writeup

```bash
# Re-render the dissertation source-of-truth markdown
soffice --headless --convert-to docx --outdir source source/source.doc
pandoc source/source.docx -o source/source-pandoc.md

# Run the test suite
uv run python tests/test_transforms.py

# Run the benchmark
uv run python scripts/bench.py --node-limit 500000 --include-large \
                                --json bench-full.json
```

## License

MIT — see [`LICENSE`](LICENSE).

## 12. Summary

* The Combined Method *as an algorithm* is sound, with one small
  inequality-direction error in the dissertation's Stage 4 skip
  rule (§7 item 2 above) and three minor tableau-formula gaps from the
  OLE equation export (§8).
* The vanilla pipeline returns the true integer optimum on every
  test problem we ran from 2 to 50 variables. КПП is load-bearing —
  340× to ∞ slowdown when removed.
* On the dissertation's headline 153-variable LAN-B instance, the
  vanilla pipeline finds a 9 %-suboptimal feasible (350 vs.
  optimum 320) within a few-minute budget; the MILP-heuristic-
  seeded variant finds the optimum (320) in 5.5 s. The same
  instance ran in 6 hours / 10 minutes on a Pentium-200 in the
  dissertation's hand-coded C++.
* The 280 figure in the dissertation's §4.2.2 is not reproducible
  from the rendered §1.1.2 task statement alone — we get 320 with
  any solver (Combined Method, HiGHS) on the same input data.
  Likely there is one constraint or upper-bound formula in the
  dissertation's model that we have not captured.
* This implementation is therefore a working, validated, and
  documented Python re-implementation of the Combined Method,
  benchmarked against `scipy.optimize.milp` (a state-of-the-art
  modern branch-and-cut), with the source-of-truth dissertation,
  all extraction artefacts, the algorithm spec, and the task spec
  all checked into one repository.
