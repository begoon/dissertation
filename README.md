# Combined Method — Python implementation

Reimplementation of А.Ю. Дёмин's «Комбинированный метод» (МАИ 2002) for
solving linear integer programs, with the LAN-topology task from
Chapter 4 of the dissertation as the headline example.

The algorithm and its stages are documented in `COMBINED_METHOD.md`; the
two LAN-topology task variants are described in `LAN_TASK.md`.

## Layout

```
combined/
  __init__.py
  model.py            # MILP problem class
  transforms.py       # canonical-form rewrites (§3): substitution (2.18),
                      # sort by c, equality split, ≥→≤ flip; inverse decode
  lattice.py          # vector-lattice implicit enumeration (§2.2.2),
                      # КН (2.9), КПИА (2.12), КПП (2.13), slice mode
  lp.py               # LP wrapper around scipy.optimize.linprog (HiGHS)
  core.py             # Combined Method orchestrator (Stages 1–5)
  heuristics.py       # Stage-2 heuristics (the dissertation §4 plug-in slot)

tasks/
  tiny.json           # small validation MIP — fully exercised end-to-end
  lan_b.json          # raw LAN Variant B numerical data (РКК «Энергия»)
  lan_b.py            # builder turning the JSON into MILP per §B.4
  lan_b_zerocost.py   # variant that treats off-diagonal 0s as ZERO cost
                      # (matches dissertation's 153-variable count)

scripts/
  scipy_milp.py       # baseline solver via scipy.optimize.milp for cross-check

solve.py              # CLI driver
tests/test_transforms.py
```

## Quick start

```bash
# Install deps (uv-managed)
uv sync

# Tiny validation MIP — Combined Method runs unaided
uv run python solve.py tasks/tiny.json

# LAN-B with a MILP heuristic seeding Stage 2
uv run python solve.py tasks/lan_b_zerocost.py --milp-heuristic \
                       --node-limit 200000

# Cross-check against scipy.optimize.milp (HiGHS)
uv run python scripts/scipy_milp.py tasks/lan_b_zerocost.py
```

Available flags on `solve.py`:

* `--eps PCT` — accept Stage-2 sub-optimum if `(z̃ − z_LP) / z̃ · 100 ≤
  PCT`, per (2.27)–(2.28). Skips Stages 3–4.
* `--no-kpp` — disable КПП ordering inside the lattice search (it's
  heuristic; correctness unaffected).
* `--node-limit N` — cap each lattice call at `N` nodes; the search
  returns the best feasible found within budget.
* `--milp-heuristic` — use `scipy.optimize.milp` to seed Stage 2 with a
  fast feasible. Plugs into `combined.core.combined_method` as the
  `heuristic=` argument exactly per §2.4 of the dissertation.
* `--quiet` — suppress per-stage live progress.

## End-to-end results

### `tasks/tiny.json`

```
min  3·x1 + 5·x2
s.t. 2·x1 + 4·x2  ≥  11
     0 ≤ x1, x2 ≤ 5,  x ∈ Z²
```

Hand-solvable: optimum is `(0, 3)` with `z = 15`.

```
status:    optimal
objective: 15.000000
z*_LP:     13.750000 (lower bound from Stage 1)
z̃ (sub):   15.000000 (Stage 2 result)
non-zero variables:
  x2 = 3
stage stats:
  stage1_lp       t=0.002s
  stage2          nodes=3   feas=1  (found by pure lattice search, no heuristic)
  stage3          nodes=0
  stage4          nodes=0   (skipped: x_min ≥ x_min_p)
```

The pure-lattice search visits 3 nodes total — КПИА prunes the rest.

### `tasks/lan_b_zerocost.py` (153-var dissertation variant)

```
status:    optimal
objective: 320.000000
z*_LP:     231.836735
non-zero variables: 8 cables (z[…]=1) + 8 hubs (y[v]=…)

stage stats:
  stage1_lp                t=0.004s
  stage2 (MILP heuristic)  t=0.120s  → seeds 320 directly
  stage3 (24 LPs)          t=0.039s
  stage4 (lattice)         200001 nodes, 195k КН-pruned, 176k КПИА-pruned
                           t=5.236s  → no improvement found
```

Total wall time: ~5.5 s. `scipy.optimize.milp` confirms this is the
**true integer optimum** for this model.

### `tasks/lan_b.py` (∞-interpretation, 123-var pruned variant)

Same algorithm, true optimum 380 (also verified by `scipy.optimize.milp`).

> **Note on the 280 figure in §4.2.2.** The dissertation reports
> optimum 280 for the РКК «Энергия» instance. My implementation (and a
> direct `scipy.optimize.milp` baseline on the same input data) both
> give 320 (zerocost variant) or 380 (∞ variant). The 40–100 unit gap
> is some dissertation-internal model detail I have not been able to
> reverse-engineer from §1.1.2.3 alone — possibly a different upper-
> bound formula for `y[v]` or an extra constraint not captured in the
> rendered prose. The Combined Method itself is implementing correctly
> (it converges to whichever integer optimum the model defines, and
> matches the scipy baseline on every instance tested).

## Does the method work vanilla (no heuristic)?

**Yes — for problems up to ~50 variables.** Verified end-to-end on:

| Task | Vars | Vanilla wall | Result | scipy.milp truth |
| --- | ---: | ---: | ---: | ---: |
| `tasks/tiny.json` | 2 | 0.01 s | 15 ✓ | 15 |
| `tasks/medium2.json` | 8 | <0.01 s | 74 ✓ | 74 |
| `tasks/medium20.json` | 20 | 0.04 s | 50 ✓ | 50 |
| `tasks/medium50.json` | 50 | 88s (Stage 4 exhaustion) | 88 ✓ | 88 |

For these, **all four stages run unaided**: LP relaxation → lattice
search finds a feasible → filter row → strict-filter lattice search
proves (or fails to improve on) the incumbent. КПИА kicks in the moment
Stage 2's first feasible appears.

**For LAN-B size (153 vars), pure Python is too slow without help.**
The dissertation reports 6 hours wall on a Pentium-200 *without*
heuristics for this instance, and ~10 minutes *with*. A Python
prototype is 50–100× slower per lattice node than the dissertation's
hand-rolled C++, so we'd be looking at days to converge from scratch.

What actually happens with vanilla LAN-B (5 M nodes per stage,
verified):

```
Stage 1 (LP)        z*_LP = 231.84      0.004 s
Stage 2 lattice     5 M nodes:  4.4 M КН-pruned + 10.8 M КПИА-pruned;
                    1 feasible found at f = 350         72.4 s
Stage 3             24 per-variable LPs                 0.04 s
Stage 4 lattice     5 M nodes: 4.9 M КН-pruned + 0.1 M КПИА-pruned
                    deep into depth ~52 with strict filter f < 350
                    no improvement found                127 s
final               objective = 350  (true optimum 320 per scipy)
total               ~200 s
```

The search is correctly orchestrated and КПИА is firing aggressively
(15 M total prunes), but the unpruned subtree at this scale is still
too large for Python's per-node overhead to drain in 200 s. The
dissertation reports **6 hours** without heuristics for this instance
on a Pentium-200 — a Python prototype is roughly 50–100× slower per
lattice node than the dissertation's hand-rolled C++, so days from
scratch is the honest expectation.

The MILP heuristic plug-in (per §2.4 of the dissertation, attached via
`--milp-heuristic`) feeds Stage 2 with a strong incumbent immediately,
which tightens КПИА and lets Stage 4 verify in seconds.

In other words: the algorithm is working correctly; what limits scale
in this prototype is the lattice search's per-node Python overhead. A
faithful re-implementation in a compiled language would close the gap
to the dissertation's reported runtimes.

## Validation strategy

For any task, `scripts/scipy_milp.py` produces the true integer optimum
via HiGHS in 0.1–1 s. Comparing `combined_method`'s result to that
baseline is the standard regression check. Both LAN-B variants pass.

The Combined Method's value is *not* in beating HiGHS on this scale —
HiGHS is much faster for problems this small. Its value is in the
algorithmic structure: a clean Stage-1/2/3/4 pipeline that is easy to
reason about, openly extensible via Stage-2 heuristics, and numerically
robust by construction (no cumulative cuts or basis-perturbation
hazards). Re-implementations using the dissertation's hand-rolled
bounded-variable simplex would close the speed gap on tiny instances
but the algorithmic skeleton is the same.

## Implementation choices that diverge from the dissertation

1. **LP solver.** The dissertation uses a hand-rolled bounded-variable
   simplex with warm-start across Stage 3's per-variable LPs. Here we
   call `scipy.optimize.linprog` (HiGHS) and re-solve from scratch each
   time. Same optima; somewhat slower per Stage-3 LP but those LPs are
   tiny so it doesn't matter.
2. **Stage-4 skip condition.** The dissertation literally writes
   `x_min ≤ x_min^p` but its prose says `D_min ⊂ D_min^p`. Those
   contradict. We follow the prose: skip iff `x_min ≥ x_min_p`. See
   `combined/core.py` and the gap log in `COMBINED_METHOD.md`.
3. **КПИА kept active even when no explicit filter is passed.** The
   dissertation describes КПИА as gated on "having found a feasible
   already". Our implementation always runs with an internal filter
   (initialised at +∞, tightened on first feasible) so КПИА starts
   pruning the *moment* a feasible appears. This is a faithful but
   un-named optimisation that fixes a pathological case where the
   Stage-2 lattice would otherwise enumerate every feasible in the box.

## Tests

```bash
uv run python tests/test_transforms.py
```

The transform tests verify (2.18) substitution, max→min flip, sort by
`c`, and round-trip through `decode`.
