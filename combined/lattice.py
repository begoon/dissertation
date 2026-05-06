"""Vector-lattice implicit enumeration (§2.2.2 of the dissertation).

Search the integer lattice [box.low, box.high] for x ∈ Z^n satisfying
A x ≤ b, minimising c · x, optionally subject to a strict filter
c · x < f_best.

Pruning rules:

* КН   (2.9)  — infeasibility: a row that cannot be repaired even by
                taking every helping forward-step to its bound is fatal.
* КПИА (2.12) — plan-aware alternative elimination: if the filter slack
                at x cannot absorb the cost of stepping on j, prune that
                step.
* КПП  (2.13) — preferred-variable ordering: heuristic ordering over the
                surviving candidates; never prunes.

Pre: c[j] ≥ 0 ∀ j and c sorted ascending. The lex-traversal guard
J_x = {j : j ≥ j_x} (eq. 2.2/2.3) only visits each lattice point once
under that ordering.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class Box:
    """Inclusive integer box [low, high] with optional fixed coordinate."""

    low: np.ndarray
    high: np.ndarray
    fixed_var: int = -1                        # -1 ⇒ no fix
    fixed_val: int = 0                         # used only when fixed_var ≥ 0

    def __post_init__(self) -> None:
        self.low = np.asarray(self.low, dtype=np.int64)
        self.high = np.asarray(self.high, dtype=np.int64)
        assert self.low.shape == self.high.shape
        assert (self.low <= self.high).all(), (self.low, self.high)


@dataclass
class Filter:
    """Strict filter c · x < f_best  (eq. 2.10)."""

    c: np.ndarray
    f_best: float


def _initial_x(box: Box) -> np.ndarray:
    x = box.low.copy()
    if box.fixed_var >= 0:
        x[box.fixed_var] = box.fixed_val
    return x


def lattice_search(
    A: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    box: Box,
    filt: Optional[Filter] = None,
    use_kpp: bool = True,
    node_limit: Optional[int] = None,
    progress_every: Optional[int] = None,        # if set, print "nodes=N best=..." every N nodes
    progress_label: str = "lattice",
) -> tuple[Optional[np.ndarray], dict]:
    """Find argmin c·x over { x ∈ Z^n : A x ≤ b, box, c·x < filt.f_best }.

    Returns (x_best, stats). x_best is None if no feasible found within
    node_limit (or at all when node_limit is None).
    """
    n = box.low.shape[0]
    assert A.shape[1] == n and b.shape == (A.shape[0],) and c.shape == (n,)
    assert (c >= -1e-12).all(), "c must be non-negative; canonicalise first"

    x = _initial_x(box)
    y_main = b - A @ x                         # row slacks; feasible iff y_main >= 0

    # Always run with a filter, even if the caller didn't provide one. This
    # lets КПИА kick in as soon as Stage 2's first feasible is found.
    if filt is None:
        filt = Filter(c=c.copy(), f_best=float("inf"))
        strict = False                         # filter starts at +inf — non-strict during Stage 2
    else:
        strict = True                          # caller-supplied filter is strict (Stage 4)

    f_best = filt.f_best
    y_filt = f_best - float(c @ x)             # may be +inf at the start

    best_x: Optional[np.ndarray] = None
    stats = {"nodes": 0, "kh_prunes": 0, "kpia_prunes": 0, "feasibles": 0}

    EPS = 1e-9

    # j_x is the largest movable index that is currently above its lower bound,
    # or -1 if none. Maintained incrementally on every step.
    def initial_jx() -> int:
        for j in range(n - 1, -1, -1):
            if j == box.fixed_var:
                continue
            if x[j] > box.low[j]:
                return j
        return -1

    j_x = initial_jx()

    def visit() -> bool:
        """Examine the current node and recurse. Returns True if node_limit hit."""
        nonlocal j_x, f_best, y_filt, best_x, y_main

        stats["nodes"] += 1
        if progress_every is not None and stats["nodes"] % progress_every == 0:
            print(f"[{progress_label}] nodes={stats['nodes']:,} "
                  f"kh={stats['kh_prunes']:,} kpia={stats['kpia_prunes']:,} "
                  f"feas={stats['feasibles']} best_f={f_best:.6f} depth={int((x > box.low).sum())}",
                  flush=True)
        if node_limit is not None and stats["nodes"] > node_limit:
            return True

        # --- Rule 1: feasibility (then back up) ---
        if (y_main >= -EPS).all():
            f = float(c @ x)
            if f < f_best - EPS:
                f_best = f
                best_x = x.copy()
                filt.f_best = f_best
                stats["feasibles"] += 1
            return False

        # --- КН (2.9) ---
        I_x = np.where(y_main < -EPS)[0]
        if I_x.size:
            # max possible increase to y[i] from all remaining helping steps:
            #   max_inc[i] = Σ_{j helps row i, j still movable from here}  -A[i,j] · room[j]
            # j is movable iff j ≥ j_x, j != fixed_var, x[j] < box.high[j].
            j_lo = max(j_x, 0)
            movable_mask = np.zeros(n, dtype=bool)
            movable_mask[j_lo:] = True
            if box.fixed_var >= j_lo:
                movable_mask[box.fixed_var] = False
            movable_mask &= x < box.high
            room = (box.high - x).astype(float)
            room[~movable_mask] = 0.0
            # contribution to row i: -A[i,j] * room[j] for j with A[i,j] < 0
            A_Ix = A[I_x]                       # (|I_x|, n)
            neg_part = np.where(A_Ix < 0, A_Ix, 0.0)
            max_inc = -(neg_part @ room)        # length |I_x|
            new_y = y_main[I_x] + max_inc
            if (new_y < -EPS).any():
                stats["kh_prunes"] += 1
                return False

        # --- candidate set J_x with КПИА (2.12) filter ---
        j_lo = max(j_x, 0)
        cands: list[int] = []
        for j in range(j_lo, n):
            if j == box.fixed_var:
                continue
            if x[j] >= box.high[j]:
                continue
            if c[j] >= y_filt - 1e-12:
                stats["kpia_prunes"] += 1
                continue
            cands.append(j)

        if not cands:
            return False

        # --- КПП (2.13) ordering ---
        if use_kpp and len(cands) > 1 and I_x.size:
            A_Ix = A[I_x]                       # (|I_x|, n)
            col_sum = A_Ix.sum(axis=0)          # (n,)
            scores = [(box.high[j] - x[j]) * col_sum[j] for j in cands]
            cands = [j for _, j in sorted(zip(scores, cands))]

        # --- recurse ---
        for j in cands:
            saved_jx = j_x
            saved_yfilt = y_filt
            x[j] += 1
            y_main -= A[:, j]
            y_filt -= c[j]
            j_x = max(saved_jx, j)
            hit_limit = visit()
            x[j] -= 1
            y_main += A[:, j]
            # the child may have tightened f_best — recompute y_filt from the
            # (possibly new) f_best rather than just undoing the −c[j].
            y_filt = filt.f_best - float(c @ x)
            j_x = saved_jx
            # Keep `saved_yfilt` linter-quiet: we deliberately recompute instead.
            del saved_yfilt
            if hit_limit:
                return True
        return False

    visit()
    return best_x, stats
