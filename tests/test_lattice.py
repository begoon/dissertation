"""Unit tests for combined.lattice (vector-lattice implicit enumeration).

Each test fixes a small, hand-traceable problem and asserts the search
produces the right point and (where load-bearing) that the right
pruning rule fired.

Convention: c is sorted ascending, c[j] ≥ 0, and rows are written in
the (2.17) "≤"-only form, matching what the orchestrator passes in.
"""

from __future__ import annotations

import numpy as np

from combined.lattice import Box, Filter, lattice_search

# ---------- finds the optimum on a known instance ----------


def test_finds_optimum_in_open_box() -> None:
    """min c·x s.t. -2x1 - 4x2 ≤ -11, 0 ≤ x ≤ 5. Optimum (0,3), z=15."""
    A = np.array([[-2.0, -4.0]])
    b = np.array([-11.0])
    c = np.array([3.0, 5.0])
    box = Box(low=np.zeros(2, dtype=np.int64), high=np.array([5, 5]))
    x, stats = lattice_search(A, b, c, box, use_kpp=True)
    assert x is not None
    assert (x == np.array([0, 3])).all()
    assert stats["feasibles"] >= 1


def test_returns_none_when_box_excludes_all_feasibles() -> None:
    """Same problem but box pinned at infeasible-only region."""
    A = np.array([[-2.0, -4.0]])
    b = np.array([-11.0])
    c = np.array([3.0, 5.0])
    # Force [0,0]..[1,2] — every point has 2x1+4x2 < 11 (infeasible).
    box = Box(low=np.zeros(2, dtype=np.int64), high=np.array([1, 2]))
    x, stats = lattice_search(A, b, c, box, use_kpp=True)
    assert x is None
    # КН should have fired since the box is fully infeasible
    assert stats["kh_prunes"] >= 1


# ---------- КН (eq. 2.9) pruning ----------


def test_kh_prunes_when_row_unrepairable() -> None:
    """Set up a single row that no remaining forward step can repair from 0."""
    # Row 0:  -1 x1 - 1 x2 ≤ -10  ⇔  x1 + x2 ≥ 10.
    # Box [0..3, 0..3]: max possible x1+x2 = 6 < 10, so КН fires at root.
    A = np.array([[-1.0, -1.0]])
    b = np.array([-10.0])
    c = np.array([1.0, 2.0])
    box = Box(low=np.zeros(2, dtype=np.int64), high=np.array([3, 3]))
    x, stats = lattice_search(A, b, c, box, use_kpp=True)
    assert x is None
    assert stats["kh_prunes"] == 1
    assert stats["nodes"] == 1  # we never descended past the root


# ---------- КПИА (eq. 2.12) pruning ----------


def test_kpia_arms_after_first_feasible() -> None:
    """With multiple feasibles, КПИА should prune the worse subtrees."""
    # Two-row LP: any x with x1 ≥ 1 OR x2 ≥ 1 within [0,3]² is "feasible".
    # Concretely:  -x1 - x2 ≤ -1.  Feasible at (1,0), (0,1), (1,1), …
    # Ascending c = (1, 2), so (1,0) is the unique optimum at z=1.
    A = np.array([[-1.0, -1.0]])
    b = np.array([-1.0])
    c = np.array([1.0, 2.0])
    box = Box(low=np.zeros(2, dtype=np.int64), high=np.array([3, 3]))
    x, stats = lattice_search(A, b, c, box, use_kpp=True)
    assert (x == np.array([1, 0])).all()
    # КПИА must have fired at least once after the first feasible
    assert stats["kpia_prunes"] >= 1


def test_strict_filter_prunes_equal_or_worse() -> None:
    """Strict filter c·x < f_best rejects solutions that match the bound."""
    # Same problem, but pass a Filter pre-pinned at f_best=1 — so an exactly-1
    # solution must NOT be returned. The search should find nothing.
    A = np.array([[-1.0, -1.0]])
    b = np.array([-1.0])
    c = np.array([1.0, 2.0])
    box = Box(low=np.zeros(2, dtype=np.int64), high=np.array([3, 3]))
    filt = Filter(c=c.copy(), f_best=1.0)
    x, _ = lattice_search(A, b, c, box, filt=filt, use_kpp=True)
    assert x is None


# ---------- КПП (eq. 2.13) ordering — affects perf, not optimum ----------


def test_kpp_off_returns_same_optimum() -> None:
    """Disabling КПП changes node counts but not the returned objective."""
    A = np.array([[-2.0, -4.0]])
    b = np.array([-11.0])
    c = np.array([3.0, 5.0])
    box = Box(low=np.zeros(2, dtype=np.int64), high=np.array([5, 5]))
    x_on, _ = lattice_search(A, b, c, box, use_kpp=True)
    x_off, _ = lattice_search(A, b, c, box, use_kpp=False)
    assert (x_on == x_off).all()
    assert (x_on == np.array([0, 3])).all()


# ---------- node limit ----------


def test_node_limit_cuts_search_short() -> None:
    """node_limit=N stops after at most N+1 nodes."""
    A = np.array([[-1.0, -1.0]])
    b = np.array([-1.0])
    c = np.array([1.0, 2.0])
    box = Box(low=np.zeros(2, dtype=np.int64), high=np.array([3, 3]))
    x, stats = lattice_search(A, b, c, box, use_kpp=True, node_limit=1)
    # Search visited at most one node beyond the limit check
    assert stats["nodes"] <= 2
    # x may or may not have been found depending on whether root was feasible


# ---------- slice mode (variable fixing) ----------


def test_slice_mode_searches_only_fixed_value() -> None:
    """Box with fixed_var=0, fixed_val=2 should never visit x[0] != 2."""
    A = np.array([[-1.0, -1.0]])
    b = np.array([-3.0])  # x1 + x2 ≥ 3
    c = np.array([1.0, 2.0])
    box = Box(
        low=np.array([2, 0]),
        high=np.array([2, 5]),
        fixed_var=0,
        fixed_val=2,
    )
    x, _ = lattice_search(A, b, c, box, use_kpp=True)
    assert x is not None
    assert x[0] == 2
    # x[1] minimal s.t. 2 + x[1] ≥ 3 ⇒ x[1] = 1
    assert x[1] == 1


# ---------- structural / edge cases ----------


def test_origin_feasible_returns_origin() -> None:
    """If the origin is already feasible, optimum is the origin (Rule 1)."""
    A = np.array([[1.0, 1.0]])  # x1 + x2 ≤ 5 — trivially satisfied at 0
    b = np.array([5.0])
    c = np.array([1.0, 2.0])
    box = Box(low=np.zeros(2, dtype=np.int64), high=np.array([3, 3]))
    x, stats = lattice_search(A, b, c, box, use_kpp=True)
    assert x is not None
    assert (x == np.zeros(2)).all()
    assert stats["feasibles"] == 1


def test_no_constraints_picks_box_low() -> None:
    """No row constraints (m=0) — every box point is feasible, optimum is low."""
    A = np.zeros((0, 2))
    b = np.zeros(0)
    c = np.array([1.0, 2.0])
    box = Box(low=np.array([1, 2]), high=np.array([3, 4]))
    x, _ = lattice_search(A, b, c, box, use_kpp=True)
    assert x is not None
    assert (x == np.array([1, 2])).all()
