"""Unit tests for combined.lp (the scipy.linprog wrapper)."""

from __future__ import annotations

import numpy as np

from combined.lp import solve_lp_min, solve_min_xj_with_filter


def test_min_with_le_constraint() -> None:
    """min x1+x2 s.t. x1+x2 ≤ 5, 0 ≤ x ≤ 5. Trivial: optimum at the origin."""
    res = solve_lp_min(
        c=np.array([1.0, 1.0]),
        A=np.array([[1.0, 1.0]]),
        b=np.array([5.0]),
        sense=["<="],
        h=np.array([5, 5]),
    )
    assert res.status == "optimal"
    assert abs(res.obj) < 1e-9
    assert np.allclose(res.x, [0.0, 0.0])


def test_min_with_ge_constraint() -> None:
    """min 3x1+5x2 s.t. 2x1+4x2 ≥ 11, 0 ≤ x ≤ 5. Optimum (0, 2.75), z=13.75."""
    res = solve_lp_min(
        c=np.array([3.0, 5.0]),
        A=np.array([[2.0, 4.0]]),
        b=np.array([11.0]),
        sense=[">="],
        h=np.array([5, 5]),
    )
    assert res.status == "optimal"
    assert abs(res.obj - 13.75) < 1e-7
    assert np.allclose(res.x, [0.0, 2.75])


def test_equality_constraint_is_respected() -> None:
    """min x1+x2 s.t. x1+x2 = 3, 0 ≤ x ≤ 5. Any point on the line scores 3."""
    res = solve_lp_min(
        c=np.array([1.0, 1.0]),
        A=np.array([[1.0, 1.0]]),
        b=np.array([3.0]),
        sense=["="],
        h=np.array([5, 5]),
    )
    assert res.status == "optimal"
    assert abs(res.obj - 3.0) < 1e-7
    assert abs(res.x.sum() - 3.0) < 1e-7


def test_infeasible_returns_infeasible() -> None:
    """x1+x2 ≥ 100 inside [0,1]^2 is infeasible — wrapper should detect it."""
    res = solve_lp_min(
        c=np.array([1.0, 1.0]),
        A=np.array([[1.0, 1.0]]),
        b=np.array([100.0]),
        sense=[">="],
        h=np.array([1, 1]),
    )
    assert res.status == "infeasible"
    assert res.x is None
    assert res.obj is None


def test_solve_min_xj_with_filter_picks_box_corner() -> None:
    """Per-variable LP used in Stage 3 of the Combined Method.

    For the standing example (min 3x1+5x2 s.t. 2x1+4x2 ≥ 11) with
    filter c·x ≤ 15, min x[1] is achieved at (2.5, 1.5) — both rows tight.
    """
    res = solve_min_xj_with_filter(
        c_filter=np.array([3.0, 5.0]),
        j=1,
        A=np.array([[-2.0, -4.0]]),  # ≥ 11 written as -≤-11
        b=np.array([-11.0]),
        sense=["<="],
        h=np.array([5, 5]),
        filter_rhs=15.0,
    )
    assert res.status == "optimal"
    assert abs(res.x[1] - 1.5) < 1e-7
    assert abs(res.x[0] - 2.5) < 1e-7
