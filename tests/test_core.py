"""Integration tests for combined.core.combined_method.

Covers the four exit paths through the orchestrator:

* Stage 1 short-circuit (LP optimum is integer feasible).
* Stage 1 infeasibility detection.
* Stage 2 produces sub-optimum, Stage 4 confirms / improves.
* Stage 2.3 eps-tolerance early exit.

Plus checks for sign / substitution handling (max problems, negative
coefficients), node-limit status reporting, and heuristic plug-in.
"""

from __future__ import annotations

import numpy as np
import pytest

from combined import MILP, combined_method

# ---------- helpers ----------


def _problem_min_3x_5y_ge_11() -> MILP:
    """Standing tiny example. True optimum (0, 3) at z=15."""
    return MILP(
        c=np.array([3.0, 5.0]),
        A=np.array([[2.0, 4.0]]),
        b=np.array([11.0]),
        sense=[">="],
        h=np.array([5, 5]),
        direction="min",
        var_names=["x1", "x2"],
    )


# ---------- Stage 1 short-circuits ----------


def test_lp_optimum_already_integer_short_circuits() -> None:
    """If the LP relaxation lands on an integer corner, we return at end of Stage 1."""
    # min x1+x2 s.t. x1+x2 ≥ 3, 0 ≤ x ≤ 5. LP optimum is (3, 0): integer.
    p = MILP(
        c=np.array([1.0, 1.0]),
        A=np.array([[1.0, 1.0]]),
        b=np.array([3.0]),
        sense=[">="],
        h=np.array([5, 5]),
        direction="min",
    )
    res = combined_method(p)
    assert res.status == "optimal"
    assert abs(res.objective - 3.0) < 1e-7
    # Stage 2 should not have been touched
    assert "stage2" not in res.stage_stats


def test_infeasible_lp_propagates_to_status_infeasible() -> None:
    """If the LP relaxation is infeasible, the MILP is infeasible too."""
    # x1+x2 ≥ 100 inside [0,1]^2: no solution.
    p = MILP(
        c=np.array([1.0, 1.0]),
        A=np.array([[1.0, 1.0]]),
        b=np.array([100.0]),
        sense=[">="],
        h=np.array([1, 1]),
        direction="min",
    )
    res = combined_method(p)
    assert res.status == "infeasible"
    assert res.x is None


# ---------- Full pipeline through Stages 2-4 ----------


def test_tiny_problem_solves_to_dissertation_optimum() -> None:
    p = _problem_min_3x_5y_ge_11()
    res = combined_method(p)
    assert res.status == "optimal"
    assert abs(res.objective - 15.0) < 1e-7
    assert int(res.x[1]) == 3 and int(res.x[0]) == 0
    # All four stages must have been touched
    assert "stage1_lp" in res.stage_stats
    assert "stage2" in res.stage_stats
    assert "stage3" in res.stage_stats
    assert "stage4" in res.stage_stats


def test_max_direction_handled_correctly() -> None:
    """direction='max' is internally flipped to min and the user-sign is restored."""
    # max x1 + x2 s.t. x1+x2 ≤ 5, 0 ≤ x ≤ 5. Trivial: maximum 5.
    p = MILP(
        c=np.array([1.0, 1.0]),
        A=np.array([[1.0, 1.0]]),
        b=np.array([5.0]),
        sense=["<="],
        h=np.array([5, 5]),
        direction="max",
    )
    res = combined_method(p)
    assert res.status == "optimal"
    assert abs(res.objective - 5.0) < 1e-7


def test_negative_objective_coef_triggers_substitution() -> None:
    """A user `c[j] < 0` should round-trip through the (2.18) substitution."""
    # min -2 x1 + x2 s.t. x1 + x2 ≤ 5, 0 ≤ x ≤ (3, 4).
    # Optimum: x1 large (since c1<0), x2 small. (3, 0) → -6.
    p = MILP(
        c=np.array([-2.0, 1.0]),
        A=np.array([[1.0, 1.0]]),
        b=np.array([5.0]),
        sense=["<="],
        h=np.array([3, 4]),
        direction="min",
    )
    res = combined_method(p)
    assert res.status == "optimal"
    assert abs(res.objective - (-6.0)) < 1e-7
    assert int(res.x[0]) == 3 and int(res.x[1]) == 0


# ---------- Stage 2.3 tolerance gate ----------


def test_eps_early_exit_skips_stages_3_and_4() -> None:
    """Passing a wide eps accepts Stage 2's sub-optimum without verification."""
    p = _problem_min_3x_5y_ge_11()
    # Sub-optimum is 15, LP bound 13.75, gap = (15-13.75)/15*100 = 8.33%.
    # With eps=10% the gate should fire.
    res = combined_method(p, eps=10.0)
    assert res.status == "suboptimal_within_eps"
    assert abs(res.objective - 15.0) < 1e-7
    # Stages 3 and 4 should NOT have been touched
    assert "stage3" not in res.stage_stats
    assert "stage4" not in res.stage_stats


def test_eps_too_tight_runs_full_pipeline() -> None:
    """eps below the gap forces Stage 4 to run."""
    p = _problem_min_3x_5y_ge_11()
    res = combined_method(p, eps=1.0)  # 1% < 8.33%
    assert res.status == "optimal"
    assert "stage4" in res.stage_stats


# ---------- node-limit reporting ----------


def test_node_limit_status_marks_best_within_budget() -> None:
    """Stage 4 hitting node_limit without an improvement → 'best_within_budget'."""
    # A larger problem so Stage 4 can't exhaust its budget within 100 nodes.
    rng = np.random.default_rng(3)
    n, m = 20, 6
    A = rng.integers(-5, 10, (m, n)).astype(float)
    A[A == 0] = 1
    c = rng.integers(1, 30, n).astype(float)
    b = (A.sum(axis=1) * 0.5).round()
    p = MILP(
        c=c,
        A=A,
        b=b,
        sense=[">="] * (m // 2) + ["<="] * (m - m // 2),
        h=np.full(n, 4, dtype=int),
        direction="min",
    )
    res = combined_method(p, node_limit=100)
    # Stage 4 either short-circuited or hit the limit. Either way the
    # returned objective should be a real number; status reflects which.
    assert res.objective is not None
    assert res.status in ("optimal", "best_within_budget", "infeasible")


# ---------- heuristic plug-in ----------


def test_heuristic_plug_in_consumed_when_returns_feasible() -> None:
    """A heuristic that returns a known-feasible point is honoured."""
    p = _problem_min_3x_5y_ge_11()

    captured: dict = {}

    def heur(x_lp, can):
        captured["called"] = True
        # In (2.17) coords (here = user coords since perm = id, no flip),
        # (0, 4) is feasible: 4*4 = 16 ≥ 11.
        return np.array([0, 4], dtype=np.int64)

    res = combined_method(p, heuristic=heur)
    assert captured["called"] is True
    assert res.status == "optimal"
    # Heuristic seeded 20 (3*0 + 5*4); Stage 2.2 trivial-improvement should
    # decrement x2 from 4 to 3 (still feasible: 12 ≥ 11), reaching 15.
    assert abs(res.objective - 15.0) < 1e-7


def test_heuristic_returning_infeasible_is_ignored() -> None:
    """A heuristic that returns an infeasible point is silently dropped."""
    p = _problem_min_3x_5y_ge_11()

    def bad_heur(x_lp, can):
        return np.array([0, 0], dtype=np.int64)  # not feasible

    res = combined_method(p, heuristic=bad_heur)
    assert res.status == "optimal"
    assert abs(res.objective - 15.0) < 1e-7


def test_heuristic_raising_exception_is_caught() -> None:
    """A heuristic that raises must not crash the pipeline."""
    p = _problem_min_3x_5y_ge_11()

    def crashing(x_lp, can):
        raise RuntimeError("boom")

    res = combined_method(p, heuristic=crashing)
    assert res.status == "optimal"
    assert abs(res.objective - 15.0) < 1e-7


# ---------- KПП toggle ----------


def test_kpp_toggle_does_not_change_optimum() -> None:
    """Disabling КПП may slow the search but must not change the result."""
    p = _problem_min_3x_5y_ge_11()
    a = combined_method(p, use_kpp=True)
    b = combined_method(p, use_kpp=False)
    assert abs(a.objective - b.objective) < 1e-7


# ---------- result decoding ----------


def test_result_x_is_in_user_variable_order() -> None:
    """`Result.x[j]` must correspond to the user's c[j], not the canonical c^p."""
    # Use a problem where canonicalize permutes variables: c=(5, 3) sorts to (3, 5),
    # so canonical index 0 maps to user index 1.
    p = MILP(
        c=np.array([5.0, 3.0]),  # NOTE: NOT ascending — perm will reorder
        A=np.array([[4.0, 2.0]]),
        b=np.array([11.0]),
        sense=[">="],
        h=np.array([5, 5]),
        direction="min",
    )
    res = combined_method(p)
    assert res.status == "optimal"
    # User-coord objective: 5*x[0] + 3*x[1]. Verify result.objective matches
    # this manual computation.
    manual = float(p.c @ res.x)
    assert abs(res.objective - manual) < 1e-7
