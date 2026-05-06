"""LP layer used by Stages 1 and 3 — wraps scipy.optimize.linprog (HiGHS).

The dissertation specifies a hand-rolled bounded-variable simplex that
keeps a warm-start tableau across Stage 3's per-variable LPs. Here we use
HiGHS via scipy and re-solve from scratch each time — functionally
equivalent (same optima), at a per-LP cost that's tiny for the problem
sizes in the dissertation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.optimize import linprog


@dataclass
class LPResult:
    status: str  # "optimal", "infeasible", "unbounded"
    x: Optional[np.ndarray] = None  # shape (n_struct,) — structural vars only
    obj: Optional[float] = None  # value of c · x in the *minimisation* sense
    raw_x: Optional[np.ndarray] = (
        None  # full vector incl. slacks (or whatever scipy returned)
    )


def solve_lp_min(
    c: np.ndarray,  # length n_struct  (no slacks)
    A: np.ndarray,  # shape (m, n_struct)  — main constraints, no slack columns
    b: np.ndarray,
    sense: list[str],
    h: np.ndarray,  # length n_struct, integer upper bounds
) -> LPResult:
    """Minimise c · x subject to A x [<=,>=,=] b, 0 ≤ x ≤ h. Continuous."""
    n = c.shape[0]
    m = len(sense)

    A_ub_rows: list[np.ndarray] = []
    b_ub: list[float] = []
    A_eq_rows: list[np.ndarray] = []
    b_eq: list[float] = []
    for i, s in enumerate(sense):
        if s == "<=":
            A_ub_rows.append(A[i])
            b_ub.append(b[i])
        elif s == ">=":
            A_ub_rows.append(-A[i])
            b_ub.append(-b[i])
        elif s == "=":
            A_eq_rows.append(A[i])
            b_eq.append(b[i])
        else:
            raise ValueError(f"unknown sense {s!r}")

    A_ub = np.array(A_ub_rows) if A_ub_rows else None
    A_eq = np.array(A_eq_rows) if A_eq_rows else None
    bounds = [(0.0, float(h[j])) for j in range(n)]

    res = linprog(
        c=c,
        A_ub=A_ub,
        b_ub=np.array(b_ub) if b_ub else None,
        A_eq=A_eq,
        b_eq=np.array(b_eq) if b_eq else None,
        bounds=bounds,
        method="highs",
    )
    if res.status == 0:
        return LPResult(
            "optimal",
            x=np.asarray(res.x),
            obj=float(res.fun),
            raw_x=np.asarray(res.x),
        )
    if res.status == 2:
        return LPResult("infeasible")
    if res.status == 3:
        return LPResult("unbounded")
    return LPResult(f"linprog_status_{res.status}: {res.message}")


def solve_min_xj_with_filter(
    c_filter: np.ndarray,  # length n_struct  — c^p (used only inside the filter row)
    j: int,
    A: np.ndarray,
    b: np.ndarray,
    sense: list[str],
    h: np.ndarray,
    filter_rhs: float,  # z̃_D^p — the filter row's RHS
) -> LPResult:
    """Solve  min x[j]  s.t. A x [sense] b  AND  c^p · x ≤ filter_rhs  AND  0 ≤ x ≤ h.

    Used per-variable in Stage 3 to compute the integer corner x_min via ceil.
    """
    n = c_filter.shape[0]
    obj = np.zeros(n)
    obj[j] = 1.0  # minimise x[j]

    # Append the filter row as a <= constraint
    A2 = np.vstack([A, c_filter[None, :]])
    b2 = np.concatenate([b, [filter_rhs]])
    sense2 = sense + ["<="]

    return solve_lp_min(obj, A2, b2, sense2, h)
