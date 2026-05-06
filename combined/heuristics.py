"""Stage-2 heuristics that plug into combined_method via the `heuristic=` arg.

A heuristic returns a feasible point (in (2.17) coordinates) for combined_method
to use as the seed for Stage 2's trivial-improvement loop, OR None to fall back
to the lattice search.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.optimize import LinearConstraint, Bounds, milp

from combined.transforms import Canonical


def scipy_milp_heuristic(
    x_opt_l: np.ndarray, can: Canonical, time_limit: float = 60.0
) -> Optional[np.ndarray]:
    """Use scipy.optimize.milp to produce an integer feasible quickly.

    This is a *heuristic* in the dissertation sense: it does not exploit the
    Combined Method's structure, but it gives Stage 2 a strong incumbent
    immediately so КПИА can prune Stage 4 hard.
    """
    n = can.c_p.shape[0]
    m = can.A_p.shape[0]
    constraints = [LinearConstraint(can.A_p, -np.inf, can.b_p)]
    bounds = Bounds(np.zeros(n), can.h_p.astype(float))
    integrality = np.ones(n)

    res = milp(
        c=can.c_p,
        constraints=constraints,
        integrality=integrality,
        bounds=bounds,
        options={"time_limit": time_limit},
    )
    if res.x is None:
        return None
    x = np.round(res.x).astype(np.int64)
    if (x < 0).any() or (x > can.h_p).any():
        return None
    if (can.A_p @ x > can.b_p + 1e-7).any():
        return None
    return x
