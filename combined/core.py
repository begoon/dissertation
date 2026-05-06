"""Combined Method orchestrator (§2.3 of COMBINED_METHOD.md).

Stages 1–5: LP relaxation → sub-optimal feasible via lattice search →
filter row + per-variable corner LPs → strict-filter lattice search → decode.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from combined.lattice import Box, Filter, lattice_search
from combined.lp import LPResult, solve_lp_min, solve_min_xj_with_filter
from combined.model import MILP
from combined.transforms import Canonical, canonicalize, decode

EPS = 1e-9
LATTICE_NODE_LIMIT = 10_000_000                 # safety net


@dataclass
class StageStats:
    nodes: int = 0
    kh_prunes: int = 0
    kpia_prunes: int = 0
    feasibles: int = 0
    wall_seconds: float = 0.0


@dataclass
class Result:
    status: str                                 # "optimal", "suboptimal_within_eps", "infeasible", "node_limit"
    x: Optional[np.ndarray] = None              # solution in user variable order, or None
    objective: Optional[float] = None           # in user min/max sense
    z_sub: Optional[float] = None               # Stage-2 sub-optimum objective (in min/(2.17) sense)
    z_lp: Optional[float] = None                # Stage-1 LP optimum (in min/(2.17) sense)
    stage_stats: dict[str, StageStats] = field(default_factory=dict)
    log: list[str] = field(default_factory=list)


# A heuristic returning either a feasible (in (2.17) coordinates) or None.
HeuristicFn = Callable[[np.ndarray, "Canonical"], Optional[np.ndarray]]


def combined_method(
    p: MILP,
    eps: Optional[float] = None,                # relative tolerance in % (eq. 2.27); if hit, stop after Stage 2
    heuristic: Optional[HeuristicFn] = None,
    use_kpp: bool = True,
    verbose: bool = False,
    node_limit: Optional[int] = None,           # per-lattice-call cap; default LATTICE_NODE_LIMIT
) -> Result:
    """Solve the user MILP via the Combined Method."""
    res = Result(status="initialising")

    def log_line(msg: str) -> None:
        res.log.append(msg)
        if verbose:
            print(msg, flush=True)
    log = res.log  # kept for back-compat; new code calls log_line()

    # ----- §3   Build canonical forms (2.17) and (2.19) -----
    can = canonicalize(p)
    n = can.n_struct
    log_line(f"[setup] n={n}, m_p={can.A_p.shape[0]}, m_l={can.A_l.shape[0]}, "
               f"flipped={int(can.flipped_vars.sum())}")

    # ----- Stage 1 ----- LP relaxation of (2.19), in (2.17) sign
    t = time.perf_counter()
    lp = solve_lp_min(can.c_p, can.A_p_for_lp(), can.b_p_for_lp(),
                      can.sense_for_lp(), can.h_p)
    s1 = StageStats()
    s1.wall_seconds = time.perf_counter() - t
    res.stage_stats["stage1_lp"] = s1

    if lp.status == "infeasible":
        log_line("[stage1] LP infeasible ⇒ MILP infeasible")
        res.status = "infeasible"
        return res
    if lp.status == "unbounded":
        log_line("[stage1] LP unbounded — should not happen with bounded h")
        res.status = "unbounded"
        return res
    if lp.status != "optimal":
        log_line(f"[stage1] LP error: {lp.status}")
        res.status = "error"
        return res

    z_lp = float(lp.obj)
    res.z_lp = z_lp
    x_opt_l = lp.x.copy()
    log_line(f"[stage1] LP optimum z*_LP = {z_lp:.6f}")

    # If LP solution happens to be integer and feasible for (2.17), we are done.
    x_round = np.round(x_opt_l).astype(np.int64)
    integral = np.allclose(x_opt_l, x_round, atol=1e-7)
    if integral and (can.A_p @ x_round <= can.b_p + 1e-7).all() and \
       (x_round >= 0).all() and (x_round <= can.h_p).all():
        log_line("[stage1] LP optimum is integer feasible — done.")
        res.status = "optimal"
        res.x = decode(x_round, can)
        sign = -1.0 if can.flipped_max_to_min else 1.0
        # Objective in user form: undo the dropped constant + sign flip.
        res.objective = sign * (float(p.c @ res.x))
        return res

    # x_min^p (eq. 2.20): floor of LP optimum, clipped to bounds
    x_min_p = np.floor(x_opt_l).astype(np.int64)
    x_min_p = np.clip(x_min_p, 0, can.h_p)
    log_line(f"[stage1] x_min_p = {x_min_p.tolist()}")

    # ----- Stage 2 ----- search for sub-optimal feasible
    t = time.perf_counter()
    box = Box(low=x_min_p.copy(), high=can.h_p.copy())
    s2 = StageStats()

    x_D_p: Optional[np.ndarray] = None

    # 2.0 — optional heuristic shortcut
    if heuristic is not None:
        try:
            x_h = heuristic(x_opt_l, can)
            if x_h is not None and _is_feasible(x_h, can):
                x_D_p = x_h.astype(np.int64)
                log_line(f"[stage2] heuristic produced feasible, c·x = {can.c_p @ x_D_p:.6f}")
        except Exception as e:                  # don't let heuristic crashes kill the run
            log_line(f"[stage2] heuristic raised {e!r}; falling back")

    # 2.1.1 — initial enumeration in [x_min^p, h^p]
    if x_D_p is None:
        x_D_p, st = lattice_search(can.A_p, can.b_p, can.c_p, box,
                                   filt=None, use_kpp=use_kpp,
                                   node_limit=node_limit if node_limit is not None else LATTICE_NODE_LIMIT,
                                   progress_every=100_000 if verbose else None,
                                   progress_label="stage2.1.1")
        _accumulate(s2, st)
        log_line(f"[stage2.1.1] lattice nodes={st['nodes']}, "
                   f"kh={st['kh_prunes']}, found={x_D_p is not None}")

    # 2.1.2 — slice expansion if no feasible found yet
    excluded: set[int] = set()
    while x_D_p is None:
        # Pick the smallest-index movable var with box.low > 0, not excluded
        k = -1
        for j in range(n):
            if j in excluded:
                continue
            if box.low[j] > 0:
                k = j
                break
        if k == -1:
            log_line("[stage2.1.2] no expansion possible ⇒ MILP infeasible")
            res.status = "infeasible"
            res.stage_stats["stage2"] = s2
            s2.wall_seconds = time.perf_counter() - t
            return res

        # Try the slice with x[k] = box.low[k] - 1
        new_v = int(box.low[k] - 1)
        slice_box = Box(low=box.low.copy(), high=box.high.copy(),
                        fixed_var=k, fixed_val=new_v)
        slice_box.low[k] = new_v               # so the lex/КН logic sees x[k] >= new_v
        slice_box.high[k] = new_v              # and ≤ new_v — i.e. exactly new_v

        x_slice, st = lattice_search(can.A_p, can.b_p, can.c_p, slice_box,
                                     filt=None, use_kpp=use_kpp,
                                     node_limit=node_limit if node_limit is not None else LATTICE_NODE_LIMIT,
                                     progress_every=100_000 if verbose else None,
                                     progress_label=f"stage2.1.2 j={k}")
        _accumulate(s2, st)
        log_line(f"[stage2.1.2] expand on j={k} to {new_v}: "
                   f"nodes={st['nodes']}, found={x_slice is not None}")

        # If КН fired at the very first node, this variable cannot be expanded further.
        if st["nodes"] <= 1 and x_slice is None:
            excluded.add(k)
            continue

        if x_slice is not None:
            x_D_p = x_slice
        # Update the cumulative search box: low[k] := new_v
        box.low[k] = new_v

    z_D_p = float(can.c_p @ x_D_p)
    log_line(f"[stage2.1] found feasible x_D_p, z = {z_D_p:.6f}")

    # 2.2 — trivial improvement (loop until no further single-var decrement helps)
    x_tilde = x_D_p.copy()
    z_tilde = z_D_p
    improved = True
    while improved:
        improved = False
        best_j = -1
        best_delta = 0
        best_gain = 0.0
        for j in range(n):
            # how much can x_tilde[j] decrease without leaving D^p?
            if x_tilde[j] == 0:
                continue
            if can.c_p[j] <= 0:
                continue
            # Try the largest decrement that keeps A x ≤ b. The delta on row i
            # from reducing x[j] by d is: y[i] += d * A[i, j] (since slack increases
            # when A[i, j] > 0). For d > 0, slacks can only stay non-negative if
            # for every i with A[i, j] < 0 we have y[i] - d * |A[i, j]| ≥ 0.
            y = can.b_p - can.A_p @ x_tilde
            d_max = int(x_tilde[j])
            for i in range(can.A_p.shape[0]):
                a = can.A_p[i, j]
                if a < 0:
                    d_cap = int(math.floor(y[i] / -a))  # y[i] + d*a ≥ 0 ⇒ d ≤ y[i]/-a
                    if d_cap < d_max:
                        d_max = d_cap
            if d_max <= 0:
                continue
            gain = can.c_p[j] * d_max
            if gain > best_gain + EPS:
                best_gain, best_j, best_delta = gain, j, d_max
        if best_j >= 0:
            x_tilde[best_j] -= best_delta
            z_tilde -= best_gain
            improved = True
    log_line(f"[stage2.2] after trivial improvement: z̃ = {z_tilde:.6f} "
               f"(improvement = {z_D_p - z_tilde:.6f})")

    s2.wall_seconds = time.perf_counter() - t
    res.stage_stats["stage2"] = s2
    res.z_sub = z_tilde

    # 2.3 — relative-tolerance early exit
    if eps is not None:
        denom = abs(z_tilde) if abs(z_tilde) > 1e-12 else 1.0
        gap_pct = (z_tilde - z_lp) / denom * 100.0
        log_line(f"[stage2.3] gap upper bound = {gap_pct:.4f}% (target {eps}%)")
        if gap_pct <= eps:
            log_line("[stage2.3] within tolerance — accepting sub-optimum")
            res.status = "suboptimal_within_eps"
            res.x = decode(x_tilde, can)
            sign = -1.0 if can.flipped_max_to_min else 1.0
            res.objective = sign * float(p.c @ res.x)
            return res

    # ----- Stage 3 ----- compute integer corner x_min for Stage-4 box
    t = time.perf_counter()
    s3 = StageStats()
    x_min = np.zeros(n, dtype=np.int64)

    # Candidate set J = vars not at LP lower bound (≈ 0). For the rest, x_min[j] = 0.
    candidates = [j for j in range(n) if x_opt_l[j] > 1e-7]
    log_line(f"[stage3] solving {len(candidates)} per-variable LPs with filter rhs = {z_tilde:.6f}")

    sense_p = ["<="] * can.A_p.shape[0]
    for j in candidates:
        sub = solve_min_xj_with_filter(can.c_p, j, can.A_p, can.b_p, sense_p,
                                       can.h_p, filter_rhs=z_tilde)
        if sub.status != "optimal":
            log_line(f"[stage3] per-var LP for j={j} status={sub.status}; assuming x_min[j]=0")
            continue
        x_min[j] = int(math.ceil(sub.x[j] - 1e-7))
    log_line(f"[stage3] x_min = {x_min.tolist()}")
    s3.wall_seconds = time.perf_counter() - t
    res.stage_stats["stage3"] = s3

    # ----- Stage 4 ----- final lattice search inside [x_min, h^p] with strict filter
    t = time.perf_counter()
    s4 = StageStats()

    # Short-circuit: if x_min ≥ x_min_p componentwise, the Stage-4 box [x_min, h^p]
    # is contained in the (already exhausted) Stage-2 box [x_min_p, h^p]. Stage 2
    # found the minimum across that box, so no strict improvement can exist.
    # NB: the dissertation literally writes "x_min ≤ x_min^p" but its own prose
    # says "D_min ⊂ D_min^p" — the two contradict each other and the ⊂ reading is
    # the algorithmically-correct one. We follow the prose.
    if (x_min >= x_min_p).all():
        log_line("[stage4] x_min ≥ x_min_p ⇒ Stage-4 box ⊂ Stage-2 box; sub-opt is optimum")
        x_opt_p = x_tilde
    else:
        if (x_min > can.h_p).any():
            log_line("[stage4] x_min exceeds h^p ⇒ box empty; sub-opt is optimum")
            x_opt_p = x_tilde
        else:
            box4 = Box(low=x_min.copy(), high=can.h_p.copy())
            filt = Filter(c=can.c_p.copy(), f_best=z_tilde)
            x_better, st = lattice_search(can.A_p, can.b_p, can.c_p, box4,
                                          filt=filt, use_kpp=use_kpp,
                                          node_limit=node_limit if node_limit is not None else LATTICE_NODE_LIMIT,
                                          progress_every=100_000 if verbose else None,
                                          progress_label="stage4")
            _accumulate(s4, st)
            log_line(f"[stage4] lattice nodes={st['nodes']}, kh={st['kh_prunes']}, "
                       f"kpia={st['kpia_prunes']}, found_better={x_better is not None}")
            if x_better is not None:
                x_opt_p = x_better
            else:
                x_opt_p = x_tilde
    s4.wall_seconds = time.perf_counter() - t
    res.stage_stats["stage4"] = s4

    # ----- Stage 5 ----- decode
    z_opt_p = float(can.c_p @ x_opt_p)
    log_line(f"[stage5] optimum z = {z_opt_p:.6f}")
    res.status = "optimal"
    res.x = decode(x_opt_p, can)
    sign = -1.0 if can.flipped_max_to_min else 1.0
    res.objective = sign * float(p.c @ res.x)

    if verbose:
        for line in log:
            print(line)

    return res


# ---------- helpers ----------

def _accumulate(stats: StageStats, st: dict) -> None:
    stats.nodes += st["nodes"]
    stats.kh_prunes += st["kh_prunes"]
    stats.kpia_prunes += st["kpia_prunes"]
    stats.feasibles += st["feasibles"]


def _is_feasible(x: np.ndarray, can: Canonical) -> bool:
    """Check x ∈ Z^n_+, x ≤ h^p, A^p x ≤ b^p."""
    if x.shape[0] != can.n_struct:
        return False
    if (x < 0).any() or (x > can.h_p).any():
        return False
    if (can.A_p @ x > can.b_p + 1e-7).any():
        return False
    return True


# ---------- patch Canonical with LP-form helpers ----------
# These methods adapt the (2.17) form to scipy's solve_lp_min, which takes
# (c, A, b, sense, h) directly. We feed it the (2.17) shape (after canonicalize
# all rows are ≤), so sense is just ["<="] * m_p.

def _A_p_for_lp(self: Canonical) -> np.ndarray:
    return self.A_p

def _b_p_for_lp(self: Canonical) -> np.ndarray:
    return self.b_p

def _sense_for_lp(self: Canonical) -> list[str]:
    return ["<="] * self.A_p.shape[0]

Canonical.A_p_for_lp = _A_p_for_lp        # type: ignore[attr-defined]
Canonical.b_p_for_lp = _b_p_for_lp        # type: ignore[attr-defined]
Canonical.sense_for_lp = _sense_for_lp    # type: ignore[attr-defined]
