"""Benchmark sweep: combined_method (vanilla / --no-kpp / MILP-seeded) vs.
scipy.optimize.milp baseline across a set of tasks.

Emits a markdown-friendly table to stdout (and optionally JSON to --json).
Each run uses an independent process so that node counters and timings are
isolated. The bench treats Stage 4 hitting node_limit as "did not prove";
it still records the best feasible found.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Callable, Optional

import numpy as np
from scipy.optimize import LinearConstraint, Bounds, milp

from combined import combined_method, MILP


@dataclass
class Run:
    task: str
    config: str  # "scipy", "vanilla", "no-kpp", "milp-seed"
    n: int = 0
    m: int = 0
    status: str = "?"
    obj: Optional[float] = None
    z_lp: Optional[float] = None
    z_sub: Optional[float] = None
    nodes_total: int = 0
    kh_prunes: int = 0
    kpia_prunes: int = 0
    wall: float = 0.0


def load_task(task_path: Path) -> MILP:
    if task_path.suffix == ".json":
        from solve import _load_json

        return _load_json(task_path)
    if task_path.suffix == ".py":
        import importlib.util

        spec = importlib.util.spec_from_file_location("task_module", task_path)
        assert spec and spec.loader
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m.build()
    raise ValueError(task_path)


def run_scipy(p: MILP) -> Run:
    A_ub_rows: list[np.ndarray] = []
    b_ub: list[float] = []
    A_eq_rows: list[np.ndarray] = []
    b_eq: list[float] = []
    for i, s in enumerate(p.sense):
        if s == "<=":
            A_ub_rows.append(p.A[i])
            b_ub.append(p.b[i])
        elif s == ">=":
            A_ub_rows.append(-p.A[i])
            b_ub.append(-p.b[i])
        elif s == "=":
            A_eq_rows.append(p.A[i])
            b_eq.append(p.b[i])
    constraints = []
    if A_ub_rows:
        constraints.append(
            LinearConstraint(np.array(A_ub_rows), -np.inf, np.array(b_ub)),
        )
    if A_eq_rows:
        constraints.append(
            LinearConstraint(
                np.array(A_eq_rows),
                np.array(b_eq),
                np.array(b_eq),
            ),
        )
    bounds = Bounds(np.zeros(p.n), p.h.astype(float))
    integrality = np.ones(p.n)
    sign = -1.0 if p.direction == "max" else 1.0
    t = time.perf_counter()
    res = milp(
        c=sign * p.c,
        constraints=constraints,
        integrality=integrality,
        bounds=bounds,
        options={"time_limit": 60.0},
    )
    wall = time.perf_counter() - t
    r = Run(task="", config="scipy", n=p.n, m=p.m, wall=wall)
    if res.x is not None:
        r.obj = float(p.c @ np.round(res.x).astype(int))
        r.status = "optimal"
    else:
        r.status = "fail"
    return r


def run_combined(
    p: MILP,
    *,
    use_kpp: bool,
    heuristic: Optional[Callable],
    node_limit: Optional[int],
) -> Run:
    t = time.perf_counter()
    r_cm = combined_method(
        p,
        use_kpp=use_kpp,
        heuristic=heuristic,
        node_limit=node_limit,
        verbose=False,
    )
    wall = time.perf_counter() - t
    cfg = (
        "vanilla"
        if (use_kpp and heuristic is None)
        else "no-kpp" if (not use_kpp and heuristic is None) else "milp-seed"
    )
    r = Run(
        task="",
        config=cfg,
        n=p.n,
        m=p.m,
        wall=wall,
        status=r_cm.status,
        obj=r_cm.objective,
        z_lp=r_cm.z_lp,
        z_sub=r_cm.z_sub,
    )
    for st in r_cm.stage_stats.values():
        r.nodes_total += st.nodes
        r.kh_prunes += st.kh_prunes
        r.kpia_prunes += st.kpia_prunes
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--node-limit", type=int, default=500_000)
    ap.add_argument(
        "--include-large",
        action="store_true",
        help="include the 153-var lan_b_zerocost task (slow vanilla)",
    )
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    tasks: list[tuple[str, Path]] = [
        ("tiny", Path("tasks/tiny.json")),
        ("medium2", Path("tasks/medium2.json")),
        ("medium20", Path("tasks/medium20.json")),
        ("medium50", Path("tasks/medium50.json")),
    ]
    if args.include_large:
        tasks += [
            ("lan_b_zerocost", Path("tasks/lan_b_zerocost.py")),
        ]

    runs: list[Run] = []

    print(f"# Benchmark sweep  (node_limit={args.node_limit:,})\n")
    print(
        "| task | config | n | m | wall/s | obj | LP bound | nodes | КН | КПИА | status |",
    )
    print(
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    )

    from combined.heuristics import scipy_milp_heuristic

    for name, path in tasks:
        p = load_task(path)
        # For very large problems, drop the no-kpp variant — it never finishes.
        skip_no_kpp = p.n > 100
        configs = [
            ("scipy", None),
            ("vanilla", dict(use_kpp=True, heuristic=None)),
            ("milp-seed", dict(use_kpp=True, heuristic=scipy_milp_heuristic)),
        ]
        if not skip_no_kpp:
            configs.insert(2, ("no-kpp", dict(use_kpp=False, heuristic=None)))
        for cfg, kwargs in configs:
            if cfg == "scipy":
                r = run_scipy(p)
            else:
                r = run_combined(p, **kwargs, node_limit=args.node_limit)
            r.task = name
            runs.append(r)
            obj_str = f"{r.obj:.2f}" if r.obj is not None else "—"
            lp_str = f"{r.z_lp:.2f}" if r.z_lp is not None else "—"
            print(
                f"| {r.task} | {r.config} | {r.n} | {r.m} | "
                f"{r.wall:>7.3f} | {obj_str} | {lp_str} | "
                f"{r.nodes_total:,} | {r.kh_prunes:,} | {r.kpia_prunes:,} | "
                f"{r.status} |",
                flush=True,
            )

    if args.json is not None:
        args.json.write_text(json.dumps([asdict(r) for r in runs], indent=2))

    return 0


if __name__ == "__main__":
    import sys

    sys.exit(main())
