"""Sanity-check baseline: solve the MILP with scipy.optimize.milp (HiGHS)."""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, Bounds, milp


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("task", type=Path)
    ap.add_argument("--time-limit", type=float, default=120.0)
    args = ap.parse_args()

    if args.task.suffix == ".json":
        from solve import _load_json

        p = _load_json(args.task)
    else:
        spec = importlib.util.spec_from_file_location("task_module", args.task)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"cannot import {args.task}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        p = module.build()

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
        options={"time_limit": args.time_limit},
    )
    wall = time.perf_counter() - t

    print(f"scipy.milp wall = {wall:.2f}s")
    print(f"status = {res.message}")
    if res.x is not None:
        x = np.round(res.x).astype(int)
        obj = float(p.c @ x) * (
            sign if p.direction == "max" else 1.0
        )  # report in user direction
        # easier: just compute c · x with original c and sign
        obj_user = (-1 if p.direction == "max" else 1) * float(sign * p.c @ x)
        print(f"objective (user sense): {float(p.c @ x):.6f}")
        nz = [(p.var_names[j], int(x[j])) for j in range(p.n) if int(x[j]) != 0]
        if nz:
            print("non-zero variables:")
            for name, val in nz:
                print(f"  {name} = {val}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
