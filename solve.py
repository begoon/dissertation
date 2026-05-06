"""CLI driver: load a task JSON or builder, run combined_method, print summary."""

from __future__ import annotations

import argparse
import importlib
import json
import sys
from pathlib import Path

import numpy as np

from combined import MILP, combined_method


def load_task(task_path: Path) -> MILP:
    """Load a task by .json (raw matrices) or .py (a builder module)."""
    if task_path.suffix == ".json":
        return _load_json(task_path)
    if task_path.suffix == ".py":
        return _load_py(task_path)
    raise ValueError(f"unsupported task file: {task_path}")


def _load_json(path: Path) -> MILP:
    data = json.loads(path.read_text())
    return MILP(
        c=np.array(data["c"], dtype=float),
        A=np.array(data["A"], dtype=float),
        b=np.array(data["b"], dtype=float),
        sense=list(data["sense"]),
        h=np.array(data["h"], dtype=np.int64),
        direction=data.get("direction", "min"),
        var_names=list(data.get("var_names", [])),
    )


def _load_py(path: Path) -> MILP:
    spec = importlib.util.spec_from_file_location("task_module", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "build"):
        raise RuntimeError(f"{path} must define a build() -> MILP function")
    return module.build()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("task", type=Path, help="path to .json or .py task")
    ap.add_argument("--eps", type=float, default=None,
                    help="relative tolerance %% for early-exit after Stage 2")
    ap.add_argument("--no-kpp", action="store_true", help="disable КПП ordering")
    ap.add_argument("--quiet", action="store_true", help="suppress per-stage log")
    ap.add_argument("--node-limit", type=int, default=None,
                    help="cap total lattice nodes per stage (returns best so far if hit)")
    ap.add_argument("--milp-heuristic", action="store_true",
                    help="seed Stage 2 with scipy.optimize.milp (fast feasible)")
    args = ap.parse_args(argv)

    p = load_task(args.task)
    print(f"task: {args.task}")
    print(f"  vars (n) = {p.n}, rows (m) = {p.m}, direction = {p.direction}")
    print()

    heuristic = None
    if args.milp_heuristic:
        from combined.heuristics import scipy_milp_heuristic
        heuristic = scipy_milp_heuristic

    res = combined_method(p, eps=args.eps, use_kpp=not args.no_kpp,
                          verbose=not args.quiet, node_limit=args.node_limit,
                          heuristic=heuristic)
    print()

    print(f"status:    {res.status}")
    if res.x is not None:
        print(f"objective: {res.objective:.6f}")
        if res.z_lp is not None:
            print(f"z*_LP:     {res.z_lp:.6f} (lower bound from Stage 1)")
        if res.z_sub is not None:
            print(f"z̃ (sub):   {res.z_sub:.6f} (Stage 2 result, in min-sense)")
        nz = [(p.var_names[j], int(res.x[j])) for j in range(p.n) if int(res.x[j]) != 0]
        if nz:
            print("non-zero variables:")
            for name, val in nz:
                print(f"  {name} = {val}")
    print()
    print("stage stats:")
    for name, st in res.stage_stats.items():
        print(f"  {name:14s}  nodes={st.nodes:>10d}  kh={st.kh_prunes:>8d}  "
              f"kpia={st.kpia_prunes:>8d}  feas={st.feasibles:>4d}  "
              f"t={st.wall_seconds:.3f}s")

    return 0 if res.status in ("optimal", "suboptimal_within_eps") else 2


if __name__ == "__main__":
    sys.exit(main())
