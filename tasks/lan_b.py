"""Build the LAN-Variant-B MILP (§1.1.2, eqs. 1.7–1.11) from the JSON data.

Decision variables, in this implementation order:

    x[d]  for d ∈ D       — user-flow on arc d (0..KPS)
    z[d]  for d ∈ D       — binary: cable laid on arc d?
    y[v]  for v ∈ V       — hub stack height at vertex v (0..h_y[v])

The arc set D is built from cable_costs: an arc (i, j) exists iff
i != j AND cable_costs[i][j] is finite (we treat 0 off-diagonal as ∞ per
the dissertation convention).

Reductions per §B.4 of LAN_TASK.md:

* drop y[v] for vertices with no users (v != S) — though we keep y[S]
  because the central commutator counts as a vertex with negative-pseudo
  user count;
* replace flow equality (1.8) by a ≥ inequality (the optimum is the
  same since over-supply costs cable but earns no benefit).

Returns a `combined.MILP` ready to feed to `combined_method`.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from combined import MILP


def build(zero_means_infinity: bool = True) -> MILP:
    """Build the LAN-B model.

    `zero_means_infinity` controls how off-diagonal `0`s in the cost matrix
    are interpreted. The dissertation's prose says zeros mean "no cable
    possible" (∞), but its variable count (153) implies it kept all
    directed pairs in the model, treating 0s as zero cost. Pass False to
    reproduce the dissertation's variable count and (per §4.2.2) optimum
    280.
    """
    data = json.loads((Path(__file__).parent / "lan_b.json").read_text())
    return _build_from_dict(data, zero_means_infinity=zero_means_infinity)


def _build_from_dict(data: dict, zero_means_infinity: bool = True) -> MILP:
    m_rooms = int(data["rooms"])
    S = int(data["central_commutator"]) - 1  # 1-based → 0-based
    KP_orig = np.array(data["users_per_room"], dtype=int)  # KP[v] for v=0..m-1
    c_k = float(data["hub_cost"])
    a = int(data["ports_per_hub"])
    cost = np.array(data["cable_costs"], dtype=float)

    # Apply the dissertation's KP[S] = -KPS trick.
    KPS = int(KP_orig.sum())
    KP = KP_orig.copy()
    KP[S] = -KPS  # ЦК absorbs all traffic

    # ---- build arc set D -----------------------------------------------------
    arcs: list[tuple[int, int]] = []
    for u in range(m_rooms):
        for v in range(m_rooms):
            if u == v:
                continue
            if zero_means_infinity and cost[u, v] == 0:
                continue  # treat 0 off-diag as no cable
            arcs.append((u, v))
    n_arcs = len(arcs)
    arc_index = {arc: idx for idx, arc in enumerate(arcs)}

    # in / out arc lists per vertex
    D_in: list[list[int]] = [[] for _ in range(m_rooms)]  # arcs INTO v
    D_out: list[list[int]] = [[] for _ in range(m_rooms)]  # arcs OUT OF v
    for d_idx, (u, v) in enumerate(arcs):
        D_out[u].append(d_idx)
        D_in[v].append(d_idx)

    # ---- variable layout -----------------------------------------------------
    # block 1: x[d] for d ∈ D                  indices [0, n_arcs)
    # block 2: z[d] for d ∈ D                  indices [n_arcs, 2 n_arcs)
    # block 3: y[v] for v ∈ V (all vertices)   indices [2 n_arcs, 2 n_arcs + m_rooms)
    n = 2 * n_arcs + m_rooms

    def x_idx(d_idx: int) -> int:
        return d_idx

    def z_idx(d_idx: int) -> int:
        return n_arcs + d_idx

    def y_idx(v: int) -> int:
        return 2 * n_arcs + v

    var_names = (
        [f"x[{u},{v}]" for (u, v) in arcs]
        + [f"z[{u},{v}]" for (u, v) in arcs]
        + [f"y[{v}]" for v in range(m_rooms)]
    )

    # ---- objective: min Σ c[d] z[d] + Σ c_k y[v] -----------------------------
    c = np.zeros(n)
    for d_idx, (u, v) in enumerate(arcs):
        c[z_idx(d_idx)] = cost[u, v]
    # y[S] does not pay the hub cost (the ЦК is not a regular concentrator), but
    # the dissertation writes Σ_{v∈V} c_k y[v] without exclusion. Since the
    # constraints will pin y[S] anyway, paying the cost is harmless if we set
    # the bound h_y[S] = 1 minimum. We exclude it from the objective to keep the
    # objective faithful to the LAN-cost interpretation.
    for v in range(m_rooms):
        if v == S:
            continue
        c[y_idx(v)] = c_k

    # ---- constraints ---------------------------------------------------------
    rows_A: list[list[float]] = []
    rows_b: list[float] = []
    rows_sense: list[str] = []

    # (1.8') Σ_{d∈D_v^+} x[d] - Σ_{d∈D_v^-} x[d]  ≥  KP[v]   ∀ v ∈ V
    for v in range(m_rooms):
        row = np.zeros(n)
        for d_idx in D_in[v]:
            row[x_idx(d_idx)] += 1.0
        for d_idx in D_out[v]:
            row[x_idx(d_idx)] -= 1.0
        rows_A.append(row.tolist())
        rows_b.append(float(KP[v]))
        rows_sense.append(">=")

    # (1.9) x[d] ≤ KPS · z[d]    ⇔    x[d] - KPS · z[d] ≤ 0    ∀ d ∈ D
    for d_idx in range(n_arcs):
        row = np.zeros(n)
        row[x_idx(d_idx)] = 1.0
        row[z_idx(d_idx)] = -float(KPS)
        rows_A.append(row.tolist())
        rows_b.append(0.0)
        rows_sense.append("<=")

    # (1.10) a · y[v] ≥ KP[v] + Σ_{d∈D_v^-} z[d] - 1
    #        ⇔   a · y[v] - Σ z[d] ≥ KP[v] - 1   (only meaningful where KP[v] >= 1)
    # We exclude S and any vertex with KP_orig[v] == 0 (no users → no hub needed).
    for v in range(m_rooms):
        if v == S:
            continue
        if KP_orig[v] == 0:
            continue
        row = np.zeros(n)
        row[y_idx(v)] = float(a)
        for d_idx in D_out[v]:
            row[z_idx(d_idx)] -= 1.0
        rhs = float(KP_orig[v] - 1)
        rows_A.append(row.tolist())
        rows_b.append(rhs)
        rows_sense.append(">=")

    # (1.11) y[v] ≤ Σ_{d∈D_v^-} x[d] + KP[v] - 1
    #        ⇔   y[v] - Σ x[d] ≤ KP[v] - 1
    for v in range(m_rooms):
        if v == S:
            continue
        if KP_orig[v] == 0:
            continue
        row = np.zeros(n)
        row[y_idx(v)] = 1.0
        for d_idx in D_out[v]:
            row[x_idx(d_idx)] -= 1.0
        rhs = float(KP_orig[v] - 1)
        rows_A.append(row.tolist())
        rows_b.append(rhs)
        rows_sense.append("<=")

    A = np.array(rows_A, dtype=float)
    b = np.array(rows_b, dtype=float)

    # ---- upper bounds h ------------------------------------------------------
    h = np.zeros(n, dtype=np.int64)
    # x[d] ≤ KPS
    for d_idx in range(n_arcs):
        h[x_idx(d_idx)] = KPS
    # z[d] ≤ 1
    for d_idx in range(n_arcs):
        h[z_idx(d_idx)] = 1
    # y[v] ≤ ⌊(KP[v] + (m-1) - 1) / (a-2)⌋ + 1 (the looser of the two bounds in §1.1.2.3)
    for v in range(m_rooms):
        if v == S:
            # ЦК needs at least one stack to absorb all the traffic. Use a generous bound.
            bound = (KPS + (m_rooms - 1) - 1) // max(a - 2, 1) + 1
        elif KP_orig[v] == 0:
            bound = 0
        else:
            bound = (int(KP_orig[v]) + (m_rooms - 1) - 1) // max(a - 2, 1) + 1
        h[y_idx(v)] = max(1, bound)  # at least 1 so the LP has room

    return MILP(
        c=c,
        A=A,
        b=b,
        sense=rows_sense,
        h=h,
        direction="min",
        var_names=var_names,
    )


if __name__ == "__main__":
    p = build()
    print(f"n = {p.n}, m = {p.m}")
    print(f"non-zero objective entries: {int((p.c != 0).sum())}")
    print(f"h range: [{p.h.min()}, {p.h.max()}]")
