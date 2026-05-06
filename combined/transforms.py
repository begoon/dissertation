"""Canonical-form transforms (§3 of COMBINED_METHOD.md).

Two canonical forms are produced from the user MILP:

* (2.17) — for the lattice search:
    min c^p · x^p,  A^p x^p ≤ b^p,  0 ≤ x^p ≤ h^p,  x^p ∈ Z^n_+,
    c^p[j] ≥ 0 ∀ j,  c^p[1] ≤ ... ≤ c^p[n].

* (2.19) — for the LP relaxation:
    same problem in equality form (no flips, no equality splits;
    slacks added per row), to feed scipy.linprog.

The two share the same variable order (post-substitution and post-sort).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from combined.model import MILP


@dataclass
class Canonical:
    """Both canonical forms plus the inverse-decode metadata."""

    # ----- shared: what was done to the input -----
    flipped_max_to_min: bool                    # whether c was negated up front
    flipped_vars: np.ndarray                    # bool, len original n; True iff x[j] = h[j] - x^p[j] applied
    perm: np.ndarray                            # int, len n; perm[j_new] = j_old (after the sort)
    h_orig: np.ndarray                          # original h, before perm

    # ----- (2.17) form for the lattice -----
    c_p: np.ndarray                             # shape (n,)
    A_p: np.ndarray                             # shape (m_p, n)
    b_p: np.ndarray                             # shape (m_p,)
    h_p: np.ndarray                             # shape (n,) integer

    # ----- (2.19) form for the LP -----
    # min c_l_for_min · x_l (we'll feed scipy as min; it's c^p padded with zeros for slacks).
    A_l: np.ndarray                             # shape (m, n + m)  — n structural + m slacks (one per original row)
    b_l: np.ndarray                             # shape (m,)
    h_l: np.ndarray                             # shape (n + m,) — slacks have effectively-infinite upper bound
    c_l_for_min: np.ndarray                     # shape (n + m,) — c^p in first n positions, 0 in slack positions
    row_sense: list[str]                        # original senses of the m rows (for reference)
    n_struct: int                               # = original n (= len(c_p))


def canonicalize(p: MILP) -> Canonical:
    """Apply (2.18) substitution + sort + equality split + ≥→≤ flip; build A^l for the LP.

    Returns a Canonical bundling everything Stages 1–4 will use.
    """
    n = p.n
    m = p.m

    # 1. Force min by flipping c if needed.
    flipped_max_to_min = p.direction == "max"
    c0 = -p.c if flipped_max_to_min else p.c.copy()
    A0 = p.A.copy()
    b0 = p.b.copy()
    h0 = p.h.copy()

    # 2. Substitution (2.18): for j with c0[j] < 0, x[j] = h[j] - x^p[j].
    #    In ANY linear expression Σ a_ij x[j] = ... this maps to
    #      Σ a_ij (h[j] - x^p[j]) = Σ a_ij h[j] - Σ a_ij x^p[j].
    #    For the objective: c'[j] = -c0[j], constant term Σ c0[j] h[j] dropped.
    #    For row i (any sense): a'[i,j] = -a0[i,j], b'[i] = b0[i] - a0[i,j] h[j].
    flipped_vars = c0 < 0
    if flipped_vars.any():
        # Update b: for every row i and every flipped j, subtract a0[i,j] * h[j]
        delta_b = A0[:, flipped_vars] @ h0[flipped_vars]
        b0 = b0 - delta_b
        # Negate columns of A and the corresponding c entries
        A0[:, flipped_vars] = -A0[:, flipped_vars]
        c0[flipped_vars] = -c0[flipped_vars]

    # 3. Sort variables ascending by c0.
    perm = np.argsort(c0, kind="stable")        # perm[j_new] = j_old
    c_sorted = c0[perm]
    A_sorted = A0[:, perm]
    h_sorted = h0[perm]

    # ----- LP form (2.19): keep mixed senses, add one slack per row, A_l = [A_sorted | S]
    # For sense '<=':  Σ a x + s = b,  s ≥ 0, no upper bound (effectively)
    # For sense '>=':  Σ a x - s = b,  s ≥ 0
    # For sense '=':                                   s coefficient = 0 (or no slack at all)
    slack_coefs = np.zeros((m, m))
    slack_h = np.full(m, 1e18)                  # effectively-∞ upper bound for the LP slack
    for i, s in enumerate(p.sense):
        if s == "<=":
            slack_coefs[i, i] = 1.0
        elif s == ">=":
            slack_coefs[i, i] = -1.0
        elif s == "=":
            slack_coefs[i, i] = 0.0             # slack column is all zeros — slack is dummy
            slack_h[i] = 0.0                    # pin to zero
        else:
            raise ValueError(f"unknown sense {s!r}")
    A_l = np.hstack([A_sorted, slack_coefs])
    b_l = b0.copy()
    c_l_for_min = np.concatenate([c_sorted, np.zeros(m)])
    h_l = np.concatenate([h_sorted.astype(float), slack_h])

    # ----- (2.17) form for the lattice search: split equalities, flip ≥, all rows ≤ -----
    rows_a: list[np.ndarray] = []
    rows_b: list[float] = []
    for i, s in enumerate(p.sense):
        ai = A_sorted[i]
        bi = b0[i]
        if s == "<=":
            rows_a.append(ai); rows_b.append(bi)
        elif s == ">=":
            rows_a.append(-ai); rows_b.append(-bi)
        elif s == "=":
            rows_a.append(ai);  rows_b.append(bi)
            rows_a.append(-ai); rows_b.append(-bi)
    A_p = np.array(rows_a, dtype=float)
    b_p = np.array(rows_b, dtype=float)
    h_p = h_sorted.astype(np.int64)
    c_p = c_sorted.astype(float)

    # Sanity: c_p must be non-negative (sort doesn't change that) and ascending.
    assert (c_p >= -1e-12).all(), c_p
    assert (np.diff(c_p) >= -1e-12).all(), c_p

    return Canonical(
        flipped_max_to_min=flipped_max_to_min,
        flipped_vars=flipped_vars,
        perm=perm,
        h_orig=p.h.copy(),
        c_p=c_p,
        A_p=A_p,
        b_p=b_p,
        h_p=h_p,
        A_l=A_l,
        b_l=b_l,
        h_l=h_l,
        c_l_for_min=c_l_for_min,
        row_sense=list(p.sense),
        n_struct=n,
    )


def decode(x_p: np.ndarray, can: Canonical) -> np.ndarray:
    """Map a solution from (2.17) back to the user's variable order/sign.

    Inverse of: sort by c, then x[j] = h[j] - x^p[j] for flipped j.
    """
    n = can.n_struct
    # Undo the permutation: x_sorted[j_new] = x_orig_after_flip[perm[j_new]]
    x_sorted_to_orig = np.empty(n, dtype=x_p.dtype)
    x_sorted_to_orig[can.perm] = x_p
    # Undo the substitution: where flipped, x = h - x^p
    x = x_sorted_to_orig.copy()
    flipped = can.flipped_vars
    x[flipped] = can.h_orig[flipped] - x_sorted_to_orig[flipped]
    return x
