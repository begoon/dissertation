"""Problem representation matching the dissertation's (2.1).

Variables x ∈ N^n with 0 ≤ x ≤ h, mixed-row linear constraints, linear objective.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

Sense = Literal["<=", ">=", "="]


@dataclass
class MILP:
    """Mixed-integer linear program in the dissertation's (2.1) form.

    min/max  c · x
       s.t.  A x  {<=, >=, =}  b   (per-row sense)
             0 ≤ x ≤ h
             x ∈ Z^n_+
    """

    c: np.ndarray              # shape (n,)
    A: np.ndarray              # shape (m, n)
    b: np.ndarray              # shape (m,)
    sense: list[Sense]         # length m
    h: np.ndarray              # shape (n,) integer
    direction: Literal["min", "max"] = "min"
    var_names: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.c = np.asarray(self.c, dtype=float)
        self.A = np.asarray(self.A, dtype=float)
        self.b = np.asarray(self.b, dtype=float)
        self.h = np.asarray(self.h, dtype=np.int64)
        if not self.var_names:
            self.var_names = [f"x[{j}]" for j in range(self.n)]
        assert self.A.shape == (self.m, self.n), (self.A.shape, self.m, self.n)
        assert self.b.shape == (self.m,), self.b.shape
        assert self.h.shape == (self.n,), self.h.shape
        assert len(self.sense) == self.m
        assert len(self.var_names) == self.n

    @property
    def n(self) -> int:
        return self.c.shape[0]

    @property
    def m(self) -> int:
        return len(self.sense)
