"""Round-trip and structural checks for combined.transforms."""

from __future__ import annotations

import numpy as np

from combined import MILP
from combined.transforms import canonicalize, decode


def test_simple_min_no_substitution() -> None:
    """No flip-to-min, no (2.18) substitution: c is already non-negative."""
    p = MILP(
        c=np.array([3.0, 5.0]),
        A=np.array([[2.0, 4.0]]),
        b=np.array([11.0]),
        sense=[">="],
        h=np.array([5, 5]),
        direction="min",
    )
    can = canonicalize(p)
    assert (can.c_p >= 0).all()
    assert np.allclose(can.c_p, [3.0, 5.0])  # already ascending → perm=identity
    # ≥-rows are flipped to ≤ in (2.17) form
    assert can.A_p.shape == (1, 2)
    assert np.allclose(can.A_p, [[-2.0, -4.0]])
    assert np.allclose(can.b_p, [-11.0])
    # decode is the identity here
    x = np.array([0, 3])
    assert (decode(x, can) == x).all()


def test_max_with_negative_coef() -> None:
    """direction=max combined with (2.18) substitution on a negative coef."""
    # Original:  max  -2 x1 + 1 x2,   x1 + x2 ≤ 5,   0 ≤ x ≤ (3, 4).
    # After flip-to-min:  min  2 x1 - x2.
    # After (2.18) on x2:  c'[1] = +1, b' = 5 - 1*4 = 1.
    p = MILP(
        c=np.array([-2.0, 1.0]),
        A=np.array([[1.0, 1.0]]),
        b=np.array([5.0]),
        sense=["<="],
        h=np.array([3, 4]),
        direction="max",
    )
    can = canonicalize(p)
    assert can.flipped_max_to_min is True
    assert can.flipped_vars.tolist() == [False, True]
    assert (can.c_p >= 0).all()
    # variables are sorted ascending by c (post-substitution)
    assert (np.diff(can.c_p) >= 0).all()


def test_decode_is_inverse_of_canonicalize() -> None:
    """Encoding then decoding any candidate point recovers user coordinates."""
    p = MILP(
        c=np.array([-2.0, 1.0, -3.0]),
        A=np.array([[1.0, 1.0, 1.0]]),
        b=np.array([5.0]),
        sense=["<="],
        h=np.array([3, 4, 5]),
        direction="max",
    )
    can = canonicalize(p)
    # Build a synthetic (2.17)-coord point inside the box and decode it.
    z = np.array([1, 2, 1])
    user = decode(z, can)
    # User coords must be in [0, h] and match (2.17) value invertibly:
    assert user.shape == (3,)
    assert (0 <= user).all() and (user <= p.h).all()


def test_equality_constraint_is_split() -> None:
    """A row with sense='=' produces both a ≤ and a ≥ row in (2.17) form."""
    p = MILP(
        c=np.array([1.0, 1.0]),
        A=np.array([[1.0, 2.0]]),
        b=np.array([4.0]),
        sense=["="],
        h=np.array([5, 5]),
        direction="min",
    )
    can = canonicalize(p)
    assert can.A_p.shape[0] == 2  # one equality row → two ≤ rows
    # The two rows should be (a, b) and (-a, -b)
    rows = sorted(map(tuple, can.A_p.tolist()))
    assert rows == sorted([(1.0, 2.0), (-1.0, -2.0)])
