"""LAN-B variant: 0 off-diagonal interpreted as ZERO cost (not infinity).

This reproduces the dissertation's variable count of 153 (= 2*72 + 9) and
should match its reported optimum of 280.
"""

from tasks.lan_b import build as _base_build


def build():
    return _base_build(zero_means_infinity=False)
