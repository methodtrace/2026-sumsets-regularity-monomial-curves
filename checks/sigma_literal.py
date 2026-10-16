"""The counterexamples do not depend on the central-interval convention for sigma.

Definition 2.2 reads the structure decomposition (1.1) with its central interval present
(hd - c2 >= c1). Under the literal reading the central interval [c1, hd - c2] may be empty:
level h satisfies the literal decomposition when hA is the union of S1 cap [0, c1 - 2],
[c1, hd - c2] and hd - (S2 cap [0, c2 - 2]). Let sigma_lit be the least h0 such that every
level h >= h0 satisfies it. Every level h >= sigma satisfies it (the central interval is then
present), so sigma_lit is one more than the largest level below sigma that fails, and
sigma_lit <= sigma.

For every set with reg < sigma in the ranges of Table 1 (the four-element sets of Theorem 4.5
with I, J, K <= 11, which include the six sets with d = 131, and the band family of
Theorem 5.3 with m, n <= 18), this script computes sigma_lit by direct sumsets and asserts
reg < sigma_lit, so each set is a counterexample to sigma <= reg under either reading.
"""
import os
import sys
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from curves.core import invariants, semigroup_members


def reverse_bits(x, width):
    return int(format(x, "0%db" % width)[::-1], 2) if x else 0


def sigma_literal(A, sigma):
    """sigma_lit by direct sumsets at levels 1, ..., sigma - 1."""
    d = A[-1]
    inv = invariants(A)
    c1, c2 = inv["c1"], inv["c2"]
    lim = sigma * d
    s1 = semigroup_members(A[1:], lim)
    s2 = semigroup_members([d - a for a in A[:-1]], lim)
    low1 = s1 & ((1 << max(c1 - 1, 0)) - 1)          # S1 cap [0, c1 - 2]
    low2 = s2 & ((1 << max(c2 - 1, 0)) - 1)          # S2 cap [0, c2 - 2]
    cur = 1
    last_fail = 0
    for h in range(1, sigma):
        nxt = 0
        for a in A:
            nxt |= cur << a
        cur = nxt
        hd = h * d
        if c2 - 2 > hd:                               # hd - (S2 cap [0, c2 - 2]) has negatives
            last_fail = h
            continue
        upper = reverse_bits(low2, hd + 1)            # {hd - y : y in S2 cap [0, c2 - 2]}
        middle = 0
        if hd - c2 >= c1:
            middle = ((1 << (hd - c2 + 1)) - 1) ^ ((1 << c1) - 1)
        if low1 | middle | upper != cur:
            last_fail = h
    return last_fail + 1


def run(sets):
    least_gap = None
    smaller = 0
    out = []
    for A in sets:
        inv = invariants(A)
        assert inv["acm"]
        reg, sigma = inv["reg"], inv["N"]
        assert reg < sigma
        lit = sigma_literal(A, sigma)
        assert reg < lit <= sigma, (A, reg, lit, sigma)
        least_gap = lit - reg if least_gap is None else min(least_gap, lit - reg)
        if lit < sigma:
            smaller += 1
        out.append((A, reg, lit, sigma))
    return out, least_gap, smaller


family = []
for I in range(1, 12):
    for J in range(1, 12):
        for K in range(1, 12):
            if gcd(I, J) != 1 or gcd(I, K) != 1 or gcd(J, K) != 1:
                continue
            d = I * J + K * (I + J)
            if I * J * K >= 2 * d - 1:
                family.append((0, K * I, (J + K) * I, d))
family.sort(key=lambda A: (A[-1], A))
res, gap, smaller = run(family)
print(f"Theorem 4.5, I, J, K <= 11: {len(res)} sets with reg < sigma; reg < sigma_lit for all; "
      f"least sigma_lit - reg {gap}; sigma_lit < sigma for {smaller}")
for A, reg, lit, sigma in res:
    if A[-1] == 131:
        print(f"  {A}: reg {reg}, sigma_lit {lit}, sigma {sigma}")
print(f"  sets with sigma_lit < sigma: d = {sorted(set(A[-1] for A, reg, lit, sigma in res if lit < sigma))}")

band = []
for m in range(1, 19):
    for n in range(1, 19):
        if gcd(m, n) != 1:
            continue
        Q = 2 * m + 1
        d = m * m + n * Q
        R = m + n - 1
        E = n * m * (m + 1)
        if R - 2 + -(-(E + 2) // d) > R:
            band.append((0, n * m, n * (m + 1), n * Q, d))
band.sort(key=lambda A: (A[-1], A))
res, gap, smaller = run(band)
print(f"Theorem 5.3, m, n <= 18: {len(res)} sets with reg < sigma; reg < sigma_lit for all; "
      f"least sigma_lit - reg {gap}; sigma_lit < sigma for {smaller}")
A, reg, lit, sigma = res[0]
print(f"  least d: {A}: reg {reg}, sigma_lit {lit}, sigma {sigma}")
