"""The band-shaped family of Theorem 5.3 in P^4:
A = {0, nm, n(m + 1), nQ, d}, Q = 2m + 1, d = m^2 + nQ, gcd(m, n) = 1.

With R = m + n - 1, by direct computation for all coprime m, n <= 18:
the curve is ACM; the homogenized Apery set is
{(n(mi + (m + 1)j), d*max(i, j) - n(mi + (m + 1)j)) : 0 <= i, j <= R, |i - j| <= m};
reg = R; c1 = nQR - d + 1; c2 = m^2 R + nm(m + 1) - d + 1; E = nm(m + 1);
and sigma = R - 2 + ceil((E + 2)/d) (for an ACM curve sigma = N).
Also reports the member with reg < sigma of least d, confirmed by direct computation for every
coprime pair (m, n) with d <= 101 (m <= 9, n <= 33). Then tabulates the Pell sequence
(3, 2), (m, n) -> (3m + 4n, 2m + 3n) and the subfamily (m, n) = (2j + 1, j) from the
proved formulas."""
import os
import sys
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from curves.core import apery_points, invariants


def formulas(m, n):
    Q = 2 * m + 1
    d = m * m + n * Q
    R = m + n - 1
    E = n * m * (m + 1)
    return d, R, n * Q * R - d + 1, m * m * R + E - d + 1, E, R - 2 + -(-(E + 2) // d)


count = below = 0
first = None
for m in range(1, 19):
    for n in range(1, 19):
        if gcd(m, n) != 1:
            continue
        Q = 2 * m + 1
        d, R, c1, c2, E, sigma = formulas(m, n)
        A = (0, n * m, n * (m + 1), n * Q, d)
        inv = invariants(A)
        assert inv["acm"] and inv["reg"] == R
        band = sorted((n * (m * i + (m + 1) * j), d * max(i, j) - n * (m * i + (m + 1) * j))
                      for i in range(R + 1) for j in range(R + 1) if abs(i - j) <= m)
        assert apery_points(A) == band
        assert (inv["c1"], inv["c2"], inv["E"], inv["N"]) == (c1, c2, E, sigma)
        count += 1
        if R < sigma:
            below += 1
            if first is None or d < first[0]:
                first = (d, m, n, A, R, sigma)
print(f"coprime pairs m, n <= 18: {count}; all formulas confirmed; reg < sigma in {below}")
d, m, n, A, R, sigma = first
small = 0
for m2 in range(1, 10):
    for n2 in range(1, 34):
        d2 = m2 * m2 + n2 * (2 * m2 + 1)
        if gcd(m2, n2) != 1 or d2 > d:
            continue
        inv2 = invariants((0, n2 * m2, n2 * (m2 + 1), n2 * (2 * m2 + 1), d2))
        assert inv2["acm"]
        assert (inv2["reg"] < inv2["N"]) == ((m2, n2) == (m, n)), (m2, n2)
        small += 1
print(f"least d with reg < sigma: d = {d}, (m, n) = ({m}, {n}), A = {A}, reg {R}, sigma {sigma}; "
      f"the only one among all {small} coprime pairs with d <= {d}")
print("Pell sequence (m, n): d, reg, c1, c2, sigma")
m, n = 3, 2
for _ in range(6):
    assert m * m - 2 * n * n == 1
    d, R, c1, c2, E, sigma = formulas(m, n)
    print(f"  ({m}, {n}): {d}, {R}, {c1}, {c2}, {sigma}")
    m, n = 3 * m + 4 * n, 2 * m + 3 * n
print("subfamily (m, n) = (2j + 1, j): reg = 3j, sigma = 3j + floor(j/2) - 1, checked for j <= 9999")
for j in range(1, 10000):
    d, R, c1, c2, E, sigma = formulas(2 * j + 1, j)
    assert d == 8 * j * j + 7 * j + 1 and R == 3 * j and E == 4 * j ** 3 + 6 * j ** 2 + 2 * j
    assert sigma == 3 * j + j // 2 - 1
print(f"  j = 8: A = {(0, 8 * 17, 8 * 18, 8 * 35, 569)}, reg 24, sigma 27")
inv = invariants((0, 136, 144, 280, 569))
assert inv["acm"] and inv["reg"] == 24 and inv["N"] == 27
