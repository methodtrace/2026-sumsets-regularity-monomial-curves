"""The P^3 family with reg < sigma: A = {0, KI, (J + K)I, IJ + K(I + J)}, I, J, K pairwise coprime.

  * by direct computation for all such triples with I, J, K <= 11 (so d <= 299): ACM,
    reg = I + J + K - 2, E = IJK, c1 + c2 = reg*d + IJK - 2d + 2,
    sigma - reg = ceil((IJK + 2)/d) - 2, and reg < sigma exactly when IJK >= 2d - 1;
  * the odd subfamily (n, n + 1, n + 2): by direct computation for odd n <= 15, and by exact
    integer arithmetic for odd n <= 9999, c1 + c2 = 10n^3 + 18n^2 + 2n, d = 3n^2 + 6n + 2,
    reg = 3n + 1 and N = floor(10n/3).
"""
import os
import sys
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from curves.core import invariants

count = below = 0
smallest = None
for I in range(1, 12):
    for J in range(1, 12):
        for K in range(1, 12):
            if gcd(I, J) != 1 or gcd(I, K) != 1 or gcd(J, K) != 1:
                continue
            d = I * J + K * (I + J)
            inv = invariants((0, K * I, (J + K) * I, d))
            R = I + J + K - 2
            assert inv["acm"] and inv["reg"] == R and inv["E"] == I * J * K
            assert inv["c1"] + inv["c2"] == R * d + I * J * K - 2 * d + 2
            assert inv["N"] - R == -(-(I * J * K + 2) // d) - 2
            assert (R < inv["N"]) == (I * J * K >= 2 * d - 1)
            count += 1
            if R < inv["N"]:
                below += 1
                if smallest is None or d < smallest[0]:
                    smallest = (d, (0, K * I, (J + K) * I, d), (I, J, K))
print(f"triples checked: {count}; with reg < sigma: {below}; smallest d: {smallest[0]}, e.g. {smallest[1]} from (I, J, K) = {smallest[2]}")
for n in range(1, 16, 2):
    d = 3 * n * n + 6 * n + 2
    inv = invariants((0, n * (n + 2), n * (2 * n + 3), d))
    assert inv["acm"] and inv["reg"] == 3 * n + 1 and inv["c1"] + inv["c2"] == 10 * n ** 3 + 18 * n ** 2 + 2 * n
    assert inv["N"] == (10 * n) // 3
for n in range(1, 10000, 2):
    d = 3 * n * n + 6 * n + 2
    c = 10 * n ** 3 + 18 * n ** 2 + 2 * n
    assert c == (3 * n + 1) * d + n * (n + 1) * (n + 2) - 2 * d + 2
    assert -(-c // d) == (10 * n) // 3
print("odd subfamily: formulas confirmed by direct computation for odd n <= 15 and by integer arithmetic for odd n <= 9999")
for n in (1, 3, 5, 7, 9, 99, 999, 9999):
    print(f"  n = {n}: reg {3 * n + 1}, sigma {(10 * n) // 3}")
