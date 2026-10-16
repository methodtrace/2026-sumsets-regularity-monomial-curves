"""The defect of an ACM curve in any P^r: E = M1 + M2 - reg*d >= 0 and sigma - reg = ceil((E + 2)/d) - 2,
where M1, M2 are the largest Apery elements of S1, S2 with respect to d.

Checks the identity for every ACM set with four elements and d <= 60 and five elements and d <= 50,
the conductors against direct reachability, and reports the largest value of E/(d(reg + 2)) with
and without a lower bound on reg (in P^3, 9E <= d(reg + 2))."""
import os
import sys
from fractions import Fraction
from itertools import combinations
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from curves.core import invariants, conductor, set_gcd

for k, dmax, rmin in ((4, 60, 0), (4, 60, 4), (5, 50, 0), (5, 50, 6)):
    count = 0
    best = None
    for d in range(k - 1, dmax + 1):
        for mid in combinations(range(1, d), k - 2):
            if set_gcd(mid + (d,)) != 1:
                continue
            A = (0,) + mid + (d,)
            inv = invariants(A)
            if not inv["acm"]:
                continue
            count += 1
            E, R, N = inv["E"], inv["reg"], inv["N"]
            assert E >= 0 and N - R == -(-(E + 2) // d) - 2
            if rmin == 0:
                assert inv["c1"] == conductor(A[1:]) and inv["c2"] == conductor([d - a for a in A[:-1]])
            if R < rmin:
                continue
            q = Fraction(E, d * (R + 2))
            if best is None or q > best[0]:
                best = (q, A, R, E)
    q, A, R, E = best
    print(f"{k} elements, d <= {dmax}, reg >= {rmin}: ACM sets {count}; largest E/(d(reg + 2)) = {q} = {float(q):.4f} at {A} (reg {R}, E {E})")
