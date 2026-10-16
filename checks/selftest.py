"""Self-test: invariants from the homogenized Apery set agree with direct sumsets.

For every four-element set with d <= 40 and every five-element set with d <= 24: tau, c1 and c2
from curves.core.invariants equal tau from comparing hA with its intersection form level by
level and the conductors from direct reachability."""
import os
import sys
from itertools import combinations
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from curves.core import invariants, brute_invariants, set_gcd

for k, dmax in ((4, 40), (5, 24)):
    count = 0
    for d in range(k - 1, dmax + 1):
        for mid in combinations(range(1, d), k - 2):
            if set_gcd(mid + (d,)) != 1:
                continue
            A = (0,) + mid + (d,)
            fast = invariants(A)
            slow = brute_invariants(A)
            assert (fast["tau"], fast["c1"], fast["c2"]) == (slow["tau"], slow["c1"], slow["c2"]), A
            count += 1
    print(f"{k}-element sets with d <= {dmax}: {count}; Apery-set invariants equal direct sumsets in every case")
