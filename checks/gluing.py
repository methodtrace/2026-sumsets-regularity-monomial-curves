"""Gluing. For four-element B' = {0, a, b, d'} and g in {2, 3} with gcd(g, d') = 1, let
A = gB' with d' added (so d = gd'). Checks, against direct computation:
  * c1 = g c1' + (g - 1)(d' - 1) and c2 = g c2' + (g - 1)(d - d' - 1);
  * if B' is not ACM: tau(A) = tau(B') + g - 1 and
    (tau - 1)d - (c1 + c2) = g((tau' - 1)d' - (c1' + c2')) + 2(g - 1);
  * if B' is ACM: A is ACM, reg(A) = reg(B') + g - 1, and N(A) = N(B') + g - 1, or N(B') + g - 2
    when c1' + c2' = 1 mod d'.
Range: all B' with d' <= 40."""
import os
import sys
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from curves.core import invariants, set_gcd

dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 40
acm_cases = exceptional = non_cases = 0
for dp in range(3, dmax + 1):
    for a in range(1, dp - 1):
        for b in range(a + 1, dp):
            if set_gcd((a, b, dp)) != 1:
                continue
            base = invariants((0, a, b, dp))
            for g in (2, 3):
                if gcd(g, dp) != 1:
                    continue
                d = g * dp
                A = tuple(sorted({0, g * a, g * b, d, dp}))
                inv = invariants(A)
                assert inv["c1"] == g * base["c1"] + (g - 1) * (dp - 1)
                assert inv["c2"] == g * base["c2"] + (g - 1) * (d - dp - 1)
                if base["acm"]:
                    assert inv["acm"] and inv["reg"] == base["reg"] + g - 1
                    if (base["c1"] + base["c2"]) % dp == 1:
                        exceptional += 1
                        assert inv["N"] == base["N"] + g - 2
                    else:
                        assert inv["N"] == base["N"] + g - 1
                    acm_cases += 1
                else:
                    assert not inv["acm"] and inv["tau"] == base["tau"] + g - 1
                    lhs = (inv["tau"] - 1) * d - (inv["c1"] + inv["c2"])
                    assert lhs == g * ((base["tau"] - 1) * dp - (base["c1"] + base["c2"])) + 2 * (g - 1)
                    non_cases += 1
print(f"gluings with d' <= {dmax}: {acm_cases} of ACM sets ({exceptional} with c1' + c2' = 1 mod d') and {non_cases} of non-ACM sets; every assertion holds")
inv = invariants((0, 80, 131, 150, 262))
print(f"example: {{0, 80, 131, 150, 262}} (from {{0, 40, 75, 131}}, g = 2): ACM {inv['acm']}, reg {inv['reg']}, sigma {inv['N']}")
