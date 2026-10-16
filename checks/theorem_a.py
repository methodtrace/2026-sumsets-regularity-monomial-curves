"""Non-ACM curves in P^3 (Theorem A and its equality case), and the fan.

For every four-element set A = {0, alpha, beta, d} with gcd 1 and d <= DMAX (default 150):
  * fan: the conductor formulas from the fan equal the conductors; the curve is non-ACM exactly
    when nu < mu; for non-ACM sets tau >= T - 3 with T = (s_mu - r_mu) + (s_mu - s_{mu+1} + p_{mu+1} - p_mu),
    and the full maximum over nu + 1 <= i <= mu of deg g_i + max(deg f_i, deg f_{i-1}) - 2 in the
    proof of Lemma 3.2 equals tau + 1 (the regularity of the ideal);
  * Theorem A: for non-ACM sets reg = tau and (tau - 1)d - (c1 + c2 + d) >= 0;
  * equality holds exactly for {0, b - 1, b + 1, 2b} with b even;
  * the number of non-ACM sets with tau = N + 2 but positive unrounded slack;
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from curves.core import invariants, set_gcd
from curves.fan import fan_data

dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 150
total = nonacm = rounded_only = 0
equality = []
min_slack = None
for d in range(3, dmax + 1):
    for a in range(1, d - 1):
        for b in range(a + 1, d):
            if set_gcd((a, b, d)) != 1:
                continue
            A = (0, a, b, d)
            inv = invariants(A)
            total += 1
            xi1, eta1, xi2, eta2, s, p, r, rp, v, nu, mu = fan_data(a, b, d)
            fr1 = -d + xi1 * (s[v] - 1) + eta1 * (p[v + 1] - 1) - min(xi1 * s[v + 1], eta1 * p[v])
            fr2 = -d + xi2 * (s[mu] - 1) + eta2 * (p[mu + 1] - 1) - min(xi2 * s[mu + 1], eta2 * p[mu])
            assert fr1 + 1 == inv["c1"] and fr2 + 1 == inv["c2"], A
            assert (nu < mu) == (not inv["acm"]), A
            if inv["acm"]:
                continue
            nonacm += 1
            T = (s[mu] - r[mu]) + (s[mu] - s[mu + 1] + p[mu + 1] - p[mu])
            assert inv["tau"] >= T - 3, A
            # the full maximum in the proof of BGM Thm 2.8 (regularity of the ideal) is tau + 1
            deg_g = [s[k] + max(0, -r[k]) + max(0, -rp[k]) for k in range(len(s))]
            deg_f = [s[k] - s[k + 1] + p[k + 1] - p[k] for k in range(len(s) - 1)]
            bgm = max(deg_g[k] + max(deg_f[k], deg_f[k - 1]) for k in range(nu + 1, mu + 1)) - 2
            assert bgm == inv["tau"] + 1, A
            assert inv["reg"] == inv["tau"], A
            slack = (inv["tau"] - 1) * d - (inv["c1"] + inv["c2"] + d)
            assert slack >= 0, A
            min_slack = slack if min_slack is None else min(min_slack, slack)
            if slack == 0:
                equality.append(A)
            elif inv["tau"] == inv["N"] + 2:
                rounded_only += 1
expected = [(0, b - 1, b + 1, 2 * b) for b in range(2, dmax // 2 + 1, 2)]
assert equality == expected, equality
print(f"four-element sets with d <= {dmax}: {total}, of which non-ACM: {nonacm}")
print("fan conductor formulas, non-ACM criterion nu < mu, and tau >= T - 3: all hold")
print("for every non-ACM set the full maximum in the proof of Lemma 3.2 equals tau + 1")
print(f"Theorem A: reg = tau and (tau - 1)d >= c1 + c2 + d for every non-ACM set; least slack {min_slack}")
print(f"equality exactly for {{0, b - 1, b + 1, 2b}}, b even: {len(equality)} sets, the first {equality[:3]}")
print(f"non-ACM sets with tau = N + 2 and positive unrounded slack: {rounded_only}")
