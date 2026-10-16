"""ACM curves in P^3: the six-parameter description, regularity and conductors, Theorem B, and the
sets with reg < sigma.

For every ACM four-element set A = {0, alpha, beta, d} with d <= DMAX (default 131):
  * the homogenized Apery set is {i*B + j*A' : (i, j) in D} for B = (beta, d - beta),
    A' = (alpha, d - alpha) and D = {0 <= i < l, 0 <= j < h, i < l - w or j < h - y};
    l, h, w, y, u, v are recovered from the Apery set and satisfy d = lh - wy,
    alpha = wu + lv, beta = hu + yv and the stated inequalities;
  * reg = l + h - min(w, y) - 2 and the conductor formulas at the two corners;
  * after reflecting so that w >= y: c1 + c2 = (reg - 2)d + 2 + E with
    E = max(0, y(d - alpha) - w(d - beta)), 0 <= E <= y(l - w)(h - w), 9E <= d(reg + 2);
  * 9N <= 10 reg - 2, and sigma - reg = N - reg = ceil((E + 2)/d) - 2.
Also 9N <= 10 reg - 7 (Remark 4.7) and no set with reg < sigma below d = 131.
Reports the least value of 10 reg - 2 - 9N, the sets attaining it, and every set with reg < sigma.
Here (l, h, w, y) are the parameters (l1, l2, w1, w2) of Theorem 4.1.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from curves.core import apery_points, invariants, set_gcd


def parameters(a, b, d, pts):
    """Recover (l, h, w, y, u, v) from the Apery points of an ACM set {0, a, b, d}."""
    coords = {}
    for x, yy in pts:
        s = (x + yy) // d
        i = (x - s * a) // (b - a)
        assert i * b + (s - i) * a == x
        coords[(i, s - i)] = x
    D = set(coords)
    l = min(i for i in range(d + 2) if (i, 0) not in D)
    h = min(j for j in range(d + 2) if (0, j) not in D)
    notch = [(i, j) for i in range(l) for j in range(h) if (i, j) not in D]
    res_j = {coords[(0, j)] % d: j for j in range(h)}
    res_i = {coords[(i, 0)] % d: i for i in range(l)}
    yp = res_j[(l * b) % d]
    wp = res_i[(h * a) % d]
    if notch:
        a0 = min(i for i, j in notch)
        b0 = min(j for i, j in notch)
        w, y = l - a0, h - b0
        assert len(notch) == w * y and (w, y) == (wp, yp)
    else:
        w, y = wp, yp
        assert w * y == 0
    u = (l * b - y * a) // d
    v = (h * a - w * b) // d
    assert l * b - y * a == u * d and h * a - w * b == v * d
    expected = {(i, j) for i in range(l) for j in range(h) if i < l - w or j < h - y}
    assert D == expected
    return l, h, w, y, u, v


dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 131
acm = 0
min_slack = None
below = []
for d in range(3, dmax + 1):
    for a in range(1, d - 1):
        for b in range(a + 1, d):
            if set_gcd((a, b, d)) != 1:
                continue
            A = (0, a, b, d)
            inv = invariants(A)
            if not inv["acm"]:
                continue
            acm += 1
            R, c1, c2, N, E = inv["reg"], inv["c1"], inv["c2"], inv["N"], inv["E"]
            l, h, w, y, u, v = parameters(a, b, d, apery_points(A))
            assert d == l * h - w * y and a == w * u + l * v and b == h * u + y * v
            assert 0 <= w < min(l, h) and 0 <= y < min(l, h) and 1 <= u <= l - y and 0 <= v <= h - w
            assert R == l + h - min(w, y) - 2
            assert c1 == (l - 1) * b + (h - 1) * a - min(w * b, y * a) - d + 1
            assert c2 == (l - 1) * (d - b) + (h - 1) * (d - a) - min(w * (d - b), y * (d - a)) - d + 1
            if w >= y:
                W, Yv, al, be, H, L = w, y, a, b, h, l
            else:
                W, Yv, al, be, H, L = y, w, d - b, d - a, l, h
            E_formula = max(0, Yv * (d - al) - W * (d - be))
            assert E_formula == E and c1 + c2 == (R - 2) * d + 2 + E
            assert 0 <= E <= Yv * (L - W) * (H - W)
            assert 9 * E <= d * (R + 2)
            assert 9 * N <= 10 * R - 2
            assert 9 * N <= 10 * R - 7
            assert N - R == -(-(E + 2) // d) - 2
            slack = 10 * R - 2 - 9 * N
            if min_slack is None or slack < min_slack[0]:
                min_slack = (slack, A, R, N)
                minimizers = [A]
            elif slack == min_slack[0]:
                minimizers.append(A)
            if R < N:
                below.append((A, R, N))
print(f"ACM four-element sets with d <= {dmax}: {acm}; every assertion holds")
print(f"least value of 10 reg - 2 - 9N: {min_slack[0]}, at {min_slack[1]} (reg {min_slack[2]}, N {min_slack[3]})")
print(f"sets attaining it: {len(minimizers)}, with d in {sorted(set(A[3] for A in minimizers))}: {minimizers}")
assert all(A[3] >= 131 for A, _, _ in below)
if dmax >= 131:
    assert len(minimizers) == 6 and {A[3] for A in minimizers} == {31}
    assert len(below) == 6 and all(A[3] == 131 for A, _, _ in below)
print(f"sets with reg < sigma: {len(below)}")
for A, R, N in below:
    print(f"  {A}: reg {R}, sigma {N}")
