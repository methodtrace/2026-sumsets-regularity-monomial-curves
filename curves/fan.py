"""The fan of Bermejo, Gimenez and Morales for a four-element set {0, alpha, beta, d}."""
from math import gcd


def fan_data(alpha, beta, d):
    """For 0 < alpha < beta < d, return (xi1, eta1, xi2, eta2, s, p, r, rp, v, nu, mu) with
    (xi1, eta1) = (alpha, beta), (xi2, eta2) = (d - alpha, d - beta) and the fan vectors
    (s_i, p_i), i = -1, ..., l + 1, stored at list positions 0, ..., l + 2."""
    assert 0 < alpha < beta < d
    xi1, eta1 = alpha, beta
    xi2, eta2 = d - alpha, d - beta
    g = gcd(xi1, d)
    s = [d // g]
    p = [0]
    s0 = next(t for t in range(d // g) if (t * xi1 - g * eta1) % d == 0)
    s.append(s0)
    p.append(g)
    while s[-1] != 0:
        q = -(-s[-2] // s[-1])
        s.append(q * s[-1] - s[-2])
        p.append(q * p[-1] - p[-2])
    r = [(si * xi1 - pi * eta1) // d for si, pi in zip(s, p)]
    rp = [(si * xi2 - pi * eta2) // d for si, pi in zip(s, p)]
    n = len(s)
    v = max(k for k in range(n) if r[k] > 0)
    nu = max(k for k in range(n) if r[k] >= 0)
    mu = max(k for k in range(n) if rp[k] > 0)
    return xi1, eta1, xi2, eta2, s, p, r, rp, v, nu, mu
