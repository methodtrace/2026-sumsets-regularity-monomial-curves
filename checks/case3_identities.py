"""Case 3 in the proof of Theorem A: the four expressions D_11, D_12, D_21, D_22.

Notation of the paper: X = r_nu, Y = -r_mu, e = r'_mu - 1, f = -r'_{mu+1}, w = s_mu - s_{mu+1},
t = p_{mu+1} - p_mu, pbar = p_mu, sbar = s_{mu+1}, with
  e + 1 = sbar + w - pbar + Y,   d = pbar*w + sbar*t + t*w,
  xi1 = (X + Y)pbar - Y t,       eta1 = (X + Y)sbar + (X + 2Y)w,
  xi2 = (e + 1 + f)pbar + (e + 1)t,  eta2 = (e + 1 + f)sbar + f w.
For each choice of the two maxima, d(pbar + e) - 2 - (choice) equals the stated nonnegative form
minus (sbar + w)(xi1 + xi2 - d), and xi1 + xi2 - d = pbar(X + 2Y + f + sbar - pbar - t), which
vanishes by relation (3.4). Exact polynomial arithmetic over the integers, standard library only."""


def const(c):
    return {(): c} if c else {}


def var(name):
    return {((name, 1),): 1}


def add(a, b, sign=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + sign * v
        if out[k] == 0:
            del out[k]
    return out


def mul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            exps = dict(ka)
            for n, x in kb:
                exps[n] = exps.get(n, 0) + x
            k = tuple(sorted(exps.items()))
            out[k] = out.get(k, 0) + va * vb
            if out[k] == 0:
                del out[k]
    return out


def total(*terms):
    out = {}
    for t_ in terms:
        out = add(out, t_)
    return out


def prod(*factors):
    out = const(1)
    for f_ in factors:
        out = mul(out, f_)
    return out


pbar, sbar, w, t, X, Y, f = (var(n) for n in ("pbar", "sbar", "w", "t", "X", "Y", "f"))
e = total(sbar, w, mul(pbar, const(-1)), Y, const(-1))                 # e + 1 = sbar + w - pbar + Y
e1 = add(e, const(1))
d = total(prod(pbar, w), prod(sbar, t), prod(t, w))
xi1 = add(mul(add(X, Y), pbar), mul(Y, t), -1)
eta1 = add(mul(add(X, Y), sbar), mul(total(X, Y, Y), w))
xi2 = add(mul(total(e1, f), pbar), mul(e1, t))
eta2 = add(mul(total(e1, f), sbar), mul(f, w))
rhs = add(mul(d, add(pbar, e)), const(-2))
s_mu = add(sbar, w)            # s_mu
p_nu = add(pbar, t, -1)        # p_nu
choices = {
    "D_11": add(mul(xi1, s_mu), mul(xi2, sbar)),
    "D_12": add(mul(xi1, s_mu), mul(eta2, pbar)),
    "D_21": add(mul(eta1, p_nu), mul(xi2, sbar)),
    "D_22": add(mul(eta1, p_nu), mul(eta2, pbar)),
}
m2 = const(-2)
stated = {
    "D_11": total(prod(pbar, w, e), prod(pbar, w, Y), prod(pbar, w, f), prod(sbar, t, add(Y, const(-1))),
                  prod(t, w, e), prod(t, w, Y), m2),
    "D_12": total(prod(pbar, w, e), prod(pbar, w, Y), prod(sbar, t, e), prod(sbar, t, Y),
                  prod(t, w, e), prod(t, w, Y), m2),
    "D_21": total(prod(pbar, w, e), prod(pbar, w, f), prod(sbar, t, total(X, Y, const(-1))),
                  prod(t, w, e), prod(t, w, X), prod(const(2), t, w, Y), m2),
    "D_22": total(prod(pbar, w, e), prod(sbar, t, e), prod(sbar, t, X), prod(sbar, t, Y),
                  prod(t, w, e), prod(t, w, X), prod(const(2), t, w, Y), m2),
}
relation = add(add(xi1, xi2), d, -1)
assert relation == mul(pbar, total(X, mul(const(2), Y), f, sbar, mul(pbar, const(-1)), mul(t, const(-1))))
print("xi1 + xi2 - d = pbar(X + 2Y + f + sbar - pbar - t), which is zero by relation (3.4)")
for name in ("D_11", "D_12", "D_21", "D_22"):
    diff = add(add(rhs, choices[name], -1), stated[name], -1)
    assert diff == mul(mul(s_mu, const(-1)), relation), name
    print(f"{name}: d(pbar + e) - 2 - (choice) = stated form - (sbar + w)(xi1 + xi2 - d)")
