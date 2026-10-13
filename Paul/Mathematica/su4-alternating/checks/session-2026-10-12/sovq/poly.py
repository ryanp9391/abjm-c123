"""mpmath polynomial helpers: coefficient lists low -> high."""
import itertools
import mpmath as mp


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[k] if k < len(a) else 0) + (b[k] if k < len(b) else 0) for k in range(n)]


def pscale(a, c):
    return [c * x for x in a]


def pmul(a, b):
    out = [mp.mpf(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def pshift(a, c):
    """coefficients of p(u + c)"""
    n = len(a)
    out = [mp.mpf(0)] * n
    for k, ak in enumerate(a):
        if ak == 0:
            continue
        for j in range(k + 1):
            out[j] += ak * mp.binomial(k, j) * c ** (k - j)
    return out


def peval(a, u):
    r = mp.mpf(0)
    for x in reversed(a):
        r = r * u + x
    return r


def ptrim(a, tol=mp.mpf(10) ** -30):
    s = max(abs(x) for x in a) if a else 0
    b = list(a)
    while len(b) > 1 and abs(b[-1]) <= tol * s:
        b.pop()
    return b


def pdivmod(a, b):
    a = list(a)
    q = [mp.mpf(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b):
        c = a[-1] / b[-1]
        k = len(a) - len(b)
        q[k] = c
        for j in range(len(b)):
            a[k + j] -= c * b[j]
        a.pop()
    return q, a


def pdet(M):
    """determinant of a square matrix of polynomials (Leibniz)."""
    n = len(M)
    tot = [mp.mpf(0)]
    for p in itertools.permutations(range(n)):
        s = 1
        for i in range(n):
            for j in range(i + 1, n):
                if p[i] > p[j]:
                    s = -s
        term = [mp.mpf(s)]
        for i in range(n):
            term = pmul(term, M[i][p[i]])
        tot = padd(tot, term)
    return tot


def proots(a):
    a = ptrim(a)
    if len(a) <= 1:
        return []
    return mp.polyroots(list(reversed(a)), maxsteps=200, extraprec=200)


def from_roots(rs, lead=1):
    p = [mp.mpf(lead)]
    for r in rs:
        p = pmul(p, [-r, mp.mpf(1)])
    return p
