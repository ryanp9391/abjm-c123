"""Exact rational matrices from several GF(p) evaluations (CRT + rational reconstruction)."""
import numpy as np
from fractions import Fraction as Fr
from math import isqrt
import sympy

PRIMES = []
_p = 16777213
while len(PRIMES) < 80:
    PRIMES.append(_p)
    _p = sympy.prevprime(_p)


def ratrec_big(a, m):
    a %= m
    bound = isqrt(m // 2)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)


def exact_from_primes(fn, shape_hint=None, max_primes=60, min_primes=2, verbose=False):
    """fn(p) -> int64 array mod p. Returns object array of Fractions."""
    X = None
    m = 1
    prev = None
    for k, p in enumerate(PRIMES[:max_primes]):
        a = np.asarray(fn(p), dtype=object) % p
        if X is None:
            X = a.copy()
            m = p
        else:
            inv = pow(m % p, -1, p)
            t = ((a - X % p) * inv) % p
            X = X + m * t
            m = m * p
        if k + 1 >= min_primes:
            flat = X.ravel()
            rec = [ratrec_big(int(x), m) for x in flat]
            if all(r is not None for r in rec):
                cur = np.array(rec, dtype=object).reshape(X.shape)
                if prev is not None and np.all(cur == prev):
                    if verbose:
                        print('  stable after', k + 1, 'primes')
                    return cur
                prev = cur
    raise RuntimeError('no stable reconstruction')
