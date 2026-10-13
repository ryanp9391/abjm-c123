"""CRT + rational reconstruction over several primes."""
import numpy as np, pickle, math
from fractions import Fraction as Fr
_res8 = None
def res8():
    global _res8
    if _res8 is None: _res8 = pickle.load(open('dual_mod8.pkl', 'rb'))
    return _res8

def crt_pair(a1, m1, a2, m2):
    return (a1 + m1 * ((a2 - a1) * pow(m1, -1, m2) % m2)) % (m1 * m2), m1 * m2

def ratrec_big(a, m):
    a %= m
    r0, r1 = m, a; s0, s1 = 0, 1
    B = math.isqrt(m // 2)
    while r1 > B:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)

def reconstruct(vals_by_prime, nprimes=None):
    """vals_by_prime: dict p -> int residue (scalar). Uses all but last prime for reconstruction, last one to verify."""
    ps = list(vals_by_prime.keys())
    if nprimes: ps = ps[:nprimes]
    a, m = vals_by_prime[ps[0]] % ps[0], ps[0]
    for p in ps[1:-1]:
        a, m = crt_pair(a, m, vals_by_prime[p] % p, p)
    r = ratrec_big(a, m)
    if r is None: return None
    pl = ps[-1]
    if (r.numerator - r.denominator * vals_by_prime[pl]) % pl != 0: return None
    return r
