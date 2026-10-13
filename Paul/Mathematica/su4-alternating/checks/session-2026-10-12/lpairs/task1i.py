# Task 1i: higher per-pair degree.  Monomials prod_p X1_p^i X2_p^j X3_p^k with per-pair degree <= d, pairs 1..l, at L = 3, l = 2.
# Is <W_a'^2 (x) W_a| in the span?  (also: is it in the span of pair-1-only monomials of degree <= D?)
import numpy as np, itertools, time, sys
from fractions import Fraction as Fr
import eng, lp, modlin
from eng import P
h = Fr(1, 2); y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
Lp = 3; thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)]
ch = lp.chain_abjm(Lp, thF, h, z, 'p')
a = Fr(1, 3); ref = lp.Wstate(ch, [a] * Lp)
def op(p, k, v):
    uF = ch.th[2 * p - 1]
    if k == 1: return ch.apply_t(1, uF, v)
    if k == 2: return ch.apply_t(2, uF + h, v)
    return ch.apply_t(3, uF + h, v, reg=True)
def pair_monos(d):
    return [m for m in itertools.product(range(d + 1), repeat=3) if sum(m) <= d]
kap = Fr(3, 10); ap = a * (1 - kap * z[3]) / (1 - kap * z[0])
tgt2 = lp.Wstate(ch, [ap, ap, a])
for d in (1, 2, 3):
    t0 = time.time()
    # build all monomials pair1^m1 pair2^m2 recursively with caching
    cache = {}
    def get(m1, m2):
        key = (m1, m2)
        if key in cache: return cache[key]
        if sum(m2) > 0:
            k = [i for i in range(3) if m2[i] > 0][0]; mm = list(m2); mm[k] -= 1
            v = op(2, k + 1, get(m1, tuple(mm)))
        elif sum(m1) > 0:
            k = [i for i in range(3) if m1[i] > 0][0]; mm = list(m1); mm[k] -= 1
            v = op(1, k + 1, get(tuple(mm), m2))
        else: v = ref
        cache[key] = v; return v
    PM = pair_monos(d)
    mons = [get(m1, m2) for m1 in PM for m2 in PM]
    rk, cons, sol, null = modlin.solve(mons, tgt2)
    print('d = %d: %d monomials, rank %d, target W_a\'^2 W_a in span: %s (%.0f s)' % (d, len(mons), rk, cons, time.time() - t0))
    mons1 = [get(m1, (0, 0, 0)) for m1 in pair_monos(2 * d)]
    rk1, cons1, _, _ = modlin.solve(mons1, tgt2)
    print('        pair-1 only, degree <= %d: %d monomials, rank %d, target in span: %s' % (2 * d, len(mons1), rk1, cons1))
