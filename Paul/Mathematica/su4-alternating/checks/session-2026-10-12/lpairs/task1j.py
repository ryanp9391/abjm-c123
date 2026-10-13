# Task 1j: intersection of span(pair-1/pair-2 monomials, per-pair degree <= d) with {W-type states on pairs 1,2} (x) W_a, L = 3
import numpy as np, itertools, time
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
e = np.eye(4, dtype=np.int64)
u = [np.multiply.outer(e[0], e[3]), (np.multiply.outer(e[3], e[3]) - np.multiply.outer(e[0], e[0])) % P, (-np.multiply.outer(e[3], e[0])) % P]
wa = eng.prod_state(eng.W(a), 'p')
U = [np.multiply.outer(np.multiply.outer(u[i], u[j]) % P, wa).reshape(-1) % P for i in range(3) for j in range(3)]
for d in (1, 2):
    PM = [m for m in itertools.product(range(d + 1), repeat=3) if sum(m) <= d]
    mons = [get(m1, m2).reshape(-1) for m1 in PM for m2 in PM]
    cols = mons + [(-x) % P for x in U]
    rk, cons, sol, null = modlin.solve(cols, np.zeros_like(cols[0]))
    print('per-pair degree <= %d: %d monomials; intersection with W-type(pairs 1,2) (x) W_a: dim %d' % (d, len(mons), len(null)))
    for nv in null: print('    d_ij:', [str(modlin.ratrec(int(x))) for x in nv[len(mons):]])
