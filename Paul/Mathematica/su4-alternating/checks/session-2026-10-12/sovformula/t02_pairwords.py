"""Task 1b: the 16 regularised pair-1 SoV words on <0| (companion, ABJM point pair-by-pair, thF1 != thF2):
   O1 t1(thF1)^k, O1 in {1, t1(u0), W = t2'(u0) - t1'(u0) t1(thF), Z = t3'(u0)}, u0 = thF1 + h.
   Which combinations localise on pair 1, i.e. lie in C^16 (x) (e1 e4)^{L-1}?  exact Fractions."""
import numpy as np, sys, itertools, time
from frlin import nullspace
from pgbase import *
from dual import D, val, der
L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
h = Fr(1, 2)
thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)][:L]
th = []
for f in thF: th += [f + 2 * h, f]
zgen = len(sys.argv) > 2 and sys.argv[2] == 'gen'
z = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)] if zgen else [Fr(3, 2), Fr(5, 4), Fr(5, 4), Fr(4, 7)]
G = companion(z)
ch = Chain(['A', 'F'] * L, th, h, G)
e = lambda i: [Fr(1) if k == i else Fr(0) for k in range(4)]
vac = prod_state([e(0), e(3)] * L)
u0, uF = th[1] + h, th[1]
def tder(k, psi):
    """(psi t_k(u0), psi t_k'(u0)) via dual numbers"""
    pd = np.vectorize(lambda x: D(x, 0), otypes=[object])(psi)
    r = ch.apply_t(k, D(u0, Fr(1)), pd)
    return np.vectorize(val, otypes=[object])(r), np.vectorize(der, otypes=[object])(r)
def O(name, psi):
    if name == '1': return psi
    if name == 't1': return ch.apply_t(1, u0, psi)
    if name == 'W':
        _, d2 = tder(2, psi); _, d1 = tder(1, psi)
        return d2 - ch.apply_t(1, uF, d1)        # psi t2'(u0) - psi t1'(u0) t1(thF)
    if name == 'Z':
        v3, d3 = tder(3, psi); assert all(x == 0 for x in v3.reshape(-1)); return d3
names = ['1', 't1', 'W', 'Z']
words = {}
T0 = time.time()
for n in names:
    print("word", n, flush=True); s = O(n, vac)
    for k in range(4):
        words[(n, k)] = s
        s = ch.apply_t(1, uF, s)
keys = list(words)
print('words computed in %.1f s' % (time.time() - T0), flush=True)
M = np.array([words[k].reshape(16, -1) for k in keys])        # (16 words, 16 pair1, rest)
rest = prod_state([e(0), e(3)] * (L - 1)).reshape(-1)
ridx = int(np.nonzero(rest)[0][0])
# component outside pair1 (x) vac_rest: all rest-columns except ridx must vanish
Mout = np.array([words[k].reshape(16, -1)[:, [c for c in range(len(rest)) if c != ridx]].reshape(-1) for k in keys]).T
ns = nullspace(Mout.tolist(), 16)
print('L =', L, 'z =', z, ': #pair-1 SoV words = 16; localising combinations: dim', len(ns))
Pin = np.array([words[k].reshape(16, -1)[:, ridx] for k in keys])   # (16, 16) pair-1 parts
for v in ns:
    comb = {keys[i]: v[i] for i in range(16) if v[i] != 0}
    st = sum((v[i] * Pin[i] for i in range(16) if v[i] != 0), np.zeros(16, dtype=object))
    st = st.reshape(4, 4)
    print('  combo', {f'{a}*t1F^{b}': str(c) for (a, b), c in comb.items()})
    print('     pair-1 state (A,F,coef):', [(i + 1, j + 1, str(st[i, j])) for i in range(4) for j in range(4) if st[i, j] != 0])
print('rank of pair-1 parts of all 16 words:', 16 - len(nullspace(Pin.T.tolist(), 16)))
