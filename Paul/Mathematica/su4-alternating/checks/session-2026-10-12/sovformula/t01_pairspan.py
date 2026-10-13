"""Task 1a: localised words on the companion SoV vacuum <0| = (e1 e4)^L, ABJM point (thA_n = thF_n + 2h), exact Fractions.
Checks <0| X = lambda_X * (pair-1 state) (x) (e1 e4)^{L-1} and prints the pair-1 states."""
import numpy as np, sys
from pgbase import *
L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
h = Fr(1, 2)
thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)][:L]
th = []
for f in thF: th += [f + 2 * h, f]
y = Fr(5, 4); z = [Fr(3, 2), y, y, Fr(4, 7)]
if len(sys.argv) > 2 and sys.argv[2] == 'gen': z = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]
G = companion(z)
ch = Chain(['A', 'F'] * L, th, h, G)
def lam(k, u0):
    r = Fr(1)
    for s in range(2, 2 * L):
        t = th[s]
        if s % 2 == 0:
            for j in range(k): r *= (u0 - t - j * h)
        else:
            r *= (u0 - t - k * h)
            for j in range(k - 1): r *= (u0 - t - j * h)
    return r
e = lambda i: [Fr(1) if k == i else Fr(0) for k in range(4)]
vac = prod_state([e(0), e(3)] * L)
u1, u2 = th[1], th[1] + h
X = {'1': vac}
X['X1'] = ch.apply_t(1, u1, vac)
X['X2'] = ch.apply_t(2, u2, vac)
X['X3'] = ch.apply_t(3, u2, vac, reg=True)
X['X3X1'] = ch.apply_t(1, u1, X['X3'])
X['t1(u0)'] = ch.apply_t(1, u2, vac)
X['t1(u0)X1'] = ch.apply_t(1, u2, X['X1'])
lamv = {'1': 1, 'X1': lam(1, u1), 'X2': lam(2, u2), 'X3': lam(3, u2), 'X3X1': lam(3, u2) * lam(1, u1)}
rest = prod_state([e(0), e(3)] * (L - 1))
ch1 = Chain(['A', 'F'], th[:2], h, G)
v1 = prod_state([e(0), e(3)])
P = {'1': v1, 'X1': ch1.apply_t(1, u1, v1), 'X2': ch1.apply_t(2, u2, v1), 'X3': ch1.apply_t(3, u2, v1, reg=True)}
P['X3X1'] = ch1.apply_t(1, u1, P['X3'])
for k in ['1', 'X1', 'X2', 'X3', 'X3X1']:
    pred = lamv[k] * np.multiply.outer(P[k], rest)
    ok = all(x == 0 for x in (X[k] - pred).reshape(-1))
    nzp = [(i, j, P[k][i, j]) for i in range(4) for j in range(4) if P[k][i, j] != 0]
    print(f'{k:5s} localises as lambda * pair word: {ok};  pair state (A idx, F idx, coeff):', [(i + 1, j + 1, str(c)) for i, j, c in nzp])
# X2 = t1(u0) X1 identity
print('X2 == t1(u0) X1 on <0| :', all(x == 0 for x in (X['X2'] - X['t1(u0)X1']).reshape(-1)))
# does t1(u0) alone localise on pair 1?
r = X['t1(u0)']; flat = r.reshape(16, -1)
print('t1(u0) on <0|: rank of pair1 x rest matrix =', np.linalg.matrix_rank(np.array(flat, dtype=float)))
