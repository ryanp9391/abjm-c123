"""Joint eigenstates of t1, t2, t3 (diagonal frame, weight blocks, mpmath 50 digits) and their eigenvalue polynomials.
Output: eigen.pkl with per state: weight, block indices, eigenvector (diag frame, mp), reduced tau-bar coefficient lists."""
import time, itertools, pickle, sys
import numpy as np
import mpmath as mp
from fractions import Fraction as Fr
from sov_l4 import *

mp.mp.dps = 50
h, th, z = H0, TH0, Z0
chd = Chain(G=Gdiag(z))  # exact, diagonal twist


def weight(idx):
    i1, j2, i3, j4 = idx
    w = [0] * 4
    w[i1] += 1; w[i3] += 1; w[j2] -= 1; w[j4] -= 1
    return tuple(w)


basis = list(itertools.product(range(4), repeat=4))
blocks = {}
for n, idx in enumerate(basis):
    blocks.setdefault(weight(idx), []).append(n)
print('# weight blocks', len(blocks), 'max size', max(len(b) for b in blocks.values()))

# interpolation points for the reduced tau-bar_a (degree 4)
upts = [Fr(k, 3) + Fr(1, 17) for k in range(-2, 3)]
t0 = time.time()
mats = {(a, u): chd.t(a, u, reduced=True) for a in (1, 2, 3) for u in upts}
# combination for diagonalisation
u0, u1, u2 = Fr(3, 11), Fr(-5, 7), Fr(2, 9)
Mcomb = chd.t(1, u0) + Fr(7, 10) * chd.t(3, u1, reduced=True) + Fr(-3, 13) * chd.t(2, u2, reduced=True)
print('exact matrices built', time.time() - t0)


def tomp(x):
    return mp.mpf(x.numerator) / x.denominator


def block_mp(M, idx):
    return mp.matrix([[tomp(M[r, c]) for c in idx] for r in idx])


# check block structure (exact): t's preserve weight
for (a, u), M in list(mats.items())[:3]:
    off = 0
    for w, idx in blocks.items():
        others = [n for n in range(DIM) if n not in set(idx)]
        sub = M[np.ix_(idx, others)]
        off += sum(1 for x in sub.flat if x != 0)
    print('off-block nonzeros', a, u, off)

states = []
t0 = time.time()
for w, idx in sorted(blocks.items()):
    B = block_mp(Mcomb, idx)
    E, ER = mp.eig(B)
    n = len(idx)
    # sanity: simple spectrum within block
    for k in range(n):
        v = ER[:, k]
        # normalise: largest component 1
        p = max(range(n), key=lambda r: abs(v[r]))
        v = v / v[p]
        taus = {}
        res = mp.mpf(0)
        for a in (1, 2, 3):
            vals = []
            for u in upts:
                Mb = block_mp(mats[(a, u)], idx)
                tv = Mb * v
                lam = tv[p] / v[p]
                res = max(res, mp.norm(tv - lam * v))
                vals.append(lam)
            # interpolate degree 4 polynomial: coefficients c0..c4
            Vm = mp.matrix([[tomp(u) ** j for j in range(5)] for u in upts])
            coef = mp.lu_solve(Vm, mp.matrix(vals))
            taus[a] = [coef[j] for j in range(5)]
        states.append({'weight': w, 'idx': idx, 'vec': [v[r] for r in range(n)], 'taubar': taus, 'res': res, 'combEig': E[k]})
    print(w, n, 'maxres %.2e' % float(max(s['res'] for s in states[-n:])), 't=%.1f' % (time.time() - t0), flush=True)

# global simplicity of the joint spectrum
ev = [s['combEig'] for s in states]
gap = min(abs(ev[i] - ev[j]) for i in range(len(ev)) for j in range(i))
print('# states', len(states), 'min gap of combination spectrum', mp.nstr(gap, 5))
pickle.dump({'states': states, 'upts': upts, 'blocks': blocks}, open('eigen.pkl', 'wb'))
