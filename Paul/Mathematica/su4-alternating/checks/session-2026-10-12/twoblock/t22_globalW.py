"""SoV support of <0| rho(1 + b E14^diag) (global, = HW class of W^L) and two-block (b, b')."""
import numpy as np
from core import *
from tbtools import *
X = ModCtx(); p = X.p; R = X.R
z = z0
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
vz4 = X.vm(e(4), X.VT); l1 = X.vm(e(1), X.Vi); lam1 = X.r(lam[0])
e1, e4 = np.array(e(1)), np.array(e(4))
def pairb(b):
    b = X.r(b); return [(e1 + b * lam1 % p * vz4) % p, (e4 - b * l1) % p]
def show(name, vs):
    c = X.coeffs(X.bra(vs)); r = X.ranks(c)
    c4 = c.reshape(4, 4, 4, 4); used = [sorted(set(np.nonzero(c4)[s])) for s in range(4)]
    print(f'{name:45s}', r, ' slot values used', [list(map(int, u)) for u in used])
show('block1 b=3/7, block2 vac (physical, sparsest)', pairb(Fr(3, 7)) + [e1, e4])
show('block1 vac, block2 b=3/7', [e1, e4] + pairb(Fr(3, 7)))
show('global b=3/7 (W^L HW class)', pairb(Fr(3, 7)) + pairb(Fr(3, 7)))
show('two-block b=3/7, b\'=-2/5', pairb(Fr(3, 7)) + pairb(Fr(-2, 5)))
# companion images: is E14^diag nilpotent rank-1 in companion frame? print the A-site generator
Ec = (np.outer(l1, vz4) * lam1) % p  # v -> v + b (v.l1) lam1? check: v E = (v.l1) vz4 ; for e1: e1.l1 = lam1
print('A-site generator check e1 E = lam1 vz4:', np.array_equal(X.vm(e1, np.outer(l1, vz4) % p), (lam1 * vz4) % p))
