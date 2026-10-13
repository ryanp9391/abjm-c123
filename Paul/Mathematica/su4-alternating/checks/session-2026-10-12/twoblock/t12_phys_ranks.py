import numpy as np, pickle
from fractions import Fraction as Fr
import sympy
from core import *
from exact import *
R8 = res8(); ps = list(R8.keys()); p = ps[0]
# recompute mod-p coefficient tensors for the gamma-pieces (avoid reconstruction failures)
z = z0
Vs = sympy.Matrix([[sympy.Rational(x.numerator, x.denominator) for x in row] for row in Vand(z)]); Vi = Vs.inv()
lam = [Fr(1) / np.prod([z[j] - z[m] for m in range(4) if m != j]) for j in range(4)]
def tofr(x): return Fr(int(sympy.fraction(x)[0]), int(sympy.fraction(x)[1]))
VT = [[tofr(Vs[j, i]) for j in range(4)] for i in range(4)]
Vinv = [[tofr(Vi[i, j]) for j in range(4)] for i in range(4)]
def vm(v, A): return [sum(Fr(v[i]) * A[i][j] for i in range(4)) for j in range(4)]
A2c = lambda v: vm(v, VT); F2c = lambda u: vm(u, Vinv)
one = [1, 1, 1, 1]; e1 = [1, 0, 0, 0]; e4 = [0, 0, 0, 1]
R = Ring('mod', p); Y = R8[p][1]
def cvec(vs):
    b = kronvecs(R, [R.arr(v) for v in vs]); return R.mm(b[None, :], Y)[0]
pair2 = [A2c(lam), F2c(one)]
pieces = {
  'g0': cvec([A2c(lam), F2c(one)] + pair2),
  'g1': (cvec([A2c(lam), F2c([1 / lam[0], 0, 0, 0])] + pair2) - cvec([A2c(e4), F2c(one)] + pair2)) % p,
  'g2': (-cvec([A2c(e4), F2c([1 / lam[0], 0, 0, 0])] + pair2)) % p,
}
for k, c in pieces.items():
    c4 = c.reshape(4, 4, 4, 4)
    print(k, 'nnz', np.count_nonzero(c), ' a3 values used', sorted(set(np.nonzero(c4)[2])), ' k4 values', sorted(set(np.nonzero(c4)[3])))
    for k4 in range(3):
        F = c4[:, :, 0, k4]
        print('   k4=%d  F (a1 x k2) rank %d  support:' % (k4, modrank(F, p)), [(a, b) for a in range(4) for b in range(4) if F[a, b]])
    print('   pair-cut rank', modrank(c4.reshape(16, 16), p), ' site ranks', [modrank(np.moveaxis(c4, s, 0).reshape(4, 64), p) for s in range(4)])
pickle.dump(pieces, open('phys_pieces_modp.pkl', 'wb'))
