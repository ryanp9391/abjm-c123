import sympy as sp, numpy as np, time
from fractions import Fraction as Fr
from lib import *
t0 = time.time()
ths = sp.symbols('t1 t2 t3 t4')
h = sp.Rational(1, 2)
L = 4
G = np.zeros((4, 4), dtype=object)
for i in range(4):
    for j in range(4):
        G[i, j] = sp.Rational(Z4[i].numerator, Z4[i].denominator) if i == j else sp.Integer(0)
class SC(Chain):
    pass
ch = Chain.__new__(Chain)
ch.n = L; ch.types = ['A', 'F', 'A', 'F']; ch.th = list(ths); ch.h = h; ch.dt = object
ch.G = G; ch.tw = {w: twist_matrix(w, G) for w in (1, 3)}; ch.dim = 4 ** L; ch._lax = {}
ch.conv = lambda x: x
a0 = sp.Rational(1, 3)
def kron(vs):
    out = [sp.Integer(1)]
    for v in vs: out = [x * y for x in out for y in v]
    return np.array(out, dtype=object)
v = [1, 0, 0, a0]; u = [-a0, 0, 0, 1]
psi0 = kron([v, u, v, u])
u1 = ths[0] - h
y = ch.apply(3, u1, psi0)
# compare component with psi0 rho1(Lam^-1): component index of e1 x e4 x e1 x e4 (site1 idx0): position
idx = 0 * 64 + 3 * 16 + 0 * 4 + 3
c3 = sp.factor(sp.together(y[idx] / (psi0[idx] / sp.Rational(Z4[0].numerator, Z4[0].denominator))))
print('c3 =', c3, time.time() - t0)
# check another component ratio equal
idx2 = 3 * 64 + 0 * 16 + 3 * 4 + 0
print('check ratio equal:', sp.simplify(y[idx2] / (psi0[idx2] / sp.Rational(Z4[3].numerator, Z4[3].denominator)) - c3) == 0)
