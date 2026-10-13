import sympy as sp, numpy as np, time
from lib import *
t0 = time.time()
ths = sp.symbols('t1 t2 t3 t4'); h = sp.Rational(1, 2); L = 4
za = sp.symbols('z1 z2 z3 z4')
a = sp.symbols('a')
G = np.zeros((4, 4), dtype=object)
for i in range(4):
    for j in range(4):
        G[i, j] = za[i] if i == j else sp.Integer(0)
ch = Chain.__new__(Chain)
ch.n = L; ch.types = ['A', 'F', 'A', 'F']; ch.th = list(ths); ch.h = h; ch.dt = object
ch.G = G; ch.tw = {w: twist_matrix(w, G) for w in (1,)}; ch.dim = 4 ** L; ch._lax = {}
ch.conv = lambda x: x
def kron(vs):
    out = [sp.Integer(1)]
    for v in vs: out = [x * y for x in out for y in v]
    return np.array(out, dtype=object)
v = [1, 0, 0, a]; u = [-a, 0, 0, 1]
psi0 = kron([v, u, v, u])
y = ch.apply(1, ths[1], psi0)
# pair-1 x pair-2 reshape: components
M = y.reshape(16, 16)
w0 = kron([v, u])
# pick components: pair1 (A=1,F=4) index 3 ; (A=2,F=2) index 5 ; with pair2 index of (1,4)=3
al = sp.factor(sp.together(M[3, 3] / (w0[3] * za[3])))      # alpha: coefficient of v (x) u Lam at (1,4): v1*u4*z4 = z4
sg = sp.factor(sp.together(M[5, 3] / w0[3]))                  # sigma: (2,2) entry
print('alpha =', al); print('sigma =', sg); print('sigma/alpha =', sp.factor(sp.together(sg / al)), time.time() - t0)
