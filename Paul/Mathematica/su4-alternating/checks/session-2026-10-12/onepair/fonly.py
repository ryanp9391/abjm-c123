# F-only degree-1 rotations of both pairs: exact symbolic expansion in the 4-site SoV words (exchange term)
import sympy as sp, numpy as np
from itertools import combinations
import fastchain as fc
th = sp.symbols('theta1:5'); h = sp.Symbol('h'); cs = sp.symbols('c1:5')
class SymRing:
    mode = 'sym'
    def c(self, x): return sp.sympify(x)
    def red(self, a):
        if isinstance(a, np.ndarray): return np.vectorize(sp.expand, otypes=[object])(a)
        return sp.expand(a)
    def zeros(self, n): return np.full(n, sp.Integer(0), dtype=object)
R = SymRing()
ch = fc.FastChain.__new__(fc.FastChain)
ch.types = ['A', 'F', 'A', 'F']; ch.L = 4; ch.R = R; ch.th = list(th); ch.h = h; ch.dim = 256
ch.G = [[cs[0], -cs[1], cs[2], -cs[3]], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]]
def local(psi, s, a, b, u):
    t = psi.reshape((4,) * 4)
    out = np.full(t.shape, sp.Integer(0), dtype=object)
    if a == b: out = t * (u - th[s])
    idx_src = [slice(None)] * 4; idx_dst = [slice(None)] * 4
    if ch.types[s] == 'A': idx_src[s] = a; idx_dst[s] = b; f = h
    else: idx_src[s] = b; idx_dst[s] = a; f = -h
    out[tuple(idx_dst)] = out[tuple(idx_dst)] + t[tuple(idx_src)] * f
    return out.reshape(-1)
ch.local = local
def t1(psi, u):
    TG = fc.FastChain.braTG(ch, psi, u)
    return R.red(sum(TG[i][i] for i in range(4)))
idx = lambda i1, j1, i2, j2: ((i1 * 4 + j1) * 4 + i2) * 4 + j2
vac = R.zeros(256); vac[idx(0, 3, 0, 3)] = 1
w = {}
w[(0, 0)] = vac
w[(1, 0)] = t1(vac, th[1])
w[(0, 1)] = t1(vac, th[3])
w[(1, 1)] = t1(w[(1, 0)], th[3])
for k, v in w.items():
    nz = {(i // 64, (i // 16) % 4, (i // 4) % 4, i % 4): sp.factor(v[i]) for i in range(256) if v[i] != 0}
    print('word (k2,k4) =', k, ':', {tuple(x + 1 for x in key): val for key, val in nz.items()})
g2, d2, g4, d4 = sp.symbols('gamma2 delta2 gamma4 delta4')
e1 = [1, 0, 0, 0]; uu2 = [0, 0, d2, g2]; uu4 = [0, 0, d4, g4]
bra = np.array([e1[i1] * uu2[j1] * e1[i2] * uu4[j2] for i1 in range(4) for j1 in range(4) for i2 in range(4) for j2 in range(4)], dtype=object)
x = sp.symbols('x00 x10 x01 x11')
eqs = [sp.expand(x[0] * w[(0, 0)][i] + x[1] * w[(1, 0)][i] + x[2] * w[(0, 1)][i] + x[3] * w[(1, 1)][i] - bra[i]) for i in range(256)]
sol = sp.solve([e for e in eqs if e != 0], x, dict=True)[0]
t = {(i, j): th[i - 1] - th[j - 1] for i in range(1, 5) for j in range(1, 5)}
Y2 = h * t[1, 2] * t[2, 3] * (t[2, 4] - h)
Y24 = h**2 * t[1, 2] * t[1, 4] * t[2, 3] * t[3, 4] * (t[2, 4] - h) * (t[2, 4] + h)
pred = {x[0]: g2 * g4, x[1]: (d2 * g4 - h * g2 * d4 / t[2, 4]) / Y2, x[2]: g2 * d4 / (h * t[1, 4] * t[2, 4] * t[3, 4]), x[3]: d2 * d4 / Y24}
print('F-only exchange formula coefficients match:', all(sp.simplify(sol[k] - pred[k]) == 0 for k in x))
