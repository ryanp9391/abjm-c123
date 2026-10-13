"""symR6.py -- symbolic (sympy) pair-level words at z = (z1, y, y, z4) and the kappa-family coefficients for the A-frozen references (R5: v = e1; R6: v = (0,1,-1,0)) and for W_a (calibration)."""
import sympy as sp
from fractions import Fraction as Fr
import numpy as np
from rm_base import Chain, diagm, h
z1, y, z4, kap = sp.symbols('z1 y z4 kappa')
hh = sp.Rational(1, 2)
G = np.zeros((4, 4), dtype=object); G[:] = 0
for i, z in enumerate([z1, y, y, z4]): G[i, i] = z
import pgbase
ch = pgbase.Chain(['A', 'F'], [sp.Integer(1), sp.Integer(0)], hh, G)
def ops():
    res = []
    for k, u, reg in [(1, sp.Integer(0), False), (2, hh, False), (3, hh, True)]:
        M = []
        for a in range(4):
            for b in range(4):
                psi = np.zeros((4, 4), dtype=object); psi[:] = sp.Integer(0); psi[a, b] = sp.Integer(1)
                M.append([sp.expand(x) for x in ch.apply_t(k, u, psi, reg=reg).reshape(16)])
        res.append(sp.Matrix(M))
    return res
X1, X2, X3 = ops()
def words(v, u):
    r = sp.Matrix([[a * b for a in v for b in u]])
    return [r, r * X1, r * X2, r * X3, r * X3 * X1]
def target(v, u):
    zz = [z1, y, y, z4]; f = [1 - kap * z for z in zz]
    vv = [v[j] * f[j] for j in range(4)]
    uu = [u[j] * sp.prod([f[m] for m in range(4) if m != j]) for j in range(4)]
    return sp.Matrix([[a * b for a in vv for b in uu]]) / (1 - kap * y) ** 2
def solve(v, u, label):
    ws = words(v, u); T = target(v, u)
    c = sp.symbols('c1:6')
    eqs = [sp.together(sum(c[i] * ws[i][j] for i in range(5)) - T[j]) for j in range(16)]
    eqs = [sp.numer(sp.together(e)) for e in eqs if e != 0]
    sol = sp.solve(eqs, c, dict=True)
    print(label)
    if not sol: print('  no solution'); return None
    for s in sol:
        for k in c:
            val = s.get(k, k)
            print('  ', k, '=', sp.factor(val))
    return sol
u1, u2, u4, a = sp.symbols('u1 u2 u4 a')
# calibration W_a
sW = solve([1, 0, 0, a], [-a, 0, 0, 1], 'W_a')
# K_kappa formula (pair level): s^2, kap s/h^2, kap s/(2 h^4 y), kap s/(6 h^5 y^2), kap^2/(24 h^7 y^2)
s_ = 1 - kap * (z1 + z4)
Kf = [s_**2, kap * s_ / hh**2, kap * s_ / (2 * hh**4 * y), kap * s_ / (6 * hh**5 * y**2), kap**2 / (24 * hh**7 * y**2)]
if sW:
    print('  matches K_kappa formula:', [sp.simplify(sW[0].get(sp.symbols('c%d' % (i + 1)), None) - Kf[i]) == 0 for i in range(5)])
sol6 = solve([0, 1, -1, 0], [u1, u2, u2, u4], 'R6: v = (0,1,-1,0), u = (u1,u2,u2,u4)')
sol5 = solve([1, 0, 0, 0], [0, u2, sp.Symbol('u3'), u4], 'R5: v = e1, u = (0,u2,u3,u4)')
