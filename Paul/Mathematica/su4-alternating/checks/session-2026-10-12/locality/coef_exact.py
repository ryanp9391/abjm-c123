import sys, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
from lib import *
L = int(sys.argv[1]); thp = int(sys.argv[2])
th = {1: TH6, 2: TH6b, 3: TH6c}[thp]
ch = build_chain(L, dt=Fraction); ch.th = th[:L]
h = H; z = Z4; a0 = Fr(1, 3)
def kron_obj(vs):
    out = np.array([Fr(1)], dtype=object)
    for v in vs: out = np.array([x * y for x in out for y in v], dtype=object)
    return out
v = [Fr(1), 0, 0, a0]; u = [-a0, 0, 0, Fr(1)]
psi0 = kron_obj([v, u] * (L // 2))
t0 = time.time()
w3 = ch.apply(3, th[0] - h, psi0); w31 = ch.apply(1, th[1], w3); w1 = ch.apply(1, th[1], psi0)
rest = kron_obj([v, u] * (L // 2 - 1)) if L > 2 else np.array([Fr(1)], dtype=object)
S16 = [Fr(1) if i // 4 == i % 4 else Fr(0) for i in range(16)]
S = np.array([x * y for x in S16 for y in rest], dtype=object)
basis = [psi0, w3, w31, w1, S]
names = ['psi0', 't3(th1-h)', 't3(th1-h)t1(th2)', 't1(th2)', 'S x w0']
print('built words', time.time() - t0, flush=True)
rows = sorted(set(i for b in basis for i in range(len(b)) if b[i] != 0))
A = sp.Matrix([[sp.Rational(b[i].numerator, b[i].denominator) for b in basis] for i in rows])
def solve(vec):
    bvec = sp.Matrix([sp.Rational(Fr(vec[i]).numerator, Fr(vec[i]).denominator) for i in rows])
    sol = (A.T * A).LUsolve(A.T * bvec)
    res = A * sol - bvec
    # also check components outside rows vanish
    outside = [i for i in range(len(vec)) if i not in set(rows) and vec[i] != 0]
    return sol, (res.norm() == 0 and not outside)
# T1(a') = t0 + a' t1 + a'^2 t2 : use a' = 0,1,2 and solve polynomial
def T1(a1):
    return kron_obj([[Fr(1), 0, 0, a1], [-a1, 0, 0, Fr(1)]] + [v, u] * (L // 2 - 1))
sols = [solve(T1(Fr(k))) for k in (0, 1, 2)]
print('T1 solvable (a1 = 0,1,2):', [s[1] for s in sols])
c0, c1, c2 = sols[0][0], (4 * sols[1][0] - 3 * sols[0][0] - sols[2][0]) / 2, (sols[2][0] - 2 * sols[1][0] + sols[0][0]) / 2
ap = sp.symbols('ap')
for k, nm in enumerate(names):
    print(f'  coefficient of {nm:20s} in T1(a\') = {sp.factor(c0[k] + ap * c1[k] + ap**2 * c2[k])}')
# T2 (kappa): numerator polynomial after multiplying by (1-k z1)(1-k z4)
kk = sp.symbols('kappa')
def T2(k):
    vD = [v[i] * (1 - k * z[i]) for i in range(4)]; uD = [u[i] / (1 - k * z[i]) for i in range(4)]
    return kron_obj([vD, uD] + [v, u] * (L // 2 - 1))
ks = [Fr(0), Fr(1, 3), Fr(-1, 2)]
sols2 = [solve(T2(k)) for k in ks]
print('T2 solvable:', [s[1] for s in sols2])
# coefficient c(k) is rational in k: c(k) * (1-k z1)(1-k z4) polynomial degree 2: interpolate with 3 points
def poly(vals):
    out = []
    for idx in range(5):
        pts = [(sp.Rational(k.numerator, k.denominator), vals[j][0][idx] * (1 - sp.Rational(k.numerator, k.denominator) * sp.Rational(z[0].numerator, z[0].denominator)) * (1 - sp.Rational(k.numerator, k.denominator) * sp.Rational(z[3].numerator, z[3].denominator))) for j, k in enumerate(ks)]
        out.append(sp.factor(sp.interpolate(pts, kk)))
    return out
for nm, pl in zip(names, poly(sols2)):
    print(f'  coefficient of {nm:20s} in T2(kappa) x (1-kappa z1)(1-kappa z4) = {pl}')
