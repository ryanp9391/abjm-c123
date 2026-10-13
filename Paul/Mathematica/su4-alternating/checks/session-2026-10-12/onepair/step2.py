# eigenstates, transfer-matrix eigenvalues, Q-functions at the rational point (exact for product states)
import sympy as sp, pickle, time
from chain import *
U = sp.Symbol('U')
hN = sp.Rational(1, 2); th1N, th2N = sp.Rational(2, 7), sp.Rational(-3, 5)
zN = [sp.Rational(3, 2), sp.Rational(5, 4), sp.Rational(2, 3), sp.Rational(4, 7)]
G = companion(esym(zN))
ch = Chain(['A', 'F'], [th1N, th2N], hN, G)
t0 = time.time()
# t_a(U) as full 16x16 operators (bra action: psi*M), polynomial in U
tops = {a: ch.op_t(a, U).applyfunc(sp.expand) for a in range(1, 5)}
print('t ops built', time.time() - t0)
# kets: operators act on kets as M*psi (column), since <phi|M|psi> = (phi M) psi
V = lambda z: sp.Matrix([z**3, z**2, z, 1])
def pvec(a):
    others = [zN[m] for m in range(4) if m != a]
    poly = sp.Poly(sp.prod([U - zz for zz in others]), U)
    return sp.Matrix(poly.all_coeffs())
kets = {}
for a in range(4):
    for b in range(4):
        if a != b:
            kets[(a, b)] = sp.kronecker_product(pvec(a), V(zN[b]))
# check eigen
def eigval(M, psi):
    w = M * psi
    n = next(i for i in range(16) if psi[i] != 0)
    lam = sp.expand(w[n] / psi[n])
    ok = (w - lam * psi).applyfunc(sp.expand) == sp.zeros(16, 1)
    return lam, ok
res = {}
for key, psi in kets.items():
    taus = {}
    for a in range(1, 5):
        lam, ok = eigval(tops[a], psi)
        assert ok, (key, a)
        taus[a] = sp.Poly(lam, U)
    res[key] = taus
print('product states are eigenstates of t1..t4: all ok')
for key in res: print(key, [str(res[key][a].as_expr().factor()) for a in (1, 4)])
pickle.dump({'tops': tops, 'kets': kets, 'res': res}, open('step2.pkl', 'wb'))
