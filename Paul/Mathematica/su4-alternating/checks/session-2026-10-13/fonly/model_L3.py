import itertools, sympy as sp
from model_cherednik import model
h = sp.Rational(1, 2); p, s = sp.symbols('p s')
def coeffs(th):
    L = 3
    basis, idx, Y = model(L, th, h); n = len(basis)
    vac = sp.zeros(1, n); vac[0, 0] = 1
    ks = list(itertools.product(range(3), repeat=L))
    rows = []
    for k in ks:
        v = vac
        for si, kk in enumerate(k):
            for _ in range(kk): v = v * Y[si]
        rows.append(v)
    A = sp.Matrix.vstack(*rows)
    w = [p, -s, 1]
    tgt = sp.Matrix([[sp.prod([w[l] for l in t]) for t in basis]])
    c = (A.T).LUsolve(tgt.T)
    return ks, [sp.expand(ci) for ci in c], A
th = [sp.Rational(1, 3), sp.Rational(-2, 5), sp.Rational(3, 7)]
ks, c, A = coeffs(th)
cd = dict(zip(ks, c))
# symmetry in k
bad = 0
for k in ks:
    for perm in itertools.permutations(range(3)):
        kk = tuple(k[i] for i in perm)
        if sp.expand(cd[k] - cd[kk]) != 0: bad += 1
print('number of asymmetric (k, perm) pairs:', bad)
# express in monomial symmetric polynomials of the Y's: list distinct coefficient classes
cls = {}
for k in ks:
    key = tuple(sorted(k, reverse=True)); cls.setdefault(key, set()).add(sp.srepr(cd[k]))
print({k: len(v) for k, v in cls.items()})
import pickle; pickle.dump({str(k): str(v) for k, v in cd.items()}, open('L3_coeffs_th.pkl', 'wb'))
for key in sorted(cls):
    print(key, sp.factor(cd[key]))
