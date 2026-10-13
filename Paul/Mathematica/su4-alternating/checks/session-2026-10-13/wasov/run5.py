"""Vacuum Q-functions of twist zz (exact, sympy), bra-side dictionary, mirror (dual) dictionary with state/twist independent constants,
and its use for tauTilde_y(vac_w).  Compare with notebook cells 16-19 of Wa-sov-expansion-L4.wb."""
import itertools, sympy as sp
from sympy import Rational as R
from fractions import Fraction as Fr
h = R(1, 2); th = [R(2, 7), R(-3, 5), R(8, 11), R(-11, 8)]
types = ['A', 'F', 'A', 'F']
Z = [R(3, 2), R(5, 4), R(2, 3), R(4, 7)]; W = [R(7, 5), R(9, 8), R(3, 7), R(5, 9)]; Z2 = [R(11, 7), R(5, 3), R(2, 9), R(13, 5)]
u = sp.symbols('u')
mu = [2, 0, 0, -2]
def V(i, x):   # vacuum eigenvalue of T_ii on the ket vacuum (i = 1..4)
    if i == 1: return (x - th[0] + h) * (x - th[1]) * (x - th[2] + h) * (x - th[3])
    if i in (2, 3): return (x - th[0]) * (x - th[1]) * (x - th[2]) * (x - th[3])
    return (x - th[0]) * (x - th[1] - h) * (x - th[2]) * (x - th[3] - h)
def tauvac(zz, a, x):
    if a == 0: return sp.Integer(1)
    return sum(sp.prod([zz[I[k]] * V(I[k] + 1, x - k * h) for k in range(a)]) for I in itertools.combinations(range(4), a))
def vacq(zz, i):
    d = 2 + mu[i - 1]
    cols = []
    for k in range(d + 1):
        expr = sum((-1) ** a * zz[i - 1] ** (-a) * tauvac(zz, a, u) * (u - a * h) ** k * sp.prod([V(4, u - b * h) for b in range(a, 4)]) for a in range(5))
        cols.append(sp.Poly(sp.expand(expr), u).all_coeffs()[::-1])
    n = max(len(c) for c in cols)
    Mx = sp.Matrix([[c[r] if r < len(c) else 0 for c in cols] for r in range(n)])
    ns = Mx.nullspace()
    assert len(ns) == 1, len(ns)
    v = ns[0] / ns[0][-1]
    return sum(v[k] * u ** k for k in range(d + 1))
def Fvac(k):
    F = [th[1], th[3]]; A = [th[0], th[2]]
    if k == 1: return sp.Integer(1)
    if k == 2: return sp.prod([u - t for t in F])
    if k == 3: return sp.prod([(u - t) * (u - t - h) for t in F])
    return sp.prod([(u - t) * (u - t - h) * (u - t - 2 * h) for t in F]) * sp.prod([u - t + h for t in A])
def vacQ(zz):
    q = {(i,): vacq(zz, i) for i in range(1, 5)}
    for k in (2, 3, 4):
        for I in itertools.combinations(range(1, 5), k):
            P = sp.expand(sp.Matrix([[zz[i - 1] ** (k - 1 - s) * q[(i,)].subs(u, u - s * h) for s in range(k)] for i in I]).det())
            qq, rr = sp.div(sp.Poly(P, u), sp.Poly(Fvac(k), u))
            assert rr.is_zero, (I, rr)
            q[I] = qq.as_expr()
    return q
triples = list(itertools.combinations(range(1, 5), 3))
def wd(zz, q, I, x, t):
    return sp.prod([zz[i - 1] for i in I]) ** ((x - t) / h) * q[I].subs(u, x)
def dualconsts(zz):
    q = vacQ(zz); out = {}
    for s in range(4):
        t = th[s]
        if types[s] == 'A':
            out[(s, 1)] = [sp.nsimplify(tauvac(zz, 1, t) * wd(zz, q, I, t - h, t) / wd(zz, q, I, t - 2 * h, t)) for I in triples]
        else:
            mm = sp.Matrix([[wd(zz, q, triples[l], t + (j - 1) * h, t) for j in range(4)] for l in range(3)])
            D0 = mm[:, [1, 2, 3]].det()
            for n in (1, 2, 3):
                Dn = mm[:, [j for j in range(4) if j != n]].det()
                out[(s, n)] = [sp.nsimplify(tauvac(zz, n, t + n * h) * D0 / Dn)]
    return out, q
cz, qz = dualconsts(Z); cw, qw = dualconsts(W); c2, q2 = dualconsts(Z2)
print('A-site constants independent of the triple:', all(len(set(v)) == 1 for k, v in cz.items() if types[k[0]] == 'A'))
print('constants at twist z, w, w2 equal:', all(cz[k][0] == cw[k][0] == c2[k][0] for k in cz))
for k in sorted(cz): print(' const', k, cz[k][0])
# bra-side dictionary on the vacuum Q's of twist W
def Qt(zz, q, i, t, j): return zz[i - 1] ** (-j) * sp.prod([V(4, t - b * h) for b in range(j, 5)]) * q[(i,)].subs(u, t - j * h)
def braerr(zz, q):
    e = 0
    for s in (1, 3):
        for lab in range(4):
            for j in range(1, 5):
                e = max(e, abs(sp.nsimplify((Qt(zz, q, j, th[s], 0) / Qt(zz, q, j, th[s], 1)) ** lab / tauvac(zz, 1, th[s]) ** lab - 1)))
    for s in (0, 2):
        for lab in range(4):
            for tri in triples:
                D = lambda pts: sp.Matrix([[Qt(zz, q, i, th[s], jj) for jj in pts] for i in tri]).det()
                val = D([jj for jj in (1, 2, 3, 4) if jj != lab + 1]) / D([2, 3, 4])
                e = max(e, abs(sp.nsimplify(val / tauvac(zz, lab, th[s] - h) - 1)))
    return e
print('bra-side dictionary error (twist w):', braerr(W, qw), ' (twist z):', braerr(Z, qz))
# mirror dictionary error: constants from z used at w, w2
def dualsite(zz, q, s, n):
    t = th[s]
    if n == 0: return sp.Integer(1)
    if types[s] == 'A':
        return (cz[(s, 1)][0] / sp.prod(Z) * sp.prod(zz) * wd(zz, q, triples[0], t - 2 * h, t) / wd(zz, q, triples[0], t - h, t)) ** n
    mm = sp.Matrix([[wd(zz, q, triples[l], t + (j - 1) * h, t) for j in range(4)] for l in range(3)])
    return cz[(s, n)][0] * (sp.prod(zz) / sp.prod(Z)) ** n * mm[:, [j for j in range(4) if j != n]].det() / mm[:, [1, 2, 3]].det()
def tauy_Q(zz, q, y):
    return sp.prod([dualsite(zz, q, s, y[s]) for s in range(4)])
def tauy_ev(zz, y):
    return sp.prod([(1 if y[s] == 0 else (tauvac(zz, 1, th[s]) ** y[s] if types[s] == 'A' else tauvac(zz, y[s], th[s] + y[s] * h))) for s in range(4)])
errs = []
for zz, q in ((W, qw), (Z2, q2), (Z, qz)):
    errs.append(max(abs(sp.nsimplify(tauy_Q(zz, q, y) / tauy_ev(zz, y) - 1)) for y in itertools.product(range(4), repeat=4)))
print('mirror dictionary (constants from vac(z)) vs eigenvalue product, all 256 y, at w, w2, z:', errs)
print('--- ratios of constants (w or w2 over z) and powers of prod(z):')
pz = sp.prod(Z); pw = sp.prod(W); p2 = sp.prod(Z2)
for k in sorted(cz):
    rw = sp.nsimplify(cw[k][0] / cz[k][0]); r2 = sp.nsimplify(c2[k][0] / cz[k][0])
    # find integer m with rw == (pw/pz)^m
    found = [m for m in range(-6, 7) if sp.simplify(rw - (pw / pz) ** m) == 0 and sp.simplify(r2 - (p2 / pz) ** m) == 0]
    print(k, rw, r2, 'power of detG ratio:', found)
