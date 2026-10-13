"""Exact bilinear forms K(g_s, g_s') = <B|Psi>/<0|Psi> for two rotated sites (rest at SoV vacuum), as polynomials in the per-site
SoV values; tests of 'one row per site' determinant structure (det of the 2x2 form must be separable)."""
import numpy as np, pickle, sympy as sp
from fractions import Fraction as Fr
from core import *
from exact import crt_pair, ratrec_big
Ys = pickle.load(open('Y_mod20.pkl', 'rb')); primes = list(Ys.keys())
th = th0; h = h0
def idx(c): return ((c[0] * 4 + c[1]) * 4 + c[2]) * 4 + c[3]
def coeffs_of_component(comp):
    """exact c_x for the product basis bra e_comp (0-based component tuple); returns dict label->Fraction"""
    row = idx(comp); out = {}
    for x in range(256):
        vals = [int(Ys[pp][row, x]) for pp in primes]
        if all(v == 0 for v in vals[:3]): continue
        a_, m_ = vals[0], primes[0]
        for v, pp in zip(vals[1:-1], primes[1:-1]): a_, m_ = crt_pair(a_, m_, v % pp, pp)
        r = ratrec_big(a_, m_)
        assert r is not None and (r.numerator - r.denominator * vals[-1]) % primes[-1] == 0, 'reconstruction'
        out[wordsIdx[x]] = r
    return out
# per-site SoV symbols
A1 = sp.symbols('t0_1 t1_1 t2_1 t3_1'); A3 = sp.symbols('t0_3 t1_3 t2_3 t3_3'); r2, r4 = sp.symbols('r2 r4')
def Wsym(x):
    return (1 if x[0] == 0 else A1[x[0]]) * r2 ** x[1] * (1 if x[2] == 0 else A3[x[2]]) * r4 ** x[3]
def Kpoly(comp):
    return sp.expand(sum(sp.Rational(v.numerator, v.denominator) * Wsym(x) for x, v in coeffs_of_component(comp).items()))
S = lambda f: sp.Rational(f.numerator, f.denominator)
th12, th13, th14, th23, th24, th34 = [S(th[i] - th[j]) for i, j in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]]
hS = S(h)
vac = (0, 3, 0, 3)
def form(siteA, compsA, siteB, compsB):
    """2-site bilinear form: matrix K[i][j] for components compsA at siteA, compsB at siteB (0-based), rest vacuum"""
    Kmat = []
    for ca in compsA:
        row = []
        for cb in compsB:
            c = list(vac); c[siteA] = ca; c[siteB] = cb
            row.append(Kpoly(tuple(c)))
        Kmat.append(row)
    return sp.Matrix(Kmat)
print('=== (a) F-only, sites 2,4: u_s = gamma_s e4 + delta_s e3 ===')
KF = form(1, [3, 2], 3, [3, 2])   # rows: site-2 (e4: gamma, e3: delta); cols: site-4
print(KF)
n2 = hS * th12 * th23 * (th24 - hS); n4 = hS * th14 * th24 * th34
rho2, rho4 = sp.symbols('rho2 rho4')
KFn = sp.simplify(KF.subs({r2: n2 * rho2, r4: n4 * rho4}))
print('in rho_s = r_s/n_s, n2 = h th12 th23 (th24-h), n4 = h th14 th24 th34:'); print(sp.factor(KFn))
g2, d2, g4, d4 = sp.symbols('gamma2 delta2 gamma4 delta4')
Kc = sp.expand((sp.Matrix([[g2, d2]]) * KFn * sp.Matrix([g4, d4]))[0])
coord = (g2 + d2 * rho2) * (g4 + d4 * rho4) - hS * g2 * d4 * rho2 / th24 - hS * d2 * d4 * rho2 * rho4 / (th24 + hS)
print("coordinator's formula (with rho2 in the gamma2 delta4 exchange term) exact:", sp.simplify(Kc - coord) == 0)
coord2 = (g2 + d2 * rho2) * (g4 + d4 * rho4) - hS * g2 * d4 * rho4 / th24 - hS * d2 * d4 * rho2 * rho4 / (th24 + hS)
print("Vandermonde reading (rho4 in that term) exact:", sp.simplify(Kc - coord2) == 0)
dF = sp.factor(sp.expand(KFn.det()))
print('det of the 2x2 form:', dF, ' -> separable f(rho2) g(rho4)?', sp.factor(dF).is_Mul and all(len(f.free_symbols) <= 1 for f in sp.Mul.make_args(sp.factor(dF))))
print('=== (a2) F-only k<=2: u_s in span(e2,e3,e4): generic rank of the 3x3 form (<=2 needed for any one-row-per-site determinant) ===')
KF3 = form(1, [3, 2, 1], 3, [3, 2, 1])
print('rank at random r:', KF3.subs({r2: sp.Rational(7, 3), r4: sp.Rational(-5, 11)}).rank())
print('=== (b) A-only, sites 1,3: v_s = alpha_s e1 + beta_s e4 (F sites e4) ===')
KA = form(0, [0, 3], 2, [0, 3])
print(KA)
dA = sp.factor(sp.expand(KA.det()))
print('det:', dA)
print('=== (b2) A-only v_s in span(e1,e2,e3) (BPS with e4): rank of 3x3 form at random tau ===')
KA3 = form(0, [0, 1, 2], 2, [0, 1, 2])
rnd = {A1[1]: sp.Rational(3, 7), A1[2]: sp.Rational(-2, 5), A1[3]: sp.Rational(5, 9), A3[1]: sp.Rational(-7, 4), A3[2]: sp.Rational(1, 6), A3[3]: sp.Rational(9, 5), r2: sp.Rational(2, 3), r4: sp.Rational(-3, 7)}
print('rank:', KA3.subs(rnd).rank())
print('=== (c) mixed A-F pairs (A: span(e1,e4), F: span(e4,e3)) ===')
for (sa, sf) in [(0, 1), (0, 3), (2, 1), (2, 3)]:
    K = form(sa, [0, 3], sf, [3, 2])
    d = sp.factor(sp.expand(K.det()))
    print(f' A site {sa+1}, F site {sf+1}: K =', K.tolist(), '\n    det =', d)
pickle.dump({'KF': KF, 'KA': KA}, open('bilinear_forms.pkl', 'wb'))
