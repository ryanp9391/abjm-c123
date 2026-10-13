"""step0a.py -- Step 0a of the plan, EXACT (Fractions, pgbase.Chain, L = 2 = A F A F, generic theta and z):
 <W_a^L| t_k^{G}(u) = tau_vac,k(u) <W_a^L|  (k = 1,2,3)  for G_a = Lambda + a (z1 - z4) E41  (and NOT for the other sign), and <W_a^L| is no eigenbra of t^Lambda for a != 0.
 Also: the E14 twist, the companion-frame images V G V^-1 and V^-1 G V, G_a = exp(a E41) Lambda exp(-a E41)."""
import sys, os
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'session-2026-10-12', 'lpairs')))
import numpy as np
from fractions import Fraction as Fr
import pgbase as pg
import sympy as sp
h = Fr(1, 2); th = [Fr(2, 7), Fr(-3, 5), Fr(8, 11), Fr(-11, 8)]; z = [Fr(3, 2), Fr(5, 4), Fr(2, 3), Fr(4, 7)]; a = Fr(1, 3)
sites = ['A', 'F', 'A', 'F']
def Wbra(al, Lp=2):
    vecs = []
    for _ in range(Lp): vecs += [[1, 0, 0, al], [-al, 0, 0, 1]]
    return pg.prod_state([[Fr(x) for x in v] for v in vecs])
def vacbra(Lp=2):
    return pg.prod_state([[Fr(x) for x in v] for v in [[1, 0, 0, 0], [0, 0, 0, 1]] * Lp])
def mkG(extra):
    G = pg.diagm(z).copy()
    for (i, j), v in extra.items(): G[i, j] = G[i, j] + v
    return G
def ratio_check(bra, ch, k, u):
    out = ch.apply_t(k, u, bra)
    nz = np.argwhere(bra != 0)[0]; idx = tuple(nz)
    lam = out[idx] / bra[idx]
    return lam, bool(np.all(out == lam * bra))
Wb = Wbra(a); vac = vacbra()
us = [Fr(1, 5), Fr(-2, 3), Fr(3, 7)]
res = {}
chL = pg.Chain(sites, th, h, pg.diagm(z))
tauvac = {(k, u): ratio_check(vac, chL, k, u)[0] for k in (1, 2, 3) for u in us}
print('tau_vac,k(u):', {k_: str(v) for k_, v in tauvac.items()})
cands = {
 'G_a+ = Lambda + a(z1-z4) E41': mkG({(3, 0): a * (z[0] - z[3])}),
 'G_a- = Lambda + a(z4-z1) E41': mkG({(3, 0): a * (z[3] - z[0])}),
 'Lambda (a != 0)': mkG({}),
 'Lambda + (z4-z1)/a E14': mkG({(0, 3): (z[3] - z[0]) / a}),
 'Lambda + (z1-z4)/a E14': mkG({(0, 3): (z[0] - z[3]) / a}),
}
for name, G in cands.items():
    ch = pg.Chain(sites, th, h, G)
    row = []
    for k in (1, 2, 3):
        for u in us:
            lam, ok = ratio_check(Wb, ch, k, u)
            row.append((k, u, ok, (lam == tauvac[(k, u)]) if ok else None))
    print('%-34s eigenbra for all (k,u)? %s ; eigenvalue = tau_vac? %s' % (name, all(r[2] for r in row), all(r[3] for r in row if r[2])) )
# vacuum control: <vac| eigenbra of t^Lambda and of t^{G_a+}
for name in ('G_a+ = Lambda + a(z1-z4) E41',):
    ch = pg.Chain(sites, th, h, cands[name]); print('control <vac| t^{G_a+}: eigenbra', all(ratio_check(vac, ch, k, u)[1] for k in (1, 2, 3) for u in us))
# the E14-twist eigenvalue (other vacuum LW = e4 x e1 per pair?)
ch = pg.Chain(sites, th, h, cands['Lambda + (z4-z1)/a E14'])
lams = [ratio_check(Wb, ch, 1, u) for u in us]
print('E14 twist: <W_a| eigenbra of t1?', [l[1] for l in lams], ' eigenvalue != tau_vac:', [l[0] != tauvac[(1, u)] for l, u in zip(lams, us)])
# group-theoretic form and companion frames
Ga = sp.Matrix(4, 4, lambda i, j: sp.Rational(str(cands['G_a+ = Lambda + a(z1-z4) E41'][i, j])))
Lam = sp.diag(*[sp.Rational(str(x)) for x in z]); asp = sp.Rational(1, 3)
E41 = sp.zeros(4, 4); E41[3, 0] = 1
g = (asp * E41).exp()
print('G_a+ == exp(a E41) Lambda exp(-a E41):', sp.simplify(g * Lam * g.inv() - Ga) == sp.zeros(4, 4))
V = sp.Matrix(4, 4, lambda i, j: sp.Rational(str(z[j])) ** (3 - i))     # columns (z_j^3, z_j^2, z_j, 1)
zs = [sp.Rational(str(x)) for x in z]
c1 = sum(zs); c2 = sum(zs[i] * zs[j] for i in range(4) for j in range(i + 1, 4)); c3 = sum(zs[i] * zs[j] * zs[k] for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4)); c4 = zs[0] * zs[1] * zs[2] * zs[3]
Gc = sp.Matrix([[c1, -c2, c3, -c4], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
print('V Lambda V^-1 == G_c(z):', V * Lam * V.inv() == Gc, ' ; V^-1 Lambda V == G_c?', V.inv() * Lam * V == Gc)
M1 = V * Ga * V.inv(); M2 = V.inv() * Ga * V
print('V G_a V^-1 =', M1.tolist()); print('   minus G_c =', (M1 - Gc).tolist(), ' rank', (M1 - Gc).rank())
print('V^-1 G_a V =', M2.tolist()); print('   minus V^-1 Lambda V =', (M2 - V.inv() * Lam * V).tolist(), ' rank', (M2 - V.inv() * Lam * V).rank())
# generic symbolic form of the companion image: V E41 V^-1 (rank one)
