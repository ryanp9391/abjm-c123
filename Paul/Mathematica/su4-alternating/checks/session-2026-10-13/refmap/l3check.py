"""l3check.py -- item 'verify at L = 3 pairs': for each passing reference (v(x)u)^3 (6 sites), exact mod p = 268435399:
 (i) the five words localise: word_w = lambda_w * (pair-1 word) (x) (v(x)u)^2;  (ii) sum_w (c_w/lambda_w) word_w = T_kappa/(1-kappa y)^2 (x) (v(x)u)^2
with c_w the pair-level coefficients (free ones set to 0).  Python fallback for notebook section 7.  Output: l3check_out.txt"""
import sys, os
from rm_base import *
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'session-2026-10-12', 'lpairs'))
import eng
from eng import modp, P as PR
def tw_u(w, z):
    out = []
    for j in range(4):
        p = Fr(1)
        for m in range(3): p *= (z[j] - w[m])
        out.append(p)
    return out
out = open('l3check_out.txt', 'w')
def Pr(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
y = Fr(5, 4); zz = [Fr(3, 2), y, y, Fr(4, 7)]
thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)]
def chain(L):
    th = []
    for n in range(L): th += [thF[n] + 1, thF[n]]
    return eng.Chain(L, th, h, zz, 'p')
def mp(vs): return np.array([modp(x) for x in vs], dtype=np.int64)
def words_chain(ch, psi):
    u1 = ch.th[1]; u2 = ch.th[1] + h
    x1 = ch.apply_t(1, u1, psi); x2 = ch.apply_t(2, u2, psi); x3 = ch.apply_t(3, u2, psi, reg=True); x31 = ch.apply_t(1, u1, x3)
    return [psi, x1, x2, x3, x31]
def lamw(ch):
    u1 = ch.th[1]; u2 = ch.th[1] + h
    l1 = eng.lam(ch, 1, u1, 1); l2 = eng.lam(ch, 2, u2, 1); l3 = eng.lam(ch, 3, u2, 1)
    return [Fr(1), l1, l2, l3, l3 * l1]
def check(L, name, v, u):
    v = vec(*v); u = vec(*u)
    ch = chain(L)
    psi = eng.prod_state([v, u] * L, 'p'); rest = eng.prod_state([v, u] * (L - 1), 'p')
    ws = words_chain(ch, psi)
    pw = words(v, u, pair_ops_z(zz))
    lam = lamw(ch)
    loc = []
    for w, p_, l in zip(ws, pw, lam):
        exp_ = (modp(l) * (np.multiply.outer(mp(p_), rest).reshape(-1) % PR)) % PR
        loc.append(bool(np.array_equal(w.reshape(-1) % PR, exp_)))
    res = []
    for kap in KAPS + [1 / zz[3]]:
        tgt = [t / (1 - kap * y) ** 2 for t in target(v, u, kap, zz)]
        c = solve_coeffs(pw, tgt)
        if c is None: res.append('no pair-level solution'); continue
        lhs = np.zeros(ws[0].size, dtype=np.int64)
        for ci, w, l in zip(c, ws, lam): lhs = (lhs + modp(ci / l) * (w.reshape(-1) % PR)) % PR
        rhs = (np.multiply.outer(mp(tgt), rest).reshape(-1)) % PR
        res.append(bool(np.array_equal(lhs, rhs)))
    Pr('L=%d' % L, name, '| localisation of 5 words:', loc, '| K-identity at kappa = 3/10, -7/5, 1/z4:', res)
refs = [
 ('W_a, a=1/3', (1, 0, 0, Fr(1, 3)), (Fr(-1, 3), 0, 0, 1)),
 ('<0_w| limit, generic w: (0,1,-1,0) (x) u_w', (0, 1, -1, 0), tw_u([Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9)], zz)),
 ('<0_w| limit, w1=z1: e1 (x) u_w', (1, 0, 0, 0), tw_u([zz[0], Fr(9, 8), Fr(3, 7), Fr(5, 9)], zz)),
 ('<0|comp limit: (0,1,-1,0) (x) (1,1,1,1)', (0, 1, -1, 0), (1, 1, 1, 1)),
 ('R2: v=(1,2,3,0), u=(0,3,-2,5)', (1, 2, 3, 0), (0, 3, -2, 5)),
 ('R5: v=e1, u=(0,2,3,5)', (1, 0, 0, 0), (0, 2, 3, 5)),
 ('R6: v=(0,2,3,0), u=(1,3,-2,4)', (0, 2, 3, 0), (1, 3, -2, 4)),
 ('closure R1: v=(1,2,3,2), u=(0,3,-2,0)', (1, 2, 3, 2), (0, 3, -2, 0)),
 ('negative control: random BPS (no K-family) -- localisation only', (1, 2, -1, 3), (Fr(5, 3), Fr(-1, 6), Fr(4, 3), 0)),
 ('negative control: Borel v=(1,2,3,0), u=(-5,1,1,3) [p.q=5]', (1, 2, 3, 0), (-5, 1, 1, 3)),
]
for L in (2, 3):
    for r in refs: check(L, *r)
out.close()
