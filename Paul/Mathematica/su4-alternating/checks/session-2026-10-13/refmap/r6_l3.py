"""r6_l3.py -- chain-level (L = 2, 3 pairs, exact mod p) check of the A-frozen operator K^F_kappa = [1 - k e1' + k^2 e2'] + 2k X1 + 8k^2 X2 (pair-level), i.e.
c1 + (2k/lam1) X1 + (8k^2/lam2) X2 on the chain, on (0,1,-1,0) (x) (u1,u2,u2,u4)."""
import sys, os
from helpers import *
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'session-2026-10-12', 'lpairs'))
import eng
from eng import modp, P as PR
y = Fr(5, 4); zz = [Fr(3, 2), y, y, Fr(4, 7)]
thF = [Fr(-3, 5), Fr(-11, 8), Fr(-7, 9)]
e1 = zz[0] + y + zz[3]; e2 = y * zz[0] + y * zz[3] + zz[0] * zz[3]
def mp(vs): return np.array([modp(x) for x in vs], dtype=np.int64)
for L in (2, 3):
    th = []
    for n in range(L): th += [thF[n] + 1, thF[n]]
    ch = eng.Chain(L, th, h, zz, 'p')
    u1 = th[1]; u2 = th[1] + h
    l1 = eng.lam(ch, 1, u1, 1); l2 = eng.lam(ch, 2, u2, 1)
    for nm, u in (('twisted vacuum limit', tw_u([Fr(7, 5), Fr(9, 8), Fr(3, 7), Fr(5, 9)], zz)), ('companion vacuum limit', [Fr(1)] * 4), ('generic u1,u2,u4', [Fr(2), Fr(-3, 7), Fr(-3, 7), Fr(5)])):
        v = vec(0, 1, -1, 0); u = list(u)
        psi = eng.prod_state([v, u] * L, 'p'); rest = eng.prod_state([v, u] * (L - 1), 'p')
        x1 = ch.apply_t(1, u1, psi).reshape(-1) % PR; x2 = ch.apply_t(2, u2, psi).reshape(-1) % PR
        res = []
        for kap in (Fr(3, 10), Fr(-7, 5), Fr(7, 4), Fr(2, 9)):
            c1 = 1 - kap * e1 + kap**2 * e2
            lhs = (modp(c1) * (psi.reshape(-1) % PR) + modp(2 * kap / l1) * x1 + modp(8 * kap**2 / l2) * x2) % PR
            tgt = [t / (1 - kap * y) ** 2 for t in target(v, u, kap, zz)]
            rhs = np.multiply.outer(mp(tgt), rest).reshape(-1) % PR
            res.append(bool(np.array_equal(lhs, rhs)))
        print('L=%d' % L, nm, ': c1 psi + (2k/lam1) psi.X1 + (8k^2/lam2) psi.X2 == T_kappa/(1-ky)^2 (x) rest  for kappa = 3/10, -7/5, 7/4, 2/9:', res)
