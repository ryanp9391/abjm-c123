import pickle, numpy as np, sys
from ue2 import *
np.set_printoptions(linewidth=220, precision=6, suppress=True)
Lp = int(sys.argv[1]); d = pickle.load(open('scan2_L%d.pkl' % Lp, 'rb')); rows = d['rows']; a = d['a']; c = d['c']; phi = d['phi']
S = Setup(Lp)
print('odd-bridge chain: <bra_eps|Psi_eps> vs word_eig * w (exact), limit (a-a_) C_1 ... L=%d' % Lp)
for i in sorted(rows):
    n = nval(Lp, S.DYN[i]); r = rows[i][-1]; eps = r['eps']
    C1 = S.C0[i][1]
    X11 = None
    for k in (1, 2):
        o = r['odd'][k]
        kap = 1 - c * eps
        ap = (a * (c + phi[3]) / (c + phi[0]) if k == 1 else a * (c + phi[0]) / (c + phi[3]))   # a'' (A family) or a' (F family)
        lim = (ap - a) if k == 1 else (a - ap)
        ex = o['lhs'] / o['nrm']
        print('#%d %s n=%s k=%d: |lhs-rhs|=%.1e  lhs/nrm=%s  (a-a\') C1=%s  wordeig=%s' % (i, S.DYN[i], n, k, abs(o['lhs'] - o['rhs']), np.round(ex, 6), np.round(lim * C1, 6), np.round(o['wordeig'], 6)))
