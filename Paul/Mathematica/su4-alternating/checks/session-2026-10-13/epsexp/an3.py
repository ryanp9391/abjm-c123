"""an3.py L -- normalised pair variables A = tauF/(h^2 l1), B = tau2/(2 h^4 y l2), C = tau3hat/(6 h^5 y^2 l3);  K = s^2 + kappa s (A+B+C) + kappa^2 A C / 4.
First-order data A', B', C' (Richardson from the twisted scan), K' decomposition, orders in eps."""
import pickle, numpy as np, sys, math
from ue_base import *
np.set_printoptions(linewidth=240, precision=5, suppress=True)
Lp = int(sys.argv[1]); d = pickle.load(open('scan_L%d_c0.3.pkl' % Lp, 'rb')); rows = d['rows']; K1 = d['K1']; phi = d['phi']; c = d['c']; thdir = d['thdir']
S = Setup(Lp)
def ABC(eps, taus):
    ch, z = S.chain(eps, phi, thdir); y = z[1]; kap = 1 - c * eps
    cc = S.Kcoef(ch, kap, y); s = 1 - kap * (z[0] + z[3])
    return np.array([cc[1] / (kap * s) * taus[0], cc[2] / (kap * s) * taus[2 - 1 + 0] if False else cc[2] / (kap * s) * taus[1], cc[3] / (kap * s) * taus[2]]), cc, s, kap
ch0, z0 = S.chain(0.0, phi, thdir)
print('L=%d : state dyn n | K1 | A0 B0 C0 | A\' B\' C\' (Richardson) | K\' (Richardson of K(eps)/... ) | K\'(chain rule from A\',B\',C\' + explicit) ' % Lp)
for i in sorted(rows):
    n = nval(Lp, S.DYN[i])
    if n is None or n > 3: continue
    taus0 = K1[i][1]
    A0, cc0, s0, k0 = ABC(0.0, taus0)
    ee = [r['eps'] for r in rows[i]]
    g = lambda e, v: v
    vals = {r['eps']: ABC(r['eps'], r['taus'])[0] for r in rows[i]}
    e1, e2 = 2e-4, 1e-4
    gA = lambda e: (vals[e] - A0) / e
    Ap = 2 * gA(e2) - gA(e1)
    Kv = {r['eps']: r['K'] for r in rows[i]}
    K0 = K1[i][0]
    gK = lambda e: (Kv[e] - K0) / e
    Kp = 2 * gK(e2) - gK(e1)
    # chain rule: K = s^2 + kap s (A+B+C) + kap^2 A C/4 ; at eps=0: kap=1, s=-1
    A_, B_, C_ = A0
    dK_dA = -1 + 0.25 * C_; dK_dB = -1; dK_dC = -1 + 0.25 * A_
    # explicit part: d/deps at fixed A,B,C : s = 1 - kap (z1+z4) ; kap = 1 - c eps ; z = 1 + eps phi
    def Kfix(e):
        kap = 1 - c * e; z1 = 1 + e * phi[0]; z4 = 1 + e * phi[3]; s = 1 - kap * (z1 + z4)
        return s * s + kap * s * (A_ + B_ + C_) + kap ** 2 * A_ * C_ / 4
    he = 1e-6
    Kexp = (Kfix(he) - Kfix(-he)) / (2 * he)
    Kchain = dK_dA * Ap[0] + dK_dB * Ap[1] + dK_dC * Ap[2] + Kexp
    print('#%2d %s n=%s | K1=%+.5f%+.5fi | A0,B0,C0 = %s | dA,dB,dC = %s | Kp(Rich)=%+.5f%+.5fi | chain=%+.5f%+.5fi | explicit part %+.5f' % (
        i, S.DYN[i], n, K0.real, K0.imag, np.round(A0, 5), np.round(Ap, 4), Kp.real, Kp.imag, Kchain.real, Kchain.imag, Kexp))
