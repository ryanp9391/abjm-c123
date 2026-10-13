"""hf.py L states -- first-order shifts A', B', C' of the normalised pair eigenvalues by Hellmann-Feynman in the untwisted eigenbasis,
compared with the Richardson derivative from the twisted scan (scan_L{L}_c0.3.pkl)."""
import pickle, numpy as np, sys
from ue_base import *
np.set_printoptions(linewidth=220, precision=5, suppress=True)
Lp = int(sys.argv[1]); states = [int(x) for x in sys.argv[2].split(',')]
d = pickle.load(open('scan_L%d_c0.3.pkl' % Lp, 'rb')); rows = d['rows']; K1 = d['K1']; phi = d['phi']; c = d['c']; thdir = d['thdir']
S = Setup(Lp)
def ops_at(eps):
    ch, z = S.chain(eps, phi, thdir); y = z[1]; kap = 1 - c * eps
    cc = S.Kcoef(ch, kap, y); s = 1 - kap * (z[0] + z[3]); uF = ch.th[1]; u0 = uF + h
    return ch, [(cc[1] / (kap * s), 1, uF, False), (cc[2] / (kap * s), 2, u0, False), (cc[3] / (kap * s), 3, u0, True)]
bysec = {}
for i in states: bysec.setdefault(tuple(S.W0[i]), []).append(i)
for w, ids in bysec.items():
    sec = S.SECT[w]
    hh = 1e-4
    chp, opsp = ops_at(hh); chm, opsm = ops_at(-hh); ch0, ops0 = ops_at(0.0)
    O = {}
    for name, k in (('A', 0), ('B', 1), ('C', 2)):
        Op = S.sector_matrix(chp, [opsp[k]], sec); Om = S.sector_matrix(chm, [opsm[k]], sec); O0 = S.sector_matrix(ch0, [ops0[k]], sec)
        O[name] = (O0, (Op - Om) / (2 * hh))
    M0 = S.sector_matrix(ch0, GEN, sec)
    ev, R = np.linalg.eig(M0); Rinv = np.linalg.inv(R)
    for i in ids:
        psi0 = S.PSI0[i][sec]; k0 = int(np.argmax(np.abs(Rinv @ psi0))); l0 = Rinv[k0] / (Rinv[k0] @ psi0)
        hf = {nm: l0 @ O[nm][1] @ psi0 for nm in 'ABC'}
        zero = {nm: l0 @ O[nm][0] @ psi0 for nm in 'ABC'}
        # Richardson from scan
        taus0 = K1[i][1]; A0 = None
        vals = {}
        for r in rows[i]:
            ch, z = S.chain(r['eps'], phi, thdir); y = z[1]; kap = 1 - c * r['eps']; cc = S.Kcoef(ch, kap, y); s = 1 - kap * (z[0] + z[3])
            t = r['taus']; vals[r['eps']] = np.array([cc[1] / (kap * s) * t[0], cc[2] / (kap * s) * t[1], cc[3] / (kap * s) * t[2]])
        ch0_, z0_ = S.chain(0.0, phi, thdir); cc0 = S.Kcoef(ch0_, 1.0, 1.0)
        v0 = np.array([-cc0[1] * taus0[0], -cc0[2] * taus0[1], -cc0[3] * taus0[2]])
        g = lambda e: (vals[e] - v0) / e
        rich = 2 * g(1e-4) - g(2e-4)
        print('#%d %s n=%s :  HF (A\',B\',C\') = %s | scan Richardson = %s | zeroth order (HF) %s vs scan %s' % (i, S.DYN[i], nval(Lp, S.DYN[i]), np.round([hf['A'], hf['B'], hf['C']], 5), np.round(rich, 5), np.round([zero[n] for n in 'ABC'], 4), np.round(v0, 4)))
